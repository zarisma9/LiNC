#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate conformers at defined phi/psi angles from an input SMILES string.
(After setting the dihedral angles, optimize the remaining geometry to reduce atomic clashes.)
Improvement: systematic atom naming.
"""
import sys
from collections import defaultdict, deque
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms
from openbabel import pybel

def assign_atom_names(mol):
    """
    Apply systematic atom naming rules:
    - Backbone: N, CA, C, O, H, HA
    - Acetyl (N-term cap): C101, O100, C102, HAC1/HAC2/HAC3
    - NME (C-term cap): N100, C103, HNM1/HNM2/HNM3/HNM4
    - Side chain: C1, C2, ..., N1, O1, ... in BFS order 
    - Side-chain H: H1, H2 for H atoms attached to C / HN1, HO1 for H atoms attached to N or O
    """
    assigned = set()
    
    # ── 1. Backbone ───────────────────────────────────────────────
    all_matches = mol.GetSubstructMatches(Chem.MolFromSmarts('[NH1][C@@H][C]=O'))
    if not all_matches:
        # Fallback: try without chirality
        all_matches = mol.GetSubstructMatches(Chem.MolFromSmarts('[NH1][CH][C]=O'))
    
    if not all_matches:
        raise ValueError("Backbone motif [NH1][C@@H][C]=O could not be found.")
    
    N_idx, CA_idx, C_idx, O_idx = all_matches[0]
    
    for idx, bname in [(N_idx,'N'),(CA_idx,'CA'),(C_idx,'C'),(O_idx,'O')]:
        mol.GetAtomWithIdx(idx).SetProp("name", bname)
        assigned.add(idx)
    
    # Backbone hydrogens
    for nb in mol.GetAtomWithIdx(N_idx).GetNeighbors():
        if nb.GetSymbol() == 'H' and nb.GetIdx() not in assigned:
            nb.SetProp("name", "H")
            assigned.add(nb.GetIdx())
    
    for nb in mol.GetAtomWithIdx(CA_idx).GetNeighbors():
        if nb.GetSymbol() == 'H' and nb.GetIdx() not in assigned:
            nb.SetProp("name", "HA")
            assigned.add(nb.GetIdx())
    
    print(f"  Backbone: N={N_idx}, CA={CA_idx}, C={C_idx}, O={O_idx}")
    
    # ── 2. Acetyl (N-term cap) ────────────────────────────────────
    acetyl_C = None
    for nb in mol.GetAtomWithIdx(N_idx).GetNeighbors():
        if nb.GetIdx() not in assigned and nb.GetSymbol() == 'C':
            acetyl_C = nb.GetIdx()
    
    if acetyl_C is not None:
        mol.GetAtomWithIdx(acetyl_C).SetProp("name", "C101")
        assigned.add(acetyl_C)
        
        for nb in mol.GetAtomWithIdx(acetyl_C).GetNeighbors():
            if nb.GetIdx() in assigned:
                continue
            if nb.GetSymbol() == 'O':
                nb.SetProp("name", "O100")
                assigned.add(nb.GetIdx())
            elif nb.GetSymbol() == 'C':
                nb.SetProp("name", "C102")
                assigned.add(nb.GetIdx())
                h_idx = 1
                for nb2 in nb.GetNeighbors():
                    if nb2.GetSymbol() == 'H' and nb2.GetIdx() not in assigned:
                        hname = f"HAC{h_idx}"
                        nb2.SetProp("name", hname)
                        assigned.add(nb2.GetIdx())
                        h_idx += 1
        print(f"  Acetyl cap assigned")
    
    # ── 3. NME (C-term cap) ───────────────────────────────────────
    nme_N = None
    for nb in mol.GetAtomWithIdx(C_idx).GetNeighbors():
        if nb.GetIdx() not in assigned and nb.GetSymbol() == 'N':
            nme_N = nb.GetIdx()
    
    if nme_N is not None:
        mol.GetAtomWithIdx(nme_N).SetProp("name", "N100")
        assigned.add(nme_N)
        
        for nb in mol.GetAtomWithIdx(nme_N).GetNeighbors():
            if nb.GetIdx() in assigned:
                continue
            if nb.GetSymbol() == 'H':
                nb.SetProp("name", "HNM4")
                assigned.add(nb.GetIdx())
            elif nb.GetSymbol() == 'C':
                nb.SetProp("name", "C103")
                assigned.add(nb.GetIdx())
                h_idx = 1
                for nb2 in nb.GetNeighbors():
                    if nb2.GetSymbol() == 'H' and nb2.GetIdx() not in assigned:
                        hname = f"HNM{h_idx}"
                        nb2.SetProp("name", hname)
                        assigned.add(nb2.GetIdx())
                        h_idx += 1
        print(f"  NME cap assigned")
    
    # ── 4. Side chain (heavy atoms, BFS) ──────────────────────────
    sc_counter = defaultdict(int)
    heavy_atom_names = {}
    
    queue = deque()
    for nb in mol.GetAtomWithIdx(CA_idx).GetNeighbors():
        if nb.GetIdx() not in assigned and nb.GetSymbol() != 'H':
            queue.append(nb.GetIdx())
    
    visited = set(assigned)
    
    while queue:
        idx = queue.popleft()
        if idx in visited:
            continue
        visited.add(idx)
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() == 'H':
            continue
        elem = atom.GetSymbol()
        sc_counter[elem] += 1
        sc_name = elem + str(sc_counter[elem])
        atom.SetProp("name", sc_name)
        assigned.add(idx)
        heavy_atom_names[idx] = sc_name
        for nb in atom.GetNeighbors():
            if nb.GetIdx() not in visited:
                queue.append(nb.GetIdx())
    
    # ── 5. Side-chain hydrogens ───────────────────────────────────
    for idx, hv_name in heavy_atom_names.items():
        elem = mol.GetAtomWithIdx(idx).GetSymbol()
        num_part = hv_name[len(elem):]  # C1 -> "1", C10 -> "10"
        
        if elem == 'C':
            prefix = f"H{num_part}"
        else:
            prefix = f"H{elem}{num_part}"
        
        h_neighbors = [nb for nb in mol.GetAtomWithIdx(idx).GetNeighbors()
                       if nb.GetSymbol() == 'H' and nb.GetIdx() not in assigned]
        
        for i, nb in enumerate(h_neighbors):
            if len(h_neighbors) == 1:
                hname = prefix
            else:
                hname = f"{prefix}{i+1}"
            
            nb.SetProp("name", hname)
            assigned.add(nb.GetIdx())
    
    # ── 6. Check for unassigned atoms ─────────────────────────────
    unassigned = [a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() not in assigned]
    if unassigned:
        print(f"  WARNING: unassigned atoms: {unassigned}")
    
    print(f"  Total assigned atoms: {len(assigned)}/{mol.GetNumAtoms()}")

def find_phi_psi_indices(mol):
    """
    Find the atoms defining the peptide-backbone phi/psi dihedral angles.
    """
    carbonyl_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            has_double_O = False
            for bond in atom.GetBonds():
                other = bond.GetOtherAtom(atom)
                if other.GetSymbol() == 'O' and bond.GetBondType() == Chem.BondType.DOUBLE:
                    has_double_O = True
                    break
            if has_double_O:
                carbonyl_carbons.append(atom)
    
    print(f"Found {len(carbonyl_carbons)} carbonyl carbons")
    
    for c_carbonyl in carbonyl_carbons:
        n_atoms = [nbr for nbr in c_carbonyl.GetNeighbors() if nbr.GetSymbol() == 'N']
        
        if not n_atoms:
            continue
            
        for n_i in n_atoms:
            c_alpha_candidates = []
            for nbr in n_i.GetNeighbors():
                if nbr.GetSymbol() == 'C' and nbr.GetIdx() != c_carbonyl.GetIdx():
                    c_alpha_candidates.append(nbr)
            
            if not c_alpha_candidates:
                continue
            
            for c_alpha in c_alpha_candidates:
                c_i_carbonyl = None
                for nbr in c_alpha.GetNeighbors():
                    if nbr.GetSymbol() == 'C' and nbr.GetIdx() != n_i.GetIdx():
                        for bond in nbr.GetBonds():
                            other = bond.GetOtherAtom(nbr)
                            if other.GetSymbol() == 'O' and bond.GetBondType() == Chem.BondType.DOUBLE:
                                c_i_carbonyl = nbr
                                break
                    if c_i_carbonyl:
                        break
                
                if not c_i_carbonyl:
                    continue
                
                n_ip1 = None
                for nbr in c_i_carbonyl.GetNeighbors():
                    if nbr.GetSymbol() == 'N' and nbr.GetIdx() != c_alpha.GetIdx():
                        n_ip1 = nbr
                        break
                
                if not n_ip1:
                    continue
                
                c_im1 = c_carbonyl
                
                phi_idx = (c_im1.GetIdx()+1, n_i.GetIdx()+1, c_alpha.GetIdx()+1, c_i_carbonyl.GetIdx()+1)
                psi_idx = (n_i.GetIdx()+1, c_alpha.GetIdx()+1, c_i_carbonyl.GetIdx()+1, n_ip1.GetIdx()+1)
                
                print(f"Found backbone atoms:")
                print(f"  C(i-1): {c_im1.GetIdx()+1}")
                print(f"  N(i):   {n_i.GetIdx()+1}")
                print(f"  CA(i):  {c_alpha.GetIdx()+1}")
                print(f"  C(i):   {c_i_carbonyl.GetIdx()+1}")
                print(f"  N(i+1): {n_ip1.GetIdx()+1}")
                
                return phi_idx, psi_idx
    
    raise ValueError("Atoms defining the phi/psi dihedral angles could not be found.")

def write_mol2_with_names(mol, filename):
    """
    Write an RDKit molecule to MOL2 while
    using mol.GetAtomWithIdx(i).GetProp("name") as the atom name.
    """
    from rdkit.Chem import rdmolfiles
    
    # Create a temporary MOL2 file
    temp_file = filename + ".tmp"
    try:
        rdmolfiles.MolToMol2File(mol, temp_file)
    except Exception:
        # Use OpenBabel if RDKit fails
        block = Chem.MolToMolBlock(mol)
        obmol = pybel.readstring("mol", block)
        obmol.write("mol2", temp_file, overwrite=True)
    
    # Read the MOL2 file and replace atom names
    lines = open(temp_file).readlines()
    out_lines = []
    in_atoms = False
    atom_idx = 0
    
    for line in lines:
        if line.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            out_lines.append(line)
            continue
        if line.startswith("@<TRIPOS>BOND"):
            in_atoms = False
            out_lines.append(line)
            continue
        
        if in_atoms and line.strip():
            parts = line.split()
            atom = mol.GetAtomWithIdx(atom_idx)
            
            # Use the name stored in GetProp("name")
            try:
                atom_name = atom.GetProp("name")
            except:
                # Use a default name if no atom name is available
                atom_name = f"{atom.GetSymbol()}{atom_idx+1}"
            
            # Reconstruct the MOL2 atom line
            new_line = (
                f"{int(parts[0]):>7} {atom_name:<6} {float(parts[2]):>8.4f} "
                f"{float(parts[3]):>8.4f} {float(parts[4]):>8.4f} "
                f"{parts[5]:<8} {parts[6]:>3} {parts[7]:<8} {parts[-1]:>10}\n"
            )
            out_lines.append(new_line)
            atom_idx += 1
        else:
            out_lines.append(line)
    
    # Write the final file
    open(filename, 'w').writelines(out_lines)
    
    # Remove the temporary file
    import os
    os.remove(temp_file)

def report_backbone_atoms(mol2_file, phi_idx, psi_idx):
    """Report backbone atom information required by the bash workflow."""
    ids_to_names = {}
    with open(mol2_file) as f:
        in_atoms = False
        for line in f:
            if line.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                continue
            if line.startswith("@<TRIPOS>BOND"):
                break
            if in_atoms and line.strip():
                parts = line.split()
                atom_id = int(parts[0])
                atom_name = parts[1]
                ids_to_names[atom_id] = atom_name
    
    phi_names = [ids_to_names[i] for i in phi_idx]
    psi_names = [ids_to_names[i] for i in psi_idx]
    
    print(f"Backbone phi atoms: ids {phi_idx}, names {phi_names}")
    print(f"Backbone psi atoms: ids {psi_idx}, names {psi_names}")

def get_dihedral_angle(mol, atom_indices):
    """Measure the current dihedral angle in the molecule."""
    conf = mol.GetConformer()
    idx = [i-1 for i in atom_indices]
    angle = rdMolTransforms.GetDihedralDeg(conf, *idx)
    return angle

def set_dihedral_with_verification(mol, atom_indices, target_angle, max_iterations=10):
    """
    Replace the previous simple angle-setting procedure.
    Set and restrain the target angle while relaxing the structure to reduce atomic clashes.
    """
    conf = mol.GetConformer()
    idx = [i-1 for i in atom_indices]
    
    # 1. Set the angle
    rdMolTransforms.SetDihedralDeg(conf, *idx, target_angle)
    
    # 2. Set up the force field
    mp = AllChem.MMFFGetMoleculeProperties(mol)
    ff = AllChem.MMFFGetMoleculeForceField(mol, mp)
    
    if ff is None:
        ff = AllChem.UFFGetMoleculeForceField(mol)
    
    # 3. Torsion Constraint
    force_constant = 1.0e5
    ff.MMFFAddTorsionConstraint(idx[0], idx[1], idx[2], idx[3], False, target_angle-0.1, target_angle+0.1, force_constant)
    
    # 4. Minimize
    ff.Minimize(maxIts=1000)
    
    # 5. Report the relaxed angle
    current_angle = rdMolTransforms.GetDihedralDeg(conf, *idx)
    print(f"  Relaxed structure with fixed angle: {current_angle:.2f} deg")
    return True

def main():
    if len(sys.argv) != 3:
        print("Usage: python conformer.py \"SMILES\" RESNAME", file=sys.stderr)
        sys.exit(1)
    
    smiles, resname = sys.argv[1], sys.argv[2]
    
    print(f"Input SMILES: {smiles}")
    print(f"Residue name: {resname}")
    print()
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Error: invalid SMILES.", file=sys.stderr)
        sys.exit(1)
    
    # Add hydrogens
    mol = Chem.AddHs(mol)
    print(f"Total atoms (with H): {mol.GetNumAtoms()}")
    
    # Apply atom naming before generating the 3D structure
    print("\nAssigning atom names...")
    assign_atom_names(mol)
    
    # Generate a 3D structure
    print("\nGenerating 3D structure...")
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xf00d
    result = AllChem.EmbedMolecule(mol, params)
    
    if result == -1:
        print("Warning: ETKDG failed. Retrying with useRandomCoords...", file=sys.stderr)
        result = AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
        if result == -1:
            print("Error: 3D structure generation failed.", file=sys.stderr)
            sys.exit(1)
    
    print("3D structure generated successfully.")
    
    # Initial whole-structure optimization
    print("\nOptimizing structure...")
    opt_result = AllChem.MMFFOptimizeMolecule(mol)
    if opt_result == -1:
        print("Warning: MMFF optimization failed. Retrying with UFF...", file=sys.stderr)
        AllChem.UFFOptimizeMolecule(mol)
    
    # Find phi/psi atom indices
    print("\nSearching for phi/psi atoms...")
    phi_idx, psi_idx = find_phi_psi_indices(mol)
    
    # Check initial angles
    print("\nInitial dihedral angles:")
    initial_phi = get_dihedral_angle(mol, phi_idx)
    initial_psi = get_dihedral_angle(mol, psi_idx)
    print(f"  phi = {initial_phi:.2f} deg")
    print(f"  psi = {initial_psi:.2f} deg")
    
    # Conformer-generation loop
    targets = {
        "alpha": (-60.0, -45.0),
        "beta": (-120.0, 120.0)
    }
    
    for name, (phi_target, psi_target) in targets.items():
        print(f"\n{'='*60}")
        print(f"Generating {name.upper()} conformer (phi={phi_target} deg, psi={psi_target} deg)")
        print('='*60)
        
        # Copy the molecule, including atom-name properties
        m2 = Chem.Mol(mol)
        
        # Set target angles and relax the structure
        print(f"Setting phi to {phi_target} deg and relaxing structure...")
        set_dihedral_with_verification(m2, phi_idx, phi_target)
        
        print(f"Setting psi to {psi_target} deg and relaxing structure...")
        set_dihedral_with_verification(m2, psi_idx, psi_target)
        
        # Check final angles
        final_phi = get_dihedral_angle(m2, phi_idx)
        final_psi = get_dihedral_angle(m2, psi_idx)
        
        print(f"\nFinal dihedral angles:")
        print(f"  phi = {final_phi:.2f} deg (target: {phi_target:.2f} deg)")
        print(f"  psi = {final_psi:.2f} deg (target: {psi_target:.2f} deg)")
        
        # Calculate angle errors
        phi_error = abs(final_phi - phi_target)
        psi_error = abs(final_psi - psi_target)
        if phi_error > 180: phi_error = 360 - phi_error
        if psi_error > 180: psi_error = 360 - psi_error
        
        print(f"  phi error: {phi_error:.2f} deg")
        print(f"  psi error: {psi_error:.2f} deg")
        
        # Save the MOL2 file with systematic atom names
        fn = f"{name}_conf.mol2"
        write_mol2_with_names(m2, fn)
        
        # Replace the residue name in the MOL2 file
        lines = open(fn).readlines()
        out_lines = []
        in_atoms = False
        for line in lines:
            if line.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                out_lines.append(line)
                continue
            if line.startswith("@<TRIPOS>BOND"):
                in_atoms = False
                out_lines.append(line)
                continue
            if in_atoms and line.strip():
                parts = line.split()
                new_line = (
                    f"{int(parts[0]):>7} {parts[1]:<6} {float(parts[2]):>8.4f} "
                    f"{float(parts[3]):>8.4f} {float(parts[4]):>8.4f} "
                    f"{parts[5]:<8} {parts[6]:>3} {resname:<8} {parts[-1]:>10}\n"
                )
                out_lines.append(new_line)
            else:
                out_lines.append(line)
        open(fn, 'w').writelines(out_lines)
        
        # Report backbone atom information
        report_backbone_atoms(fn, phi_idx, psi_idx)
        
        print(f"Created {fn}")
    
    print("\n" + "="*60)
    print("All tasks completed successfully.")
    print("="*60)

if __name__ == '__main__':
    main()

