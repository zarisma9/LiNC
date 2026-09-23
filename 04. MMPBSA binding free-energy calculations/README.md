# MMPBSA Analysis Workflow

This directory contains a streamlined post-processing and binding-energy analysis workflow for AMBER molecular dynamics trajectories.

The workflow performs:

1. trajectory centering with `cpptraj`
2. receptor, peptide, and lipid-anchor RMSD analysis
3. peptide-receptor minimum-distance analysis
4. MMPBSA topology generation
5. MM/GBSA and MM/PBSA calculations
6. per-residue energy decomposition

The main script is:

```text
run_MMPBSA.sh
```

---

## Workflow Overview

```text
Production trajectory
07_Prod_HMR.nc
        |
        v
cpptraj autoimage
        |
        +--> SYSTEM_ID_centered.nc
        |
        v
cpptraj structural analysis
        |
        +--> SYSTEM_ID_rmsd_bb.dat
        +--> SYSTEM_ID_rmsd_peptide.dat
        +--> SYSTEM_ID_rmsd_lipid.dat
        +--> SYSTEM_ID_mindist_ha.dat
        |
        v
ante-MMPBSA.py
        |
        +--> SYSTEM_ID_HMR_dry_complex.prmtop
        +--> SYSTEM_ID_HMR_receptor.prmtop
        +--> SYSTEM_ID_HMR_ligand.prmtop
        |
        v
MMPBSA.py.MPI
        |
        +--> FINAL_RESULTS_MMPBSA.dat
        +--> FINAL_DECOMP_MMPBSA.dat
```

All analysis files are written to:

```text
BASE_DIR/BA/
```

---

# Required Input Files

The script expects the following files inside `BASE_DIR`:

```text
SYSTEM_ID.hmass.prmtop
07_Prod_HMR.nc
```

For example:

```text
/home/guest_1/LDM/data/260304_lipid_MD/MEK_seed1/
├── MEK_seed1.hmass.prmtop
└── 07_Prod_HMR.nc
```

---

# Usage

Make the script executable:

```bash
chmod +x run_MMPBSA.sh
```

Run:

```bash
./run_MMPBSA.sh SYSTEM_ID BASE_DIR
```

Example:

```bash
./run_MMPBSA.sh \
MEK_seed1 \
/home/guest_1/LDM/data/260304_lipid_MD/MEK_seed1
```

The same parameters can also be supplied through environment variables:

```bash
SYSTEM_ID=MEK_seed1 \
BASE_DIR=/home/guest_1/LDM/data/260304_lipid_MD/MEK_seed1 \
./run_MMPBSA.sh
```

---

# Default Analysis Settings

The current default masks are:

```bash
REC_MASK="32-616"
PEP_MASK="1-31"
LIGAND_MASK=":1-31"
```

Interpretation:

```text
Residues 32-616 : receptor
Residues 1-31   : peptide / ligand
Residue 20      : lipid-anchor residue used for RMSD analysis
```

The default MMPBSA frame range is:

```bash
START_FRAME=1500
END_FRAME=2000
INTERVAL=5
```

The default number of MPI processes is:

```bash
MMPBSA_CORES=8
```

---

# Step 1. Center the Trajectory

The script generates:

```text
BA/center.in
```

and runs:

```text
cpptraj
```

using:

```text
autoimage anchor :32-616
```

The output trajectory is:

```text
BA/SYSTEM_ID_centered.nc
```

Example:

```text
BA/MEK_seed1_centered.nc
```

This centered trajectory is used for all subsequent structural and MMPBSA analyses.

---

# Step 2. Structural Analysis

The script generates:

```text
BA/rmsd.in
```

and calculates four structural metrics.

## Receptor Backbone RMSD

```text
RMSDbb
```

Selection:

```text
:32-616@N,CA,C
```

Output:

```text
SYSTEM_ID_rmsd_bb.dat
```

---

## Peptide Backbone RMSD

```text
RMSDpeptide
```

Selection:

```text
:1-31@N,CA,C
```

The `nofit` option is used after receptor alignment.

Output:

```text
SYSTEM_ID_rmsd_peptide.dat
```

---

## Lipid-Anchor RMSD

Residue 20 is analyzed independently:

```text
:20
```

Output:

```text
SYSTEM_ID_rmsd_lipid.dat
```

The current script uses all atoms of residue 20.

---

## Peptide-Receptor Minimum Distance

Heavy atoms only are used.

Peptide:

```text
:1-31&!@H=
```

Receptor:

```text
:32-616&!@H=
```

Output:

```text
SYSTEM_ID_mindist_ha.dat
```

---

# Step 3. MMPBSA Topology Generation

The script uses:

```text
ante-MMPBSA.py
```

to generate separate topology files for the complex, receptor, and ligand.

Generated files:

```text
SYSTEM_ID_HMR_dry_complex.prmtop
SYSTEM_ID_HMR_receptor.prmtop
SYSTEM_ID_HMR_ligand.prmtop
```

The default ligand mask is:

```text
:1-31
```

The default strip mask is:

```text
:WAT,HOH,TIP3,Na+,Cl-,K+,Mg2+,Zn2+
```

The script also uses:

```text
--radii mbondi3
```

---

# Step 4. MMPBSA Input Generation

The script automatically generates:

```text
BA/mmpbsa.in
```

Current settings are:

```text
&general
  startframe=1500,
  endframe=2000,
  interval=5,
  verbose=1,
  keep_files=0,
  strip_mask=':WAT,HOH,TIP3,Na+,Cl-,K+,Mg2+,Zn2+'
/
```

## GB Settings

```text
&gb
  igb=8,
  saltcon=0.150,
  surften=0.0072,
  surfoff=0.0
/
```

## PB Settings

```text
&pb
  indi=2.0,
  exdi=80.0,
  istrng=0.150,
  inp=2,
  radiopt=0
/
```

## Energy Decomposition

```text
&decomp
  idecomp=1,
  print_res='1-616',
  dec_verbose=1,
  csv_format=1
/
```

This performs per-residue decomposition over residues:

```text
1-616
```

---

# Step 5. Run MMPBSA

The script runs:

```bash
mpirun -np 8 MMPBSA.py.MPI
```

using the centered trajectory.

Main outputs:

```text
FINAL_RESULTS_MMPBSA.dat
FINAL_DECOMP_MMPBSA.dat
```

Both files are stored inside:

```text
BASE_DIR/BA/
```

---

# Changing Analysis Parameters

Most important settings can be overridden without editing the script.

## Change Receptor and Peptide Masks

Example:

```bash
REC_MASK="16-598" \
PEP_MASK="1-15" \
LIGAND_MASK=":1-15" \
./run_MMPBSA.sh SYSTEM_seed1 /path/to/SYSTEM_seed1
```

---

## Change Frame Range

Example:

```bash
START_FRAME=1000 \
END_FRAME=3000 \
INTERVAL=10 \
./run_MMPBSA.sh MEK_seed1 /path/to/MEK_seed1
```

---

## Change Number of MPI Processes

Example:

```bash
MMPBSA_CORES=16 \
./run_MMPBSA.sh MEK_seed1 /path/to/MEK_seed1
```

---

## Change Solvent/Ion Strip Mask

Example:

```bash
STRIP_MASK=':WAT,Na+,Cl-,K+' \
./run_MMPBSA.sh MEK_seed1 /path/to/MEK_seed1
```

---

# Expected Directory Structure

Before analysis:

```text
MEK_seed1/
├── MEK_seed1.hmass.prmtop
├── 07_Prod_HMR.nc
└── run_MMPBSA.sh
```

After analysis:

```text
MEK_seed1/
├── MEK_seed1.hmass.prmtop
├── 07_Prod_HMR.nc
├── run_MMPBSA.sh
└── BA/
    ├── center.in
    ├── rmsd.in
    ├── mmpbsa.in
    ├── MEK_seed1_centered.nc
    ├── MEK_seed1_rmsd_bb.dat
    ├── MEK_seed1_rmsd_peptide.dat
    ├── MEK_seed1_rmsd_lipid.dat
    ├── MEK_seed1_mindist_ha.dat
    ├── MEK_seed1_HMR_dry_complex.prmtop
    ├── MEK_seed1_HMR_receptor.prmtop
    ├── MEK_seed1_HMR_ligand.prmtop
    ├── FINAL_RESULTS_MMPBSA.dat
    └── FINAL_DECOMP_MMPBSA.dat
```

---

# Requirements

The following commands must be available in the environment:

```text
cpptraj
ante-MMPBSA.py
MMPBSA.py.MPI
mpirun
```

These are typically provided through an AMBER / AmberTools environment.

The script checks that each command is available before starting.

---

# Important Checks Before Running

Before launching the analysis, confirm the following:

```text
[ ] SYSTEM_ID is correct
[ ] SYSTEM_ID.hmass.prmtop exists
[ ] 07_Prod_HMR.nc exists
[ ] Receptor residue range is correct
[ ] Peptide residue range is correct
[ ] LIGAND_MASK matches the intended MMPBSA ligand
[ ] Residue 20 is the intended lipid-anchor residue
[ ] START_FRAME and END_FRAME fall within the available trajectory
[ ] Solvent and ion residue names match the topology
[ ] print_res='1-616' matches the intended decomposition range
```

---

# Critical Mask Check

The residue numbering used in:

```text
REC_MASK
PEP_MASK
LIGAND_MASK
print_res
```

must match the residue numbering in the AMBER topology.

For the current setup:

```text
Peptide  : residues 1-31
Receptor : residues 32-616
Anchor   : residue 20
```

If a different system has a different peptide length or receptor numbering, update the masks before running the analysis.

Incorrect masks can produce misleading RMSD, distance, receptor/ligand topology splitting, and MMPBSA decomposition results.

---

# Example Complete Command

```bash
chmod +x run_MMPBSA.sh

./run_MMPBSA.sh \
MEK_seed1 \
/home/guest_1/LDM/data/260304_lipid_MD/MEK_seed1
```

For a system with different residue ranges:

```bash
REC_MASK="16-598" \
PEP_MASK="1-15" \
LIGAND_MASK=":1-15" \
START_FRAME=1500 \
END_FRAME=2000 \
INTERVAL=5 \
MMPBSA_CORES=8 \
./run_MMPBSA.sh \
SYSTEM_seed1 \
/path/to/SYSTEM_seed1
```

---

# Main Output Files

For binding-energy analysis:

```text
FINAL_RESULTS_MMPBSA.dat
```

For residue-level energy contributions:

```text
FINAL_DECOMP_MMPBSA.dat
```

For structural stability analysis:

```text
SYSTEM_ID_rmsd_bb.dat
SYSTEM_ID_rmsd_peptide.dat
SYSTEM_ID_rmsd_lipid.dat
SYSTEM_ID_mindist_ha.dat
```
