# AMBER MD Preparation and Simulation Workflow

This directory contains a stepwise workflow for converting a Chai-1 prediction into an AMBER-ready molecular dynamics system and running minimization, equilibration, and production MD.

The workflow is:

```text
make_amber_pdb.sh
        |
        v
make_tleap_input.sh
        |
        v
make_hmass_prmtop.sh
        |
        v
make_MD_input.sh
        |
        v
run_MD.sh
```

---

## 1. `make_amber_pdb.sh`

Converts a Chai-1 CIF prediction to PDB, merges the ligand with the target protein residue using `modify_pdb.py`, and generates an AMBER-formatted PDB using `pdb4amber`.

General usage:

```bash
./make_amber_pdb.sh \
INPUT_CIF \
CONVERTED_PDB \
MODIFIED_PDB \
CHAIN \
RESNUM \
RESNAME \
FINAL_PDB
```

Example:

```bash
./make_amber_pdb.sh \
pred.model_idx_0.cif \
D4K_seed1.pdb \
D4K_final.pdb \
A \
20 \
D4K \
D4K_amber.pdb
```

This runs:

```text
CIF
 |
 v
Open Babel
 |
 v
D4K_seed1.pdb
 |
 v
modify_pdb.py
 |
 v
D4K_final.pdb
 |
 v
pdb4amber
 |
 v
D4K_amber.pdb
```

### CRITICAL CHECK BEFORE CONTINUING

**The noncanonical residue name in the generated `*_amber.pdb` file must match the residue name used in the AMBER residue parameter files.**

For example, if the parameter files are:

```text
D4K.prepi
D4K.frcmod
```

then the modified residue in:

```text
D4K_amber.pdb
```

must also be named:

```text
D4K
```

Do not continue to `tleap` until this has been verified.

A quick check can be performed with:

```bash
grep " D4K " D4K_amber.pdb
```

or by inspecting residue 20 directly:

```bash
grep " A  20 " D4K_amber.pdb
```

The target residue should appear with the intended noncanonical residue name.

Incorrect residue naming can prevent `tleap` from matching the residue to the loaded `.prepi` template and may result in missing atom types, missing parameters, incorrect residue assignment, or topology-generation failure.

Also confirm that the residue name is consistent across:

```text
*_amber.pdb
*.prepi
*.frcmod
tleap input
```

For example:

```text
D4K_amber.pdb
D4K.prepi
D4K.frcmod
```

should all refer to the same noncanonical residue definition.

---

## 2. `make_tleap_input.sh`

Generates a `tleap` input file for system construction, solvation, ion addition, and topology generation.

General usage:

```bash
./make_tleap_input.sh \
AMBERPREP \
AMBERPARAMS \
PDB_FILE \
K_COUNT \
CL_COUNT \
PRMTOP_NAME \
INPCRD_NAME
```

Example:

```bash
./make_tleap_input.sh \
D4K.prepi \
D4K.frcmod \
D4K_amber.pdb \
88 \
71 \
D4K_seed1.prmtop \
D4K_seed1.inpcrd
```

The generated `tleap.in` contains:

```text
source leaprc.protein.ff19SB
source leaprc.water.opc
source leaprc.mimetic.ff15ipq

loadamberprep D4K.prepi
loadamberparams D4K.frcmod

m = loadpdb D4K_amber.pdb

solvateOct m OPCBOX 12.0

addIonsRand m K+ 88 Cl- 71

saveAmberParm m D4K_seed1.prmtop D4K_seed1.inpcrd

quit
```

Run:

```bash
tleap -f tleap.in
```

Expected outputs:

```text
D4K_seed1.prmtop
D4K_seed1.inpcrd
```

Always inspect the `tleap` output for warnings about:

```text
Unknown residue
Missing atom type
Missing parameter
Unrecognized atom
Unmatched residue
```

Do not proceed if the noncanonical residue is not recognized correctly.

---

## 3. `make_hmass_prmtop.sh`

Generates a hydrogen-mass-repartitioned topology for 4 fs production MD.

General usage:

```bash
./make_hmass_prmtop.sh INPUT_PRMTOP OUTPUT_PRMTOP
```

Example:

```bash
./make_hmass_prmtop.sh \
D4K_seed1.prmtop \
D4K_seed1.hmass.prmtop
```

The script creates:

```text
hmass_parmed.in
```

with:

```text
parm D4K_seed1.prmtop
HMassRepartition
outparm D4K_seed1.hmass.prmtop
quit
```

and automatically runs:

```bash
parmed -i hmass_parmed.in
```

Expected output:

```text
D4K_seed1.hmass.prmtop
```

The standard topology is used for minimization, heating, and equilibration.

The HMR topology is used for the 4 fs production run.

---

## 4. `make_MD_input.sh`

Generates the seven AMBER input files used for minimization, heating, equilibration, and production.

Run:

```bash
./make_MD_input.sh
```

Generated files:

```text
01_Min1.in
02_Min2.in
03_Heat_NVT.in
04_Equil_NPT_1.in
05_Equil_NPT_2.in
06_Equil_NPT_3.in
07_Prod_HMR_4fs.in
```

Current protocol:

```text
Step 1  Minimization 1
        Solvent minimization
        Strong solute restraint

Step 2  Minimization 2
        Whole-system minimization
        Weak solute restraint

Step 3  Heating
        NVT
        0 -> 310 K

Step 4  Equilibration 1
        NPT
        Restraint weight = 2.0

Step 5  Equilibration 2
        NPT
        Restraint weight = 0.5

Step 6  Equilibration 3
        NPT
        Unrestrained

Step 7  Production
        NPT
        HMR topology
        4 fs timestep
        1 us target
```

The current restraint mask is:

```text
!:WAT,Na+,K+,Cl-
```

Check this mask before running a system with different solvent, ion, ligand, or residue naming.

---

## 5. `run_MD.sh`

Runs the complete AMBER simulation workflow.

Usage:

```bash
./run_MD.sh SYS_NAME
```

Example:

```bash
./run_MD.sh D4K_seed1
```

This automatically uses:

```text
D4K_seed1.prmtop
D4K_seed1.hmass.prmtop
D4K_seed1.inpcrd
```

The workflow is:

```text
01_Min1
   |
   v
02_Min2
   |
   v
03_Heat_NVT
   |
   v
04_Equil_NPT_1
   |
   v
05_Equil_NPT_2
   |
   v
06_Equil_NPT_3
   |
   v
07_Prod_HMR_4fs
```

Default execution commands are:

```text
CPU:
mpirun -np 4 pmemd.MPI

GPU:
pmemd.cuda
```

The number of MPI processes can be changed with:

```bash
MPI_NP=8 ./run_MD.sh D4K_seed1
```

---

# Complete Example

For a D4K system:

```bash
./make_amber_pdb.sh \
pred.model_idx_0.cif \
D4K_seed1.pdb \
D4K_final.pdb \
A \
20 \
D4K \
D4K_amber.pdb
```

### Verify the noncanonical residue name

```bash
grep " A  20 " D4K_amber.pdb
```

Confirm that the residue is named:

```text
D4K
```

Then generate the `tleap` input:

```bash
./make_tleap_input.sh \
D4K.prepi \
D4K.frcmod \
D4K_amber.pdb \
88 \
71 \
D4K_seed1.prmtop \
D4K_seed1.inpcrd
```

Run `tleap`:

```bash
tleap -f tleap.in
```

Generate the HMR topology:

```bash
./make_hmass_prmtop.sh \
D4K_seed1.prmtop \
D4K_seed1.hmass.prmtop
```

Generate AMBER MD input files:

```bash
./make_MD_input.sh
```

Run MD:

```bash
./run_MD.sh D4K_seed1
```

---

# Recommended Directory Structure

```text
MD/
├── README.md
├── modify_pdb.py
├── make_amber_pdb.sh
├── make_tleap_input.sh
├── make_hmass_prmtop.sh
├── make_MD_input.sh
├── run_MD.sh
│
├── D4K.prepi
├── D4K.frcmod
├── pred.model_idx_0.cif
│
└── D4K_seed1/
    ├── D4K_amber.pdb
    ├── D4K_seed1.prmtop
    ├── D4K_seed1.hmass.prmtop
    ├── D4K_seed1.inpcrd
    ├── 01_Min1.in
    ├── 02_Min2.in
    ├── 03_Heat_NVT.in
    ├── 04_Equil_NPT_1.in
    ├── 05_Equil_NPT_2.in
    ├── 06_Equil_NPT_3.in
    └── 07_Prod_HMR_4fs.in
```

---

# Final Checklist

Before starting production MD, confirm all of the following:

```text
[ ] Chai-1 prediction was converted successfully
[ ] Target residue and ligand were merged correctly
[ ] Noncanonical residue name in *_amber.pdb is correct
[ ] Residue name matches the .prepi residue definition
[ ] Correct .frcmod file is loaded
[ ] tleap recognizes the noncanonical residue
[ ] No missing atom types or missing parameters remain
[ ] Standard .prmtop and .inpcrd were generated successfully
[ ] HMR .prmtop was generated successfully
[ ] All seven MD input files exist
[ ] Restraint mask is appropriate for the current system
```

## Most Important Warning

> **Always verify the noncanonical residue name in the final `*_amber.pdb` before running `tleap`.**
>
> The residue name must match the residue definition used by the corresponding AMBER parameter files. A naming mismatch can cause the custom residue to be interpreted incorrectly or prevent topology generation entirely.
