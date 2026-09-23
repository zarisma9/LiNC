# LiNC Parameterization Workflow

This directory contains a stepwise workflow for generating Amber-compatible parameters for noncanonical amino acids or modified residues using Gaussian 16 and AmberTools.

The workflow starts from a SMILES string and residue name, generates alpha- and beta-like conformers, performs HF/6-31G* ESP calculations, fits RESP charges using both conformers, and finally generates `.prepi` and `.frcmod` files.

---

## Workflow Overview

```text
SMILES + RESNAME
      |
      v
conformer.py
      |
      +--> alpha_conf.mol2
      +--> beta_conf.mol2
      |
      v
make_gau_input.sh
      |
      +--> RESNAMEa.gau
      +--> RESNAMEb.gau
      |
      v
run_g16.sh
      |
      +--> RESNAMEa.out
      +--> RESNAMEb.out
      |
      v
esp_extract.sh
      |
      +--> RESNAME.esp
      +--> RESNAME.crg
      +--> RESNAME_final.ac
      |
      v
make_parameter.sh
      |
      +--> RESNAME.prepi
      +--> RESNAME.frcmod
```

---

## Files

### `conformer.py`

Generates alpha- and beta-like conformers from an input SMILES string.

Expected outputs:

```text
alpha_conf.mol2
beta_conf.mol2
```

The script also prints the backbone phi and psi atom indices used by `make_gau_input.sh`.

---

### `make_gau_input.sh`

Generates Gaussian input files for the alpha and beta conformers.

Usage:

```bash
chmod +x make_gau_input.sh

./make_gau_input.sh \
'O=C(NC)[C@H](CCCCNC(CCCCCCCCCCCCC(O)=O)=O)NC(C)=O' \
D4K
```

General format:

```bash
./make_gau_input.sh 'SMILES_CODE' RESNAME
```

Expected outputs:

```text
RESNAMEa.gau
RESNAMEb.gau
```

The current Gaussian calculation uses:

```text
HF/6-31G*
SCF=Tight
Pop=MK
opt=modredundant
```

The backbone phi and psi dihedral angles are frozen during optimization.

---

### `run_g16.sh`

Runs Gaussian 16 calculations for both conformers.

Usage:

```bash
chmod +x run_g16.sh

./run_g16.sh D4K
```

Alternatively:

```bash
RESNAME=D4K ./run_g16.sh
```

Expected outputs:

```text
D4Ka.out
D4Kb.out
```

---

### `esp_extract.sh`

Extracts ESP information from the Gaussian outputs and performs a two-stage RESP charge fitting using both conformers.

Usage:

```bash
chmod +x esp_extract.sh

./esp_extract.sh D4K
```

Alternatively:

```bash
RESNAME=D4K ./esp_extract.sh
```

Main outputs include:

```text
D4Ka.esp
D4Kb.esp
D4K.esp
D4K-step1.respin
D4K-step2.respin
D4K-step1.respout
D4K-step2.respout
D4K.crg
D4K_final.ac
```

The current workflow assumes:

```text
nmol = 2
net charge = 0
```

If the target residue has a different net charge, modify the relevant `-nc` option and charge settings before running the workflow.

---

### `make_parameter.sh`

Generates the final Amber residue library and force-field parameter files.

Usage:

```bash
chmod +x make_parameter.sh

./make_parameter.sh D4K
```

Alternatively:

```bash
RESNAME=D4K ./make_parameter.sh
```

Expected outputs:

```text
D4K.prepi
D4K.frcmod
D4K.res
mainchain.D4K
```

`prepgen` uses the main-chain definition generated inside `make_parameter.sh`.

The current omitted atom names are:

```text
C101
C102
HAC1
HAC2
HAC3
O100
N100
HNM1
HNM2
HNM3
HNM4
C103
```

These names may need to be modified if a different residue or atom-naming scheme is used.

---

## Requirements

The following software must be available in the environment:

- Python 3
- Gaussian 16
- AmberTools

Required AmberTools commands include:

```text
antechamber
espgen
respgen
resp
prepgen
parmchk2
```

`AMBERHOME` must also be correctly defined.

Example:

```bash
echo $AMBERHOME
```

Additional Python packages may be required by `conformer.py`.

---

## Example: D4K

A complete D4K parameterization can be run as:

```bash
SMILES='O=C(NC)[C@H](CCCCNC(CCCCCCCCCCCCC(O)=O)=O)NC(C)=O'
RESNAME='D4K'

./make_gau_input.sh "$SMILES" "$RESNAME"

./run_g16.sh "$RESNAME"

./esp_extract.sh "$RESNAME"

./make_parameter.sh "$RESNAME"
```

Final Amber parameter files:

```text
D4K.prepi
D4K.frcmod
```

---

## Important Notes

1. The workflow currently uses two conformers for RESP fitting: alpha and beta.
2. Gaussian ESP calculations use HF/6-31G* with Merz-Kollman electrostatic potentials.
3. The current RESP workflow assumes a neutral residue (`net charge = 0`).
4. Atom names used by `prepgen` are residue-specific and may need to be changed for other ncAAs.
5. Check Gaussian output files for successful termination before running `esp_extract.sh`.
6. Inspect the final `.prepi` and `.frcmod` files before using them in production MD simulations.

---

## Recommended Directory Structure

```text
01. Parameterization of LiNC/
├── README.md
├── conformer.py
├── make_gau_input.sh
├── run_g16.sh
├── esp_extract.sh
└── make_parameter.sh
```

Generated intermediate and output files can be kept in a separate working directory for each residue if multiple ncAAs are parameterized.

Example:

```text
work/
├── D4K/
├── M4K/
├── D8K/
└── M8K/
```

This prevents intermediate files such as `alpha_conf.mol2`, `beta_conf.mol2`, `qin`, and RESP outputs from different residues from being overwritten.
