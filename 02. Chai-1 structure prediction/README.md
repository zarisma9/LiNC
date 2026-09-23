# Chai-1 Covalent Complex Input Workflow

This directory contains a small workflow for preparing and running Chai-1 structure prediction with user-defined molecular entities and covalent-bond restraints.

The input generator creates three files in a user-specified working directory:

```text
input_sequence.fa
input_covalent_bond.restraints
run.py
```

The workflow is designed so that users can define any number of chains/entities and explicitly specify the covalent connection(s) used for Chai-1 inference.

---

## Workflow Overview

```text
User-defined entities
  |
  |-- Chain A: protein / DNA / RNA / ligand
  |-- Chain B: protein / DNA / RNA / ligand
  |-- Chain C: ...
  |
  +-- User-defined covalent bond restraint(s)
          |
          v
make_chai_input.py
          |
          +--> input_sequence.fa
          +--> input_covalent_bond.restraints
          +--> run.py
                    |
                    v
                 Chai-1
                    |
                    v
                 output/
```

---

## Files

### `make_chai_input.py`

Generates all Chai-1 input files from command-line arguments.

The script requires:

- a working directory specified with `--default-path`
- at least one `--chain`
- at least one `--bond`

No output files are created unless the required chain and covalent-bond information passes validation.

---

### `input_sequence.fa`

Contains all molecular entities used for Chai-1 inference.

Example:

```text
>protein|A
H(AIB)EGTFTSDVSSYLEGQAAGEFIAWLVRGRG
>protein|B
DAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDH...
>ligand|C
OC(CCCCCCCCCCCCC(NCCCC)=O)=O
```

Each entity is defined by:

```text
>TYPE|CHAIN_ID
SEQUENCE_OR_SMILES
```

Supported entity types in `make_chai_input.py` are:

```text
protein
dna
rna
ligand
```

---

### `input_covalent_bond.restraints`

Defines covalent connections between entities.

The file format is:

```csv
chainA,res_idxA,chainB,res_idxB,connection_type,confidence,min_distance_angstrom,max_distance_angstrom,comment,restraint_id
```

Example:

```csv
A,G20@CA,C,@C20,covalent,1.0,1.3,1.7,protein-ligand,bond1
```

The user explicitly defines:

| Field | Description |
|---|---|
| `chainA` | First chain/entity |
| `res_idxA` | Residue/atom identifier in the first entity |
| `chainB` | Second chain/entity |
| `res_idxB` | Residue/atom identifier in the second entity |
| `connection_type` | Connection type, e.g. `covalent` |
| `confidence` | Restraint confidence |
| `min_distance_angstrom` | Minimum allowed distance |
| `max_distance_angstrom` | Maximum allowed distance |
| `comment` | User-defined description |
| `restraint_id` | Unique restraint identifier |

Multiple covalent restraints can be supplied by repeating `--bond`.

---

### `run.py`

Runs Chai-1 inference using:

```text
input_sequence.fa
input_covalent_bond.restraints
```

The generated script uses the directory supplied to `make_chai_input.py` as its default working path.

The path can also be overridden when running inference:

```bash
python run.py --default-path /path/to/working_directory
```

The output directory is:

```text
output/
```

If `output/` already exists, it is removed before a new inference run.

The current inference settings are:

```python
num_trunk_recycles=3
num_diffn_timesteps=200
seed=99
device="cuda:0"
use_esm_embeddings=True
num_diffn_samples=1
```

The workflow also saves:

```text
pae.pt
pde.pt
plddt.pt
```

inside the output directory.

---

## Usage

### 1. Make the generator executable

```bash
chmod +x make_chai_input.py
```

Alternatively, run it directly with Python:

```bash
python make_chai_input.py ...
```

---

## Basic Example

The general syntax is:

```bash
python make_chai_input.py \
  --default-path WORKING_DIRECTORY \
  --chain CHAIN_ID TYPE SEQUENCE_OR_SMILES \
  --chain CHAIN_ID TYPE SEQUENCE_OR_SMILES \
  --bond CHAIN_A RES_IDX_A CHAIN_B RES_IDX_B CONNECTION_TYPE CONFIDENCE MIN_DIST MAX_DIST COMMENT RESTRAINT_ID
```

Example:

```bash
python make_chai_input.py \
  --default-path /home/guest_1/LDM/data/HSK/260528_lipid/D4K \
  --chain A protein 'H(AIB)EGTFTSDVSSYLEGQAAGEFIAWLVRGRG' \
  --chain B protein 'DAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDH...' \
  --chain C ligand 'OC(CCCCCCCCCCCCC(NCCCC)=O)=O' \
  --bond A 'G20@CA' C '@C20' covalent 1.0 1.3 1.7 protein-ligand bond1
```

This generates:

```text
/home/guest_1/LDM/data/HSK/260528_lipid/D4K/
├── input_sequence.fa
├── input_covalent_bond.restraints
└── run.py
```

---

## Adding Multiple Chains

`--chain` can be repeated as many times as needed.

Example:

```bash
python make_chai_input.py \
  --default-path /path/to/example \
  --chain A protein 'SEQUENCE_A' \
  --chain B protein 'SEQUENCE_B' \
  --chain C ligand 'SMILES_C' \
  --chain D dna 'ATGCGTATGCGT' \
  --bond A 'G20@CA' C '@C20' covalent 1.0 1.3 1.7 protein-ligand bond1
```

Each chain ID must be unique.

---

## Adding Multiple Covalent Bonds

`--bond` can also be repeated.

Example:

```bash
python make_chai_input.py \
  --default-path /path/to/example \
  --chain A protein 'SEQUENCE_A' \
  --chain B protein 'SEQUENCE_B' \
  --chain C ligand 'SMILES_C' \
  --bond A 'G20@CA' C '@C20' covalent 1.0 1.3 1.7 protein-ligand bond1 \
  --bond B 'K50@NZ' C '@C1' covalent 1.0 1.3 1.7 protein-ligand bond2
```

Every chain referenced by a bond must first be defined with `--chain`.

---

## Required Covalent-Bond Definition

At least one `--bond` argument is required.

For example, this is intentionally incomplete:

```bash
python make_chai_input.py \
  --default-path /path/to/example \
  --chain A protein 'SEQUENCE_A'
```

The script exits without creating the three output files because no covalent connection was specified.

This prevents accidental generation of an input set without the intended covalent-bond restraint.

---

## Overwriting Existing Files

By default, the generator does not overwrite:

```text
input_sequence.fa
input_covalent_bond.restraints
run.py
```

To intentionally overwrite existing files, add:

```bash
--force
```

Example:

```bash
python make_chai_input.py \
  --default-path /path/to/example \
  --chain A protein 'SEQUENCE_A' \
  --chain C ligand 'SMILES_C' \
  --bond A 'G20@CA' C '@C20' covalent 1.0 1.3 1.7 protein-ligand bond1 \
  --force
```

---

## Running Chai-1

After generating the input files:

```bash
cd /path/to/working_directory
python run.py
```

or:

```bash
python run.py --default-path /path/to/working_directory
```

The resulting directory structure is expected to look like:

```text
working_directory/
├── input_sequence.fa
├── input_covalent_bond.restraints
├── run.py
└── output/
    ├── *.cif
    ├── scores.model_idx_0.npz
    ├── pae.pt
    ├── pde.pt
    └── plddt.pt
```

---

## Requirements

The inference environment requires:

- Python 3
- Chai-1 / `chai_lab`
- PyTorch
- NumPy
- CUDA-compatible GPU for the current `cuda:0` setting

The following Python import must work:

```python
from chai_lab.chai1 import run_inference
```

---

## Recommended Directory Structure

For multiple targets, use one working directory per system:

```text
Chai1/
├── make_chai_input.py
├── D4K/
│   ├── input_sequence.fa
│   ├── input_covalent_bond.restraints
│   ├── run.py
│   └── output/
├── M4K/
│   ├── input_sequence.fa
│   ├── input_covalent_bond.restraints
│   ├── run.py
│   └── output/
└── D8K/
    ├── input_sequence.fa
    ├── input_covalent_bond.restraints
    ├── run.py
    └── output/
```

Keeping each target in a separate directory prevents input and output files from different systems from being mixed.

---

## Important Notes

1. Chain IDs must be unique.
2. Every chain referenced in a covalent restraint must be defined with `--chain`.
3. At least one covalent-bond restraint is required.
4. `confidence` must be between `0.0` and `1.0`.
5. Minimum and maximum restraint distances must be positive.
6. `min_distance_angstrom` cannot be greater than `max_distance_angstrom`.
7. Atom/residue identifiers such as `G20@CA` and `@C20` must match the intended atoms in the corresponding Chai-1 entity.
8. The existing `output/` directory is deleted before each inference run.
9. Use `--force` only when existing generated inputs are intentionally being replaced.
