# LiNC

**A four-stage computational workflow for noncanonical residue parameterization, Chai-1 structure prediction, AMBER molecular dynamics, and binding-energy analysis.**

LiNC brings together scripts for generating custom residue parameters, preparing covalently connected molecular complexes for structure prediction, running molecular dynamics (MD), and performing structural analysis and MM/GBSA–MM/PBSA calculations with per-residue energy decomposition.

This root README explains how the four workflow directories connect. The linked stage READMEs provide detailed command-line options and file descriptions. The examples below use `D4K` as the residue name and `D4K_seed1` as the simulation identifier; these are example identifiers, not automatic selections of a chemical structure or random seed.

## Workflow at a glance

| Stage | Directory and detailed guide | Main inputs | Main outputs |
| :--- | :--- | :--- | :--- |
| **01** | [Parameterization of LiNC](01.%20Parameterization%20of%20LiNC/README.md) | Capped residue SMILES and `RESNAME` | `.prepi` residue template and `.frcmod` parameters |
| **02** | [Chai-1 structure prediction](02.%20Chai-1%20structure%20prediction/README.md) | Molecular sequences/SMILES and user-defined covalent restraints | Predicted `.cif` structure and confidence-related outputs |
| **03** | [Molecular dynamics simulations](03.%20Molecular%20dynamics%20simulations/README.md) | Predicted structure and matching residue parameters | Standard/HMR topologies, coordinates, and MD trajectories |
| **04** | [MMPBSA binding free-energy calculations](04.%20MMPBSA%20binding%20free-energy%20calculations/README.md) | Production trajectory and matching topology | RMSD/distance data, MM/GBSA–MM/PBSA results, and residue decomposition |

```text
Capped residue SMILES + RESNAME          Sequences / ligand SMILES + covalent bonds
                 |                                          |
                 v                                          v
01. Parameterization of LiNC             02. Chai-1 structure prediction
                 |                                          |
                 v                                          v
      RESNAME.prepi + RESNAME.frcmod              Predicted complex (.cif)
                 |                                          |
                 +--------------------+---------------------+
                                      |
                                      v
                       03. Molecular dynamics simulations
                                      |
                       CIF -> merged PDB -> AMBER PDB
                                      |
                   VERIFY NONCANONICAL RESIDUE / TEMPLATE MATCH
                                      |
                   tleap -> standard topology -> HMR topology
                                      |
                  Minimization -> heating -> equilibration -> MD
                                      |
                                      v
                  04. MMPBSA binding free-energy calculations
                                      |
                 Autoimaging -> structural analysis -> GB/PB analysis
                                      |
                      Binding-energy and decomposition outputs
```

Stages 01 and 02 provide separate inputs to Stage 03: the Chai-1 input generator does not consume the `.prepi` or `.frcmod` files. Stage 04 consumes the trajectory and topology produced by Stage 03.

## Repository layout

```text
LiNC/
├── README.md
├── 01. Parameterization of LiNC/
│   ├── README.md
│   ├── conformer.py
│   ├── make_gau_input.sh
│   ├── run_g16.sh
│   ├── esp_extract.sh
│   ├── make_parameter.sh
│   └── parameter/
│       └── *.prepi / *.frcmod
├── 02. Chai-1 structure prediction/
│   ├── README.md
│   └── make_chai_input.py
├── 03. Molecular dynamics simulations/
│   ├── README.md
│   ├── modify_pdb.py
│   ├── make_amber_pdb.sh
│   ├── make_tleap_input.sh
│   ├── make_hmass_prmtop.sh
│   ├── make_MD_input.sh
│   └── run_MD.sh
└── 04. MMPBSA binding free-energy calculations/
    ├── README.md
    └── run_MMPBSA.sh
```

The [parameter directory](01.%20Parameterization%20of%20LiNC/parameter/) contains `.prepi`/`.frcmod` pairs named `D0Z`, `D4K`, `D4Z`, `D8K`, `DEK`, `M0Z`, `M4K`, `M4Z`, `M7K`, `M8K`, and `MEK`. Before reusing a pair, verify that its residue definition and atom naming match the intended system. A matching filename alone is not a parameter-validation check.

## Software requirements

| Stage | Required software or commands |
| :--- | :--- |
| 01 | Bash, Python 3, RDKit, Open Babel Python bindings (`openbabel.pybel`), Gaussian 16 (`g16`), `antechamber`, `espgen`, `respgen`, `resp`, `prepgen`, and `parmchk2` |
| 02 | Python 3, Chai-1 (`chai_lab`), PyTorch, and NumPy; the generated runner selects `cuda:0` |
| 03 | Bash, Python 3, `obabel`, `pdb4amber`, `tleap`, `parmed`, `mpirun`, `pmemd.MPI`, and `pmemd.cuda` |
| 04 | Bash, `cpptraj`, `ante-MMPBSA.py`, `MMPBSA.py.MPI`, and `mpirun` |

Activate the appropriate environment before each stage. Stage 01 also requires `AMBERHOME` and `${AMBERHOME}/dat/leap/parm/parm10.dat`. Software installation, environment activation, and cluster resource allocation are not performed by these workflow scripts.

The examples use `bash script.sh`, so executable permissions are not required. For direct execution with `./script.sh`, first run `chmod +x script.sh`. Run CPU/GPU calculations within an appropriate allocation on a cluster; the scripts do not submit Slurm jobs or request resources themselves.

## Working directories and identifiers

Keep generated files outside the source directories and use a separate working directory for every target and MD replicate. Several scripts use fixed intermediate or output filenames.

The following setup is an example layout for the commands in this README, not a directory structure created automatically by the repository:

```bash
git clone https://github.com/zarisma9/LiNC.git
cd LiNC

export REPO_ROOT="$(pwd)"
export WORK_ROOT="${HOME}/LiNC_runs"
export RESNAME="D4K"
export SYSTEM_ID="${RESNAME}_seed1"

export PARAM_SRC="${REPO_ROOT}/01. Parameterization of LiNC"
export CHAI_SRC="${REPO_ROOT}/02. Chai-1 structure prediction"
export MD_SRC="${REPO_ROOT}/03. Molecular dynamics simulations"
export ANALYSIS_SRC="${REPO_ROOT}/04. MMPBSA binding free-energy calculations"

export PARAM_WORK="${WORK_ROOT}/parameters/${RESNAME}"
export PRED_DIR="${WORK_ROOT}/predictions/${SYSTEM_ID}"
export MD_WORK="${WORK_ROOT}/md/${SYSTEM_ID}"

mkdir -p "$PARAM_WORK" "$MD_WORK"
```

Keep these variables available when moving between stages. Quote repository paths because the source-directory names contain spaces. Use a working path without whitespace: some generated external-tool input files insert paths without additional quoting.

`RESNAME` identifies the custom residue, whereas `SYSTEM_ID` identifies a particular simulation. The MD runner calls the latter `SYS_NAME`; the analysis runner calls it `SYSTEM_ID`. Pass the same identifier to both.

## 01. Parameterization of LiNC

The parameterization sequence is:

```text
make_gau_input.sh  ->  run_g16.sh  ->  esp_extract.sh  ->  make_parameter.sh
  calls conformer.py   runs a/b      two-conformer       .prepi + .frcmod
                       Gaussian     two-stage RESP
```

`make_gau_input.sh` runs `conformer.py` internally; a separate conformer-generation command is not required. The generator creates alpha- and beta-like conformers, then Gaussian inputs with backbone dihedral restraints. The current calculation settings include `HF/6-31G*`, `SCF=Tight`, `Pop=MK`, and `opt=modredundant`.

> [!IMPORTANT]
> **Check the backbone-log strings before generating new parameters.** In the inspected repository revision, `conformer.py` prints `Backbone phi atoms: ids` and `Backbone psi atoms: ids`, but `make_gau_input.sh` searches for `Backbone φ atoms: ids` and `Backbone ψ atoms: ids`. These strings do not match, so atom-index extraction can stop the workflow. Make the search strings agree with the conformer output in your working copy before running the commands below. This README does not patch the scripts.
>
> Source: [conformer.py](01.%20Parameterization%20of%20LiNC/conformer.py) and [make_gau_input.sh](01.%20Parameterization%20of%20LiNC/make_gau_input.sh).

Prepare an isolated parameterization directory:

```bash
cp "${PARAM_SRC}/conformer.py" \
   "${PARAM_SRC}/make_gau_input.sh" \
   "${PARAM_SRC}/run_g16.sh" \
   "${PARAM_SRC}/esp_extract.sh" \
   "${PARAM_SRC}/make_parameter.sh" \
   "$PARAM_WORK/"

cd "$PARAM_WORK"

# Review the log-string mismatch above before running.
SMILES='O=C(NC)[C@H](CCCCNC(CCCCCCCCCCCCC(O)=O)=O)NC(C)=O'

bash make_gau_input.sh "$SMILES" "$RESNAME"
bash run_g16.sh "$RESNAME"
```

Inspect both Gaussian outputs before proceeding to RESP and parameter generation:

```bash
bash esp_extract.sh "$RESNAME"
bash make_parameter.sh "$RESNAME"
```

The main handoff files are `${RESNAME}.prepi` and `${RESNAME}.frcmod`. Intermediate files include `alpha_conf.mol2`, `beta_conf.mol2`, the two Gaussian inputs/outputs, ESP files, RESP inputs/outputs, and `${RESNAME}_final.ac`.

The current scripts use neutral charge settings (`0 1` for Gaussian body extraction, `-nc 0` in the RESP preparation step, and `CHARGE 0.0` in the main-chain definition). They also use specific cap atom names in `OMIT_NAME`. Review these settings for the intended charge state and atom-naming scheme rather than treating the example as a general parameterization of any ncAA.

**Alternative: use an existing parameter pair.** After checking that the bundled pair describes the intended residue, copy it into the parameter working directory instead of regenerating it:

```bash
cp "${PARAM_SRC}/parameter/${RESNAME}.prepi" "$PARAM_WORK/"
cp "${PARAM_SRC}/parameter/${RESNAME}.frcmod" "$PARAM_WORK/"
```

## 02. Chai-1 structure prediction

`make_chai_input.py` writes `input_sequence.fa`, `input_covalent_bond.restraints`, and `run.py` into `--default-path`. Repeat `--chain` to add entities and repeat `--bond` to add connections.

```text
--chain CHAIN_ID TYPE SEQUENCE_OR_SMILES
--bond CHAIN_A RES_IDX_A CHAIN_B RES_IDX_B CONNECTION_TYPE CONFIDENCE MIN_DIST MAX_DIST COMMENT RESTRAINT_ID
```

The generator accepts the entity types `protein`, `dna`, `rna`, and `ligand`. At least one `--bond` is required; omitting it prevents generation of the input set. The later MD scripts are configured for the protein/modified-peptide workflow and are not automatically adapted to every entity type accepted by this generator.

The example below follows the Stage 02 README. Set `RECEPTOR_SEQUENCE` to the **complete** chain B sequence first; the abbreviated sequence in the stage documentation is not a complete input. Confirm that the chosen ligand SMILES, attachment atoms, and Stage 01 parameter definition describe the intended chemistry.

```bash
: "${RECEPTOR_SEQUENCE:?Set RECEPTOR_SEQUENCE to the complete chain B sequence}"

PEPTIDE_SEQUENCE='H(AIB)EGTFTSDVSSYLEGQAAGEFIAWLVRGRG'
LIGAND_SMILES='OC(CCCCCCCCCCCCC(NCCCC)=O)=O'

python3 "${CHAI_SRC}/make_chai_input.py" \
  --default-path "$PRED_DIR" \
  --chain A protein "$PEPTIDE_SEQUENCE" \
  --chain B protein "$RECEPTOR_SEQUENCE" \
  --chain C ligand "$LIGAND_SMILES" \
  --bond A 'G20@CA' C '@C20' covalent 1.0 1.3 1.7 protein-ligand bond1
```

The example restraint is written as:

```csv
chainA,res_idxA,chainB,res_idxB,connection_type,confidence,min_distance_angstrom,max_distance_angstrom,comment,restraint_id
A,G20@CA,C,@C20,covalent,1.0,1.3,1.7,protein-ligand,bond1
```

Attachment identifiers are user-supplied, not inferred or chemically validated by this wrapper. Review them for every system. The generator checks basic fields and referenced chain IDs; it does not establish that a given atom is the correct chemical attachment site.

Run prediction in the Chai-1 environment:

```bash
python3 "${PRED_DIR}/run.py"
```

The generated runner uses `num_trunk_recycles=3`, `num_diffn_timesteps=200`, `seed=99`, `num_diffn_samples=1`, `use_esm_embeddings=True`, and `device="cuda:0"`. Outputs are written to `${PRED_DIR}/output/`, including predicted CIF files, score files, `pae.pt`, `pde.pt`, and `plddt.pt`. Confirm the selected CIF filename before Stage 03; the downstream example uses `pred.model_idx_0.cif`.

> [!WARNING]
> The generator refuses to replace existing inputs unless `--force` is supplied. **The generated `run.py` behaves differently: it deletes an existing `output/` directory before inference.** Back up previous predictions before rerunning. Changing a directory name from `seed1` to `seed2` does not change the hard-coded Chai-1 random seed.

## 03. Molecular dynamics simulations

This stage combines the Stage 01 parameters with the Stage 02 structure:

```text
make_amber_pdb.sh
        |
        v
CHECK CUSTOM RESIDUE NAME AND TEMPLATE
        |
        v
make_tleap_input.sh -> tleap -f tleap.in
        |
        v
make_hmass_prmtop.sh -> make_MD_input.sh -> run_MD.sh
```

### Convert the prediction and verify the custom residue

`make_amber_pdb.sh` expects `modify_pdb.py` in the current working directory. Copy the helper and the checked parameter pair into the MD directory:

```bash
cp "${MD_SRC}/modify_pdb.py" "$MD_WORK/"
cp "${PARAM_WORK}/${RESNAME}.prepi" "$MD_WORK/"
cp "${PARAM_WORK}/${RESNAME}.frcmod" "$MD_WORK/"
cd "$MD_WORK"

bash "${MD_SRC}/make_amber_pdb.sh" \
  "${PRED_DIR}/output/pred.model_idx_0.cif" \
  "${SYSTEM_ID}.pdb" \
  "${RESNAME}_final.pdb" \
  A 20 "$RESNAME" \
  "${RESNAME}_amber.pdb"
```

The seven arguments specify the input CIF, converted PDB, modified PDB, target chain, target residue number, new residue name, and final AMBER PDB. The script runs Open Babel, `modify_pdb.py`, and `pdb4amber` in that order. It supplies `-l LIG` to the helper; the helper matches ligand residue names by prefix. Check the input carefully when multiple residues share that prefix. A missing target residue or ligand produces a warning in the helper, not a hard stop.

> [!CAUTION]
> **Do not run `tleap` until the noncanonical residue name in the final `*_amber.pdb` has been checked against the internal residue name of the intended `.prepi` template.**
>
> For the D4K example, the merged residue must be named `D4K`, and the loaded template must define `D4K`. Also check atom names and select the corresponding `.frcmod` file. Matching filenames alone does not establish that the PDB and residue definition agree.

Inspect the target atoms using the final PDB's actual chain and residue numbering:

```bash
awk '($1 == "ATOM" || $1 == "HETATM") &&
     substr($0,22,1) == "A" && (substr($0,23,4) + 0) == 20 {print}' \
  "${RESNAME}_amber.pdb"

# Inspect the template's internal residue-name record.
awk '$2 == "INT" {print; exit}' "${RESNAME}.prepi"
```

An empty target selection is not a successful check. Confirm the final numbering and inspect the intended residue before continuing. The expected residue-name field in this example is `D4K`; the template record is `D4K INT 0`.

### Build the solvated system and HMR topology

```bash
# Ion counts below are the repository example; choose them for your system.
K_COUNT=88
CL_COUNT=71

bash "${MD_SRC}/make_tleap_input.sh" \
  "${RESNAME}.prepi" "${RESNAME}.frcmod" "${RESNAME}_amber.pdb" \
  "$K_COUNT" "$CL_COUNT" \
  "${SYSTEM_ID}.prmtop" "${SYSTEM_ID}.inpcrd"

tleap -f tleap.in
```

**`make_tleap_input.sh` generates the input only.** The separate `tleap` command must finish successfully before HMR preparation. Inspect the output for unrecognized residues, missing atom types, or missing parameters.

The current generator writes the following source lines in this order; this documents the existing template rather than changing its force-field choices:

```text
source leaprc.protein.ff19SB
source leaprc.water.opc
source leaprc.mimetic.ff15ipq
```

It then uses `solvateOct m OPCBOX 12.0` and the supplied K+/Cl- counts. Review the generated `tleap.in` for the intended system. Ion counts are explicit inputs, not values calculated automatically from a target salt concentration.

After confirming the standard topology and coordinates:

```bash
bash "${MD_SRC}/make_hmass_prmtop.sh" \
  "${SYSTEM_ID}.prmtop" "${SYSTEM_ID}.hmass.prmtop"

bash "${MD_SRC}/make_MD_input.sh"
```

The HMR wrapper writes `hmass_parmed.in` and runs `parmed -i hmass_parmed.in`. The MD input generator writes seven input files with the following settings:

| Step | Input file | Protocol in the template | Topology |
| :--- | :--- | :--- | :--- |
| 1 | `01_Min1.in` | Minimization; maximum 3,000 cycles; restraint weight 10.0 | Standard |
| 2 | `02_Min2.in` | Minimization; maximum 10,000 cycles; restraint weight 2.0 | Standard |
| 3 | `03_Heat_NVT.in` | NVT heating, 0–310 K; 1 fs; 25 ps; restraint weight 5.0 | Standard |
| 4 | `04_Equil_NPT_1.in` | NPT at 310 K; 2 fs; 200 ps; restraint weight 2.0 | Standard |
| 5 | `05_Equil_NPT_2.in` | NPT at 310 K; 2 fs; 200 ps; restraint weight 0.5 | Standard |
| 6 | `06_Equil_NPT_3.in` | Unrestrained NPT at 310 K; 2 fs; 500 ps | Standard |
| 7 | `07_Prod_HMR_4fs.in` | Unrestrained NPT at 310 K; 4 fs; 1 microsecond target | HMR |

Review the hard-coded restraint mask `!:WAT,Na+,K+,Cl-` and the generated protocol before running MD. The production input specifies `dt=0.004` and `nstlim=250000000`; creating this input does not establish that the simulation has reached that duration.

### Run MD

From the MD working directory, inside the appropriate CPU/GPU allocation:

```bash
MPI_NP=4 bash "${MD_SRC}/run_MD.sh" "$SYSTEM_ID"
```

Steps 1–2 use `mpirun -np 4 pmemd.MPI` by default. Steps 3–7 use `pmemd.cuda`; only Step 7 switches to `${SYSTEM_ID}.hmass.prmtop`. `MPI_NP`, `DO_CPU`, and `DO_GPU` can be supplied as environment variables.

The main Stage 04 handoff files are `${SYSTEM_ID}.hmass.prmtop` and `07_Prod_HMR.nc`. Additional outputs include each step's log, restart, and MD trajectory files. The runner checks restart-file existence/size; inspect the simulation logs and completion status as well. It starts from Step 1 on each invocation and does not automatically resume an interrupted production run.

## 04. MMPBSA binding free-energy calculations

`run_MMPBSA.sh` takes `SYSTEM_ID` and `BASE_DIR`, expects `${BASE_DIR}/${SYSTEM_ID}.hmass.prmtop` and `${BASE_DIR}/07_Prod_HMR.nc`, and writes its analysis into `${BASE_DIR}/BA/`.

```bash
REC_MASK="32-616" \
PEP_MASK="1-31" \
LIGAND_MASK=":1-31" \
START_FRAME=1500 \
END_FRAME=2000 \
INTERVAL=5 \
MMPBSA_CORES=8 \
bash "${ANALYSIS_SRC}/run_MMPBSA.sh" "$SYSTEM_ID" "$MD_WORK"
```

This runs trajectory autoimaging, receptor-aligned RMSD and minimum-distance calculations, complex/receptor/ligand topology generation, and MM/GBSA–MM/PBSA with per-residue decomposition.

| Setting | Current value | How to change it |
| :--- | :--- | :--- |
| Receptor residues | `REC_MASK=32-616` | Environment variable; omit the leading colon |
| Peptide residues | `PEP_MASK=1-31` | Environment variable; omit the leading colon |
| MMPBSA ligand | `LIGAND_MASK=:1-31` | Environment variable; include the leading colon |
| Analysis frames | `1500–2000`, interval `5` | `START_FRAME`, `END_FRAME`, and `INTERVAL` |
| MPI processes | `8` | `MMPBSA_CORES` |
| Removed solvent/ions | `:WAT,HOH,TIP3,Na+,Cl-,K+,Mg2+,Zn2+` | `STRIP_MASK` |
| Residue used for the lipid-named RMSD output | `:20` | Edit the script; currently hard-coded |
| Decomposition residue range | `print_res='1-616'` | Edit the script; currently hard-coded |

> [!IMPORTANT]
> **Verify residue numbering against the topology used for the analysis.** Changing `PEP_MASK` does not automatically change `LIGAND_MASK`, the hard-coded residue `:20`, or `print_res='1-616'`. Check each selection and confirm that the trajectory contains the requested frames.

The current analysis uses mass-weighted receptor backbone alignment (`N,CA,C`) against the first frame. Peptide backbone RMSD and the residue-20 RMSD use `nofit` after that alignment. The `*_rmsd_lipid.dat` selection is **all atoms of residue 20**, not a side-chain-only or heavy-atom-only selection. The minimum-distance masks exclude names matching `H=`.

For MMPBSA, the ligand is the entire selected peptide (`:1-31` by default), not residue 20 alone. The residue decomposition output is a separate result from the total complex binding-energy output.

The generated analysis inputs use `mbondi3` radii, GB `igb=8` and `saltcon=0.150`, and PB `indi=2.0`, `exdi=80.0`, `istrng=0.150`, `inp=2`, and `radiopt=0`. Decomposition uses `idecomp=1`, `dec_verbose=1`, and `csv_format=1`. These settings are preserved in the script; consult the [Stage 04 guide](04.%20MMPBSA%20binding%20free-energy%20calculations/README.md) for the full input blocks.

## Main analysis outputs

```text
${MD_WORK}/
├── ${SYSTEM_ID}.prmtop
├── ${SYSTEM_ID}.hmass.prmtop
├── ${SYSTEM_ID}.inpcrd
├── 07_Prod_HMR.nc
└── BA/
    ├── center.in
    ├── rmsd.in
    ├── mmpbsa.in
    ├── ${SYSTEM_ID}_centered.nc
    ├── ${SYSTEM_ID}_rmsd_bb.dat
    ├── ${SYSTEM_ID}_rmsd_peptide.dat
    ├── ${SYSTEM_ID}_rmsd_lipid.dat
    ├── ${SYSTEM_ID}_mindist_ha.dat
    ├── ${SYSTEM_ID}_HMR_dry_complex.prmtop
    ├── ${SYSTEM_ID}_HMR_receptor.prmtop
    ├── ${SYSTEM_ID}_HMR_ligand.prmtop
    ├── FINAL_RESULTS_MMPBSA.dat
    └── FINAL_DECOMP_MMPBSA.dat
```

`FINAL_RESULTS_MMPBSA.dat` is the main MM/GBSA–MM/PBSA result file. `FINAL_DECOMP_MMPBSA.dat` contains the residue-level decomposition output. The RMSD and minimum-distance files provide the accompanying structural analysis.

## Before running a new system

**Chemical identity and naming.** Verify that the parameterized residue, predicted ligand/attachment, merged residue, and final PDB describe the intended system. Check the final noncanonical residue name against the `.prepi` definition before `tleap`, and inspect atom names and the selected `.frcmod`.

**Target-specific settings.** Review charge state, cap atom names, covalent attachment identifiers, source files in `tleap.in`, ion counts, MD restraints, and every analysis mask. Defaults in these scripts are not automatically inferred for a new system.

**Existing outputs.** Use isolated working directories and back up results before rerunning. Chai-1 removes its existing `output/`; MD commands use overwrite mode and fixed stage filenames; MMPBSA uses a fixed `BA/` directory and overwrite mode. Neither the MD nor the analysis script creates a new replicate directory automatically.

**Execution and provenance.** Verify each stage's logs before using its outputs in the next stage. Record the repository revision, installed software versions, input structures/parameters, runtime settings, and any local script changes alongside the results. File-existence checks are not a substitute for reviewing calculation completion and the intended molecular system.

---

For detailed options, use the four stage guides linked in the overview table. The scripts and examples documented here were inspected on **2026-09-23**, at repository revision `93f66e549757f6e26feed4a6cb3ec071165b1edd`.
