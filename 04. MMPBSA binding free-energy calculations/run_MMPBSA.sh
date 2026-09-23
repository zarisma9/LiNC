#!/usr/bin/env bash

set -euo pipefail

SYSTEM_ID="${1:-${SYSTEM_ID:-}}"
BASE_DIR="${2:-${BASE_DIR:-}}"

if [ -z "$SYSTEM_ID" ] || [ -z "$BASE_DIR" ]; then
    echo "Usage:"
    echo "  $0 SYSTEM_ID BASE_DIR"
    echo
    echo "Example:"
    echo "  $0 MEK_seed1 /home/guest_1/LDM/data/260304_lipid_MD/MEK_seed1"
    echo
    echo "Environment-variable usage is also supported:"
    echo "  SYSTEM_ID=MEK_seed1 BASE_DIR=/path/to/MEK_seed1 $0"
    exit 1
fi

BASE_DIR="$(cd "$BASE_DIR" && pwd)"
WORK_DIR="${BASE_DIR}/BA"

PARM_FILE="${BASE_DIR}/${SYSTEM_ID}.hmass.prmtop"
TRAJ_RAW="${BASE_DIR}/07_Prod_HMR.nc"
TRAJ_CENTERED="${WORK_DIR}/${SYSTEM_ID}_centered.nc"

REC_MASK="${REC_MASK:-32-616}"
PEP_MASK="${PEP_MASK:-1-31}"
LIGAND_MASK="${LIGAND_MASK:-:1-31}"

START_FRAME="${START_FRAME:-1500}"
END_FRAME="${END_FRAME:-2000}"
INTERVAL="${INTERVAL:-5}"

MMPBSA_CORES="${MMPBSA_CORES:-8}"

STRIP_MASK="${STRIP_MASK:-:WAT,HOH,TIP3,Na+,Cl-,K+,Mg2+,Zn2+}"

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

echo "====================================================="
echo " MMPBSA Analysis Pipeline"
echo "====================================================="
echo "System ID       : $SYSTEM_ID"
echo "Base directory  : $BASE_DIR"
echo "Work directory  : $WORK_DIR"
echo "Topology        : $PARM_FILE"
echo "Raw trajectory  : $TRAJ_RAW"
echo "Receptor mask   : $REC_MASK"
echo "Peptide mask    : $PEP_MASK"
echo "Ligand mask     : $LIGAND_MASK"
echo "Frame range     : $START_FRAME-$END_FRAME"
echo "Frame interval  : $INTERVAL"
echo "MPI cores       : $MMPBSA_CORES"
echo "====================================================="

REQUIRED_COMMANDS=(
    cpptraj
    ante-MMPBSA.py
    MMPBSA.py.MPI
    mpirun
)

for cmd in "${REQUIRED_COMMANDS[@]}"; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: $cmd is not available in PATH."
        exit 1
    fi
done

if [ ! -f "$PARM_FILE" ]; then
    echo "ERROR: Topology file not found: $PARM_FILE"
    exit 1
fi

if [ ! -f "$TRAJ_RAW" ]; then
    echo "ERROR: Production trajectory not found: $TRAJ_RAW"
    exit 1
fi

echo ">>> Step 1: Centering trajectory..."

cat > "${WORK_DIR}/center.in" <<EOF
parm ${PARM_FILE}
trajin ${TRAJ_RAW}
autoimage anchor :${REC_MASK}
trajout ${TRAJ_CENTERED}
run
quit
EOF

cpptraj -i "${WORK_DIR}/center.in"

if [ ! -s "$TRAJ_CENTERED" ]; then
    echo "ERROR: Centered trajectory was not created."
    exit 1
fi

echo ">>> Centered trajectory created: $TRAJ_CENTERED"

echo ">>> Step 2: Calculating RMSD and minimum distance..."

cat > "${WORK_DIR}/rmsd.in" <<EOF
parm ${PARM_FILE}
trajin ${TRAJ_CENTERED}

rms RMSDbb :${REC_MASK}@N,CA,C mass first out ${WORK_DIR}/${SYSTEM_ID}_rmsd_bb.dat

rms RMSDpeptide :${PEP_MASK}@N,CA,C mass first nofit out ${WORK_DIR}/${SYSTEM_ID}_rmsd_peptide.dat

rms RMSD_lipid :20 mass first nofit out ${WORK_DIR}/${SYSTEM_ID}_rmsd_lipid.dat

mindist mask1 (:${PEP_MASK}&!@H=) mask2 (:${REC_MASK}&!@H=) out ${WORK_DIR}/${SYSTEM_ID}_mindist_ha.dat name MinD_HA

run
quit
EOF

cpptraj -i "${WORK_DIR}/rmsd.in"

echo ">>> RMSD and minimum-distance calculations completed."

echo ">>> Step 3: Generating MMPBSA topology files..."

DRY_COMPLEX="${WORK_DIR}/${SYSTEM_ID}_HMR_dry_complex.prmtop"
RECEPTOR_TOP="${WORK_DIR}/${SYSTEM_ID}_HMR_receptor.prmtop"
LIGAND_TOP="${WORK_DIR}/${SYSTEM_ID}_HMR_ligand.prmtop"

ante-MMPBSA.py \
    -p "$PARM_FILE" \
    -c "$DRY_COMPLEX" \
    -r "$RECEPTOR_TOP" \
    -l "$LIGAND_TOP" \
    -n "$LIGAND_MASK" \
    -s "$STRIP_MASK" \
    --radii mbondi3

for file in "$DRY_COMPLEX" "$RECEPTOR_TOP" "$LIGAND_TOP"; do
    if [ ! -s "$file" ]; then
        echo "ERROR: Failed to create topology file: $file"
        exit 1
    fi
done

echo ">>> MMPBSA topology files created."

echo ">>> Step 4: Creating mmpbsa.in..."

cat > "${WORK_DIR}/mmpbsa.in" <<EOF
&general
  startframe=${START_FRAME}, endframe=${END_FRAME}, interval=${INTERVAL},
  verbose=1, keep_files=0,
  strip_mask='${STRIP_MASK}'
/
&gb
  igb=8, saltcon=0.150,
  surften=0.0072, surfoff=0.0
/
&pb
  indi=2.0, exdi=80.0,
  istrng=0.150, inp=2, radiopt=0
/
&decomp
  idecomp=1,
  print_res='1-616',
  dec_verbose=1,
  csv_format=1
/
EOF

echo ">>> Created: ${WORK_DIR}/mmpbsa.in"

echo ">>> Step 5: Running MMPBSA with ${MMPBSA_CORES} MPI processes..."

mpirun -np "$MMPBSA_CORES" MMPBSA.py.MPI -O \
    -i "${WORK_DIR}/mmpbsa.in" \
    -o "${WORK_DIR}/FINAL_RESULTS_MMPBSA.dat" \
    -do "${WORK_DIR}/FINAL_DECOMP_MMPBSA.dat" \
    -sp "$PARM_FILE" \
    -cp "$DRY_COMPLEX" \
    -rp "$RECEPTOR_TOP" \
    -lp "$LIGAND_TOP" \
    -y "$TRAJ_CENTERED"

if [ ! -s "${WORK_DIR}/FINAL_RESULTS_MMPBSA.dat" ]; then
    echo "ERROR: MMPBSA result file was not created."
    exit 1
fi

echo "====================================================="
echo "Pipeline completed successfully for ${SYSTEM_ID}"
echo "Results directory: ${WORK_DIR}"
echo "====================================================="
