#!/usr/bin/env bash

set -euo pipefail

SYS_NAME="${1:-${SYS_NAME:-}}"

if [ -z "$SYS_NAME" ]; then
    echo "Usage:"
    echo "  $0 SYS_NAME"
    echo "or"
    echo "  SYS_NAME=D8K_seed1 $0"
    exit 1
fi

PRMTOP="${SYS_NAME}.prmtop"
HMR_PRMTOP="${SYS_NAME}.hmass.prmtop"
INPCRD="${SYS_NAME}.inpcrd"

MPI_NP="${MPI_NP:-4}"
DO_CPU="${DO_CPU:-mpirun -np ${MPI_NP} pmemd.MPI}"
DO_GPU="${DO_GPU:-pmemd.cuda}"

INPUT_FILES=(
    "01_Min1.in"
    "02_Min2.in"
    "03_Heat_NVT.in"
    "04_Equil_NPT_1.in"
    "05_Equil_NPT_2.in"
    "06_Equil_NPT_3.in"
    "07_Prod_HMR_4fs.in"
)

echo "========================================"
echo " AMBER MD Workflow"
echo "========================================"
echo "System name   : $SYS_NAME"
echo "PRMTOP        : $PRMTOP"
echo "HMR PRMTOP    : $HMR_PRMTOP"
echo "INPCRD        : $INPCRD"
echo "CPU command   : $DO_CPU"
echo "GPU command   : $DO_GPU"
echo "Node          : ${SLURMD_NODENAME:-unknown}"
echo "Start time    : $(date)"
echo "========================================"

# -----------------------------------------------------------
# Check required files
# -----------------------------------------------------------

for file in "$PRMTOP" "$HMR_PRMTOP" "$INPCRD"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Required file not found: $file"
        exit 1
    fi
done

for file in "${INPUT_FILES[@]}"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: MD input file not found: $file"
        exit 1
    fi
done

# -----------------------------------------------------------
# Helper function
# -----------------------------------------------------------

check_restart() {
    local file="$1"
    local step="$2"

    if [ ! -s "$file" ]; then
        echo "ERROR: $step failed. Restart file was not created: $file"
        exit 1
    fi
}

# -----------------------------------------------------------
# Step 1. Minimization 1
# CPU / standard topology
# -----------------------------------------------------------

echo ">>> Step 1: Minimization 1"

$DO_CPU \
    -O \
    -i 01_Min1.in \
    -o 01_Min1.out \
    -p "$PRMTOP" \
    -c "$INPCRD" \
    -r 01_Min1.rst \
    -ref "$INPCRD"

check_restart "01_Min1.rst" "Step 1"

# -----------------------------------------------------------
# Step 2. Minimization 2
# CPU / standard topology
# -----------------------------------------------------------

echo ">>> Step 2: Minimization 2"

$DO_CPU \
    -O \
    -i 02_Min2.in \
    -o 02_Min2.out \
    -p "$PRMTOP" \
    -c 01_Min1.rst \
    -r 02_Min2.rst \
    -ref "$INPCRD"

check_restart "02_Min2.rst" "Step 2"

# -----------------------------------------------------------
# Step 3. Heating
# GPU / standard topology
# -----------------------------------------------------------

echo ">>> Step 3: Heating NVT"

$DO_GPU \
    -O \
    -i 03_Heat_NVT.in \
    -o 03_Heat_NVT.out \
    -p "$PRMTOP" \
    -c 02_Min2.rst \
    -r 03_Heat.rst \
    -x 03_Heat.nc \
    -ref 02_Min2.rst

check_restart "03_Heat.rst" "Step 3"

# -----------------------------------------------------------
# Step 4. Equilibration 1
# GPU / standard topology
# -----------------------------------------------------------

echo ">>> Step 4: Equilibration NPT 1"

$DO_GPU \
    -O \
    -i 04_Equil_NPT_1.in \
    -o 04_Equil_NPT_1.out \
    -p "$PRMTOP" \
    -c 03_Heat.rst \
    -r 04_Equil1.rst \
    -x 04_Equil1.nc \
    -ref 03_Heat.rst

check_restart "04_Equil1.rst" "Step 4"

# -----------------------------------------------------------
# Step 5. Equilibration 2
# GPU / standard topology
# -----------------------------------------------------------

echo ">>> Step 5: Equilibration NPT 2"

$DO_GPU \
    -O \
    -i 05_Equil_NPT_2.in \
    -o 05_Equil_NPT_2.out \
    -p "$PRMTOP" \
    -c 04_Equil1.rst \
    -r 05_Equil2.rst \
    -x 05_Equil2.nc \
    -ref 04_Equil1.rst

check_restart "05_Equil2.rst" "Step 5"

# -----------------------------------------------------------
# Step 6. Equilibration 3
# GPU / standard topology / unrestrained
# -----------------------------------------------------------

echo ">>> Step 6: Equilibration NPT 3"

$DO_GPU \
    -O \
    -i 06_Equil_NPT_3.in \
    -o 06_Equil_NPT_3.out \
    -p "$PRMTOP" \
    -c 05_Equil2.rst \
    -r 06_Equil3.rst \
    -x 06_Equil3.nc

check_restart "06_Equil3.rst" "Step 6"

# -----------------------------------------------------------
# Step 7. Production
# GPU / HMR topology / 4 fs
# -----------------------------------------------------------

echo ">>> Step 7: Production HMR 4 fs"

$DO_GPU \
    -O \
    -i 07_Prod_HMR_4fs.in \
    -o 07_Prod_HMR_4fs.out \
    -p "$HMR_PRMTOP" \
    -c 06_Equil3.rst \
    -r 07_Prod_HMR.rst \
    -x 07_Prod_HMR.nc

check_restart "07_Prod_HMR.rst" "Step 7"

echo "========================================"
echo "Simulation completed successfully."
echo "End time: $(date)"
echo "========================================"
