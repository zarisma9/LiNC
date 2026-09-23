#!/usr/bin/env bash

set -euo pipefail

RESNAME="${1:-${RESNAME:-}}"

if [ -z "$RESNAME" ]; then
    echo "Usage:"
    echo "  $0 RESNAME"
    echo "or"
    echo "  RESNAME=D4K $0"
    exit 1
fi

if ! command -v prepgen >/dev/null 2>&1; then
    echo "ERROR: prepgen is not available in PATH."
    exit 1
fi

if ! command -v parmchk2 >/dev/null 2>&1; then
    echo "ERROR: parmchk2 is not available in PATH."
    exit 1
fi

if [ -z "${AMBERHOME:-}" ]; then
    echo "ERROR: AMBERHOME is not set."
    exit 1
fi

INPUT_AC="${RESNAME}_final.ac"
MAINCHAIN_FILE="mainchain.${RESNAME}"
PREPI_FILE="${RESNAME}.prepi"
RES_FILE="${RESNAME}.res"
FRCMOD_FILE="${RESNAME}.frcmod"
PARM_FILE="${AMBERHOME}/dat/leap/parm/parm10.dat"

if [ ! -f "$INPUT_AC" ]; then
    echo "ERROR: $INPUT_AC not found."
    exit 1
fi

if [ ! -f "$PARM_FILE" ]; then
    echo "ERROR: $PARM_FILE not found."
    exit 1
fi

echo "========================================"
echo " Amber Parameter Generation"
echo " Target: $RESNAME"
echo "========================================"

# -----------------------------------------------------------
# STEP 1. Create mainchain definition file
# -----------------------------------------------------------

cat > "$MAINCHAIN_FILE" <<EOF
HEAD_NAME N
TAIL_NAME C
MAIN_CHAIN CA
OMIT_NAME C101
OMIT_NAME C102
OMIT_NAME HAC1
OMIT_NAME HAC2
OMIT_NAME HAC3
OMIT_NAME O100
OMIT_NAME N100
OMIT_NAME HNM1
OMIT_NAME HNM2
OMIT_NAME HNM3
OMIT_NAME HNM4
OMIT_NAME C103
PRE_HEAD_TYPE C
POST_TAIL_TYPE N
CHARGE 0.0
EOF

echo ">>> Created: $MAINCHAIN_FILE"

# -----------------------------------------------------------
# STEP 2. Run prepgen
# -----------------------------------------------------------

echo ">>> Running prepgen..."

prepgen \
    -i "$INPUT_AC" \
    -o "$PREPI_FILE" \
    -f prepi \
    -m "$MAINCHAIN_FILE" \
    -rn "$RESNAME" \
    -fr "$RES_FILE"

if [ ! -s "$PREPI_FILE" ]; then
    echo "ERROR: prepgen failed to create $PREPI_FILE."
    exit 1
fi

echo ">>> Created: $PREPI_FILE"

# -----------------------------------------------------------
# STEP 3. Run parmchk2
# -----------------------------------------------------------

echo ">>> Running parmchk2..."

parmchk2 \
    -i "$INPUT_AC" \
    -f ac \
    -o "$FRCMOD_FILE" \
    -a Y \
    -p "$PARM_FILE"

if [ ! -s "$FRCMOD_FILE" ]; then
    echo "ERROR: parmchk2 failed to create $FRCMOD_FILE."
    exit 1
fi

echo ">>> Created: $FRCMOD_FILE"

echo "========================================"
echo " All processes completed successfully."
echo " Result files:"
echo "   $PREPI_FILE"
echo "   $FRCMOD_FILE"
echo "   $RES_FILE"
echo "========================================"
