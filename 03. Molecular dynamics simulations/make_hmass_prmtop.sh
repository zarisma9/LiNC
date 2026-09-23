#!/usr/bin/env bash

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage:"
    echo "  $0 INPUT_PRMTOP OUTPUT_PRMTOP"
    echo
    echo "Example:"
    echo "  $0 M0P_seed1.prmtop M0P_seed1.hmass.prmtop"
    exit 1
fi

INPUT_PRMTOP="$1"
OUTPUT_PRMTOP="$2"
PARMED_INPUT="hmass_parmed.in"

if [ ! -f "$INPUT_PRMTOP" ]; then
    echo "ERROR: Input topology file not found: $INPUT_PRMTOP"
    exit 1
fi

if ! command -v parmed >/dev/null 2>&1; then
    echo "ERROR: parmed is not available in PATH."
    exit 1
fi

cat > "$PARMED_INPUT" <<EOF
parm $INPUT_PRMTOP
HMassRepartition
outparm $OUTPUT_PRMTOP
quit
EOF

echo "========================================"
echo " HMR Topology Preparation"
echo "========================================"
echo "Input PRMTOP  : $INPUT_PRMTOP"
echo "Output PRMTOP : $OUTPUT_PRMTOP"
echo "ParmEd input  : $PARMED_INPUT"
echo "========================================"

echo ">>> Running ParmEd..."

parmed -i "$PARMED_INPUT"

if [ ! -s "$OUTPUT_PRMTOP" ]; then
    echo "ERROR: ParmEd failed to create $OUTPUT_PRMTOP"
    exit 1
fi

echo "========================================"
echo "SUCCESS"
echo "HMR topology created: $OUTPUT_PRMTOP"
echo "========================================"
