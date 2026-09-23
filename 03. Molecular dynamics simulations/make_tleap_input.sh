#!/usr/bin/env bash

set -euo pipefail

if [ "$#" -lt 7 ] || [ "$#" -gt 8 ]; then
    echo "Usage:"
    echo "  $0 AMBERPREP AMBERPARAMS PDB_FILE K_COUNT CL_COUNT PRMTOP_NAME INPCRD_NAME [TLEAP_INPUT]"
    echo
    echo "Example:"
    echo "  $0 D4K.prepi D4K.frcmod D4K_amber.pdb 88 71 D4K_seed1.prmtop D4K_seed1.inpcrd"
    echo
    echo "Optional:"
    echo "  TLEAP_INPUT defaults to tleap.in"
    exit 1
fi

AMBERPREP="$1"
AMBERPARAMS="$2"
PDB_FILE="$3"
K_COUNT="$4"
CL_COUNT="$5"
PRMTOP_NAME="$6"
INPCRD_NAME="$7"
TLEAP_INPUT="${8:-tleap.in}"

if [ ! -f "$AMBERPREP" ]; then
    echo "ERROR: Amber prep file not found: $AMBERPREP"
    exit 1
fi

if [ ! -f "$AMBERPARAMS" ]; then
    echo "ERROR: Amber parameter file not found: $AMBERPARAMS"
    exit 1
fi

if [ ! -f "$PDB_FILE" ]; then
    echo "ERROR: PDB file not found: $PDB_FILE"
    exit 1
fi

if ! [[ "$K_COUNT" =~ ^[0-9]+$ ]]; then
    echo "ERROR: K_COUNT must be a non-negative integer."
    exit 1
fi

if ! [[ "$CL_COUNT" =~ ^[0-9]+$ ]]; then
    echo "ERROR: CL_COUNT must be a non-negative integer."
    exit 1
fi

cat > "$TLEAP_INPUT" <<EOF
source leaprc.protein.ff19SB
source leaprc.water.opc
source leaprc.mimetic.ff15ipq

loadamberprep $AMBERPREP
loadamberparams $AMBERPARAMS

m = loadpdb $PDB_FILE

solvateOct m OPCBOX 12.0

addIonsRand m K+ $K_COUNT Cl- $CL_COUNT

saveAmberParm m $PRMTOP_NAME $INPCRD_NAME

quit
EOF

echo "========================================"
echo " tleap input generated"
echo "========================================"
echo "Amber prep   : $AMBERPREP"
echo "Amber params : $AMBERPARAMS"
echo "PDB file     : $PDB_FILE"
echo "K+ ions      : $K_COUNT"
echo "Cl- ions     : $CL_COUNT"
echo "PRMTOP       : $PRMTOP_NAME"
echo "INPCRD       : $INPCRD_NAME"
echo "tleap input  : $TLEAP_INPUT"
echo "========================================"
echo
echo "Run with:"
echo "  tleap -f $TLEAP_INPUT"
