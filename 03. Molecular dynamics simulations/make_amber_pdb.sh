#!/usr/bin/env bash

set -euo pipefail

if [ "$#" -ne 7 ]; then
    echo "Usage:"
    echo "  $0 INPUT_CIF CONVERTED_PDB MODIFIED_PDB CHAIN RESNUM RESNAME FINAL_PDB"
    echo
    echo "Example:"
    echo "  $0 pred.model_idx_0.cif D4K_seed1.pdb D4K_final.pdb A 20 D4K D4K_amber.pdb"
    exit 1
fi

INPUT_CIF="$1"
CONVERTED_PDB="$2"
MODIFIED_PDB="$3"
CHAIN="$4"
RESNUM="$5"
RESNAME="$6"
FINAL_PDB="$7"

MODIFY_SCRIPT="modify_pdb.py"
LIGAND_NAME="LIG"

if [ ! -f "$INPUT_CIF" ]; then
    echo "ERROR: Input CIF file not found: $INPUT_CIF"
    exit 1
fi

if [ ! -f "$MODIFY_SCRIPT" ]; then
    echo "ERROR: $MODIFY_SCRIPT not found in the current directory."
    exit 1
fi

if ! command -v obabel >/dev/null 2>&1; then
    echo "ERROR: obabel is not available in PATH."
    exit 1
fi

if ! command -v pdb4amber >/dev/null 2>&1; then
    echo "ERROR: pdb4amber is not available in PATH."
    exit 1
fi

echo "========================================"
echo " Amber PDB Preparation"
echo "========================================"
echo "Input CIF      : $INPUT_CIF"
echo "Converted PDB  : $CONVERTED_PDB"
echo "Modified PDB   : $MODIFIED_PDB"
echo "Target chain   : $CHAIN"
echo "Target residue : $RESNUM"
echo "Residue name   : $RESNAME"
echo "Final PDB      : $FINAL_PDB"
echo "Ligand name    : $LIGAND_NAME"
echo "========================================"

echo ">>> Converting CIF to PDB..."
obabel -i cif "$INPUT_CIF" -o pdb -O "$CONVERTED_PDB"

if [ ! -s "$CONVERTED_PDB" ]; then
    echo "ERROR: Open Babel failed to create $CONVERTED_PDB"
    exit 1
fi

echo ">>> Modifying residue information..."
python3 "$MODIFY_SCRIPT" \
    -i "$CONVERTED_PDB" \
    -o "$MODIFIED_PDB" \
    -c "$CHAIN" \
    -r "$RESNUM" \
    -n "$RESNAME" \
    -l "$LIGAND_NAME"

if [ ! -s "$MODIFIED_PDB" ]; then
    echo "ERROR: modify_pdb.py failed to create $MODIFIED_PDB"
    exit 1
fi

echo ">>> Running pdb4amber..."
pdb4amber "$MODIFIED_PDB" > "$FINAL_PDB"

if [ ! -s "$FINAL_PDB" ]; then
    echo "ERROR: pdb4amber failed to create $FINAL_PDB"
    exit 1
fi

echo "========================================"
echo "SUCCESS"
echo "Final Amber PDB: $FINAL_PDB"
echo "========================================"
