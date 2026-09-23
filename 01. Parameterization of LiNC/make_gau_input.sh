#!/usr/bin/env bash

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 'SMILES_CODE' RESNAME"
    echo "Example:"
    echo "  $0 'O=C(NC)[C@H](CCCCNC(CCCCCCCCCCCCC(O)=O)=O)NC(C)=O' D4K"
    exit 1
fi

SMILES_CODE="$1"
RESNAME="$2"

if [ ! -f "conformer.py" ]; then
    echo "Error: conformer.py was not found in the current directory."
    exit 1
fi

if ! command -v antechamber >/dev/null 2>&1; then
    echo "Error: antechamber is not available in PATH."
    exit 1
fi

echo ">>> SMILES: $SMILES_CODE"
echo ">>> RESNAME: $RESNAME"

# 1. Run conformer.py and capture the log
LOG=$(python3 conformer.py "$SMILES_CODE" "$RESNAME")

echo "$LOG"

# 2. Extract the first backbone phi/psi atom IDs
PHI_IDS=$(echo "$LOG" \
    | grep "Backbone φ atoms: ids" \
    | head -n 1 \
    | sed 's/.*ids (\([^)]*\)).*/\1/' \
    | sed 's/,//g')

PSI_IDS=$(echo "$LOG" \
    | grep "Backbone ψ atoms: ids" \
    | head -n 1 \
    | sed 's/.*ids (\([^)]*\)).*/\1/' \
    | sed 's/,//g')

if [ -z "$PHI_IDS" ]; then
    echo "Error: PHI atom IDs could not be detected."
    exit 1
fi

if [ -z "$PSI_IDS" ]; then
    echo "Error: PSI atom IDs could not be detected."
    exit 1
fi

echo ">>> Detected PHI: $PHI_IDS"
echo ">>> Detected PSI: $PSI_IDS"

if [ ! -f "alpha_conf.mol2" ]; then
    echo "Error: alpha_conf.mol2 was not generated."
    exit 1
fi

if [ ! -f "beta_conf.mol2" ]; then
    echo "Error: beta_conf.mol2 was not generated."
    exit 1
fi

# 3. Generate Gaussian input files with Antechamber
antechamber -fi mol2 -fo gcrt -i alpha_conf.mol2 -o "${RESNAME}a.gau"
antechamber -fi mol2 -fo gcrt -i beta_conf.mol2 -o "${RESNAME}b.gau"

# 4. Rebuild Gaussian input files
for type in a b; do
    FILE="${RESNAME}${type}.gau"

    if [ ! -f "$FILE" ]; then
        echo "Error: $FILE was not generated."
        exit 1
    fi

    BODY=$(sed -n '/^[[:space:]]*0[[:space:]][[:space:]]*1[[:space:]]*$/,$p' "$FILE")

    if [ -z "$BODY" ]; then
        echo "Error: Could not find the charge/multiplicity line in $FILE."
        exit 1
    fi

    cat > "$FILE" <<EOF
--Link1--
%chk=${RESNAME}${type}
#HF/6-31G* SCF=Tight Pop=MK IOp(6/33=2) IOp(6/41=10) IOp(6/42=6) opt=modredundant

Title Line

$BODY

D $PHI_IDS F
D $PSI_IDS F

EOF

    echo ">>> Generated: $FILE"
done

echo ">>> Done."
echo ">>> Output files:"
echo "    ${RESNAME}a.gau"
echo "    ${RESNAME}b.gau"
