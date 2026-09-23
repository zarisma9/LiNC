#!/usr/bin/env bash

set -euo pipefail

RESNAME="${1:-${RESNAME:-}}"
MOL2_FILE="alpha_conf.mol2"

if [ -z "$RESNAME" ]; then
    echo "Usage:"
    echo "  $0 RESNAME"
    echo "or"
    echo "  RESNAME=D4K $0"
    exit 1
fi

REQUIRED_COMMANDS=(espgen antechamber respgen resp)

for cmd in "${REQUIRED_COMMANDS[@]}"; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: $cmd is not available in PATH."
        exit 1
    fi
done

echo "========================================"
echo " RESP Preparation and Calculation"
echo " Target: $RESNAME"
echo "========================================"

# -----------------------------------------------------------
# STEP 0. Check required input files
# -----------------------------------------------------------

REQUIRED_FILES=(
    "alpha_conf.mol2"
    "beta_conf.mol2"
    "${RESNAME}a.out"
    "${RESNAME}b.out"
)

for file in "${REQUIRED_FILES[@]}"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: $file not found."
        exit 1
    fi
done

NATOMS=$(sed -n '/@<TRIPOS>MOLECULE/{n;n;p;q}' "$MOL2_FILE" | awk '{print $1}')

if [ -z "$NATOMS" ] || ! [[ "$NATOMS" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Failed to determine atom count from $MOL2_FILE."
    exit 1
fi

echo ">>> Detected Atom Count: $NATOMS"

# -----------------------------------------------------------
# STEP 1. Generate and merge ESP files
# -----------------------------------------------------------

echo ">>> Generating ESP files..."

espgen -i "${RESNAME}a.out" -o "${RESNAME}a.esp"
espgen -i "${RESNAME}b.out" -o "${RESNAME}b.esp"

cat "${RESNAME}a.esp" "${RESNAME}b.esp" > "${RESNAME}.esp"

echo ">>> ESP Merge Complete: ${RESNAME}.esp"

# -----------------------------------------------------------
# STEP 2. Generate AC and RESPIN files
# -----------------------------------------------------------

for type in a b; do
    echo ">>> Processing Conformer: $type"

    if [ "$type" = "a" ]; then
        INPUT_MOL2="alpha_conf.mol2"
    else
        INPUT_MOL2="beta_conf.mol2"
    fi

    antechamber \
        -i "$INPUT_MOL2" \
        -fi mol2 \
        -o "${RESNAME}${type}.ac" \
        -fo ac \
        -c bcc \
        -nc 0

    respgen \
        -i "${RESNAME}${type}.ac" \
        -o "${RESNAME}-step1${type}.respin" \
        -f resp1

    respgen \
        -i "${RESNAME}${type}.ac" \
        -o "${RESNAME}-step2${type}.respin" \
        -f resp2
done

# -----------------------------------------------------------
# STEP 3. Merge and format RESPIN files
# -----------------------------------------------------------

echo ">>> Formatting RESPIN files..."

for step in 1 2; do
    FILE_A="${RESNAME}-step${step}a.respin"
    FILE_B="${RESNAME}-step${step}b.respin"
    MERGED="${RESNAME}-step${step}.respin"

    if [ ! -f "$FILE_A" ] || [ ! -f "$FILE_B" ]; then
        echo "ERROR: RESPIN input files for step $step are missing."
        exit 1
    fi

    if [ "$step" -eq 1 ]; then
        sed \
            -e 's/nmol = 1/nmol = 2/g' \
            -e '/iqopt/d' \
            "$FILE_A" > "$MERGED"
    else
        sed \
            -e 's/nmol = 1/nmol = 2/g' \
            -e '/nmol = 2,/a\ iqopt = 2,' \
            "$FILE_A" > "$MERGED"
    fi

    sed '1,/&end/d' "$FILE_B" >> "$MERGED"

    echo "" >> "$MERGED"

    for i in $(seq 1 "$NATOMS"); do
        printf "%5d\n" 2 >> "$MERGED"
        printf "%5d%5d%5d%5d\n" 1 "$i" 2 "$i" >> "$MERGED"
    done

    echo "" >> "$MERGED"

    echo ">>> Prepared Input: $MERGED"
done

# -----------------------------------------------------------
# STEP 4. Run RESP calculations
# -----------------------------------------------------------

echo "----------------------------------------"
echo " Starting RESP Calculation"
echo "----------------------------------------"

echo ">>> Running Step 1..."

resp -O \
    -i "${RESNAME}-step1.respin" \
    -o "${RESNAME}-step1.respout" \
    -e "${RESNAME}.esp" \
    -t qout-stage1

if [ ! -s "qout-stage1" ]; then
    echo "ERROR: RESP Step 1 failed."
    exit 1
fi

cp qout-stage1 qin

echo ">>> Running Step 2..."

resp -O \
    -i "${RESNAME}-step2.respin" \
    -o "${RESNAME}-step2.respout" \
    -e "${RESNAME}.esp" \
    -q qin \
    -t qout-stage2

if [ ! -s "qout-stage2" ]; then
    echo "ERROR: RESP Step 2 failed."
    exit 1
fi

# -----------------------------------------------------------
# STEP 5. Generate final AC file
# -----------------------------------------------------------

echo ">>> Finalizing..."

mv qout-stage2 "${RESNAME}.crg"

antechamber \
    -fi mol2 \
    -fo ac \
    -i "$MOL2_FILE" \
    -o "${RESNAME}_final.ac" \
    -c rc \
    -cf "${RESNAME}.crg" \
    -at amber

echo "========================================"

if [ -f "${RESNAME}_final.ac" ]; then
    echo "SUCCESS: ${RESNAME}_final.ac created."
else
    echo "ERROR: Final AC file creation failed."
    exit 1
fi

echo "========================================"

rm -f qin qout-stage1
