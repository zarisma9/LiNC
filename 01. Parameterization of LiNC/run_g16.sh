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

if ! command -v g16 >/dev/null 2>&1; then
    echo "Error: g16 is not available in PATH."
    exit 1
fi

for type in a b; do
    INPUT="${RESNAME}${type}.gau"
    OUTPUT="${RESNAME}${type}.out"

    if [ ! -f "$INPUT" ]; then
        echo "Error: $INPUT was not found."
        exit 1
    fi

    echo ">>> Running: $INPUT -> $OUTPUT"
    g16 < "$INPUT" > "$OUTPUT"
    echo ">>> Finished: $OUTPUT"
done

echo ">>> All Gaussian jobs finished."
