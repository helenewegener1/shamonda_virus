#!/usr/bin/env bash
set -e

# Cut an alignment into non-overlapping windows and infer a tree per window.
#
# Usage: ./window_cut.sh M_align.fasta [WINDOW]
#
# Writes files_for_mast/M_align_w0001_0500.fasta + .treefile, etc.

WD=${WD:-/Users/srz223/Documents/courses/MolecularEvolutionPhylogenomics/shamonda_virus/helene/out}
cd "$WD"

ALN=$1
WINDOW=${2:-500}

OUTDIR=files_for_mast
mkdir -p "$OUTDIR"

BASE=$(basename "$ALN" .fasta)

# Alignment length = length of the first sequence (all equal in an alignment)
LEN=$(awk '/^>/ {if (n) exit} !/^>/ {n += length($0)} END {print n}' "$ALN")

echo "Alignment: $ALN"
echo "Length:    $LEN columns"
echo "Window:    $WINDOW"
echo

# ------------------------------------------------------------------------------
# Cut into windows
# ------------------------------------------------------------------------------

START=1
while [ "$START" -le "$LEN" ]; do

    END=$((START + WINDOW - 1))

    # Skip a trailing remainder shorter than the window: too few sites for a tree
    if [ "$END" -gt "$LEN" ]; then
        echo "  skipping remainder ${START}-${LEN} ($((LEN - START + 1)) columns)"
        break
    fi

    OUT=$(printf "%s/%s_w%04d_%04d.fasta" "$OUTDIR" "$BASE" "$START" "$END")

    seqconverter --informat fasta --outformat fasta \
                 --keepcols ${START}-${END} \
                 -i "$ALN" > "$OUT"

    echo "  ${START}-${END} -> $OUT"

    START=$((END + 1))
done

echo

# ------------------------------------------------------------------------------
# Tree per window
# ------------------------------------------------------------------------------

for WIN in "$OUTDIR"/${BASE}_w*.fasta; do

    echo "=== $(basename "$WIN") ==="
    iqtree3 -s "$WIN" -redo \
            > "${WIN%.fasta}.log" 2>&1

done

# Collect the window trees into one file for MAST
cat "$OUTDIR"/${BASE}_w*.treefile > "${BASE}_candidates.tre"

echo
echo "Windows:   $(ls -1 "$OUTDIR"/${BASE}_w*.fasta | wc -l | tr -d ' ')"
echo "Trees:     ${BASE}_candidates.tre"


