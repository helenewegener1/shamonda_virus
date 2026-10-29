

# Cut an alignment into non-overlapping windows and infer a tree per window.
#
# Usage: ./window_cut.sh ALIGNMENT.fasta [WINDOW] [MERGE_BELOW]
#
#   WINDOW       window size in columns            (default 500)
#   MERGE_BELOW  if the trailing remainder is shorter than this, fold it
#                into the previous window instead of standing alone
#                (default 0 = always keep the remainder as its own window)
#
# Writes files_for_mast/<base>_w0001_0500.fasta + .treefile, etc.
WD=${WD:-//Users/amaqve/Documents/Amandas_PhD/Courses/molevol_phylogenomics/project_work/shamonda_virus/helene/out}
cd "$WD"
#!/usr/bin/env bash
set -euo pipefail

# Cut an alignment into (optionally overlapping) windows and infer a tree per window.
#
# Usage: ./window_cut.sh ALIGNMENT.fasta [WINDOW] [STEP] [MIN_WIN]
#
#   WINDOW   window size in columns                     (default 500)
#   STEP     distance between window starts             (default WINDOW/2)
#              STEP  < WINDOW -> overlapping windows
#              STEP == WINDOW -> non-overlapping, back-to-back
#   MIN_WIN  drop a trailing window shorter than this   (default 0 = keep it)
#
# e.g. ./window_cut.sh M.fasta 500 250   ->  1-500, 251-750, 501-1000, ...
#
# Writes files_for_mast/<base>_w0001_0500.fasta + .treefile, etc.

ALN=${1:?usage: window_cut.sh ALIGNMENT.fasta [WINDOW] [STEP] [MIN_WIN]}
WINDOW=${2:-500}
STEP=${3:-$((WINDOW / 2))}
MIN_WIN=${4:-0}

[ "$STEP" -ge 1 ] || { echo "STEP must be >= 1" >&2; exit 1; }

[ -n "${WD:-}" ] && cd "$WD"

OUTDIR=files_for_mast
mkdir -p "$OUTDIR"

BASE=$(basename "$ALN" .fasta)

# Alignment length = length of the first sequence (all equal in an alignment)
LEN=$(awk '/^>/ {if (n) exit} !/^>/ {n += length($0)} END {print n}' "$ALN")

echo "Alignment:  $ALN"
echo "Length:     $LEN columns"
echo "Window:     $WINDOW"
echo "Step:       $STEP  ($([ "$STEP" -lt "$WINDOW" ] && echo "overlapping, $((WINDOW - STEP)) columns shared" || echo "non-overlapping"))"
echo "Min window: $([ "$MIN_WIN" -gt 0 ] && echo "$MIN_WIN columns" || echo "none (trailing remainder kept)")"
echo

# Clear this base's windows from a previous run, so a different WINDOW/STEP
# cannot leave stale files that the globs below would silently pick up.
STALE=$(find "$OUTDIR" -maxdepth 1 -name "${BASE}_w*" -type f 2>/dev/null | wc -l | tr -d ' ')
if [ "$STALE" -gt 0 ]; then
    echo "Removing $STALE stale file(s) from a previous run of this alignment"
    find "$OUTDIR" -maxdepth 1 -name "${BASE}_w*" -type f -delete
    echo
fi

# ------------------------------------------------------------------------------
# Cut into windows
# ------------------------------------------------------------------------------

START=1
while [ "$START" -le "$LEN" ]; do

    END=$((START + WINDOW - 1))
    [ "$END" -gt "$LEN" ] && END=$LEN
    NCOL=$((END - START + 1))

    # A trailing window too short to carry phylogenetic signal
    if [ "$NCOL" -lt "$MIN_WIN" ]; then
        echo "  dropping ${START}-${END} ($NCOL columns, below MIN_WIN=$MIN_WIN)"
        break
    fi

    OUT=$(printf "%s/%s_w%04d_%04d.fasta" "$OUTDIR" "$BASE" "$START" "$END")

    seqconverter --informat fasta --outformat fasta \
                 --keepcols ${START}-${END} \
                 -i "$ALN" > "$OUT"

    if [ "$NCOL" -lt "$WINDOW" ]; then
        echo "  ${START}-${END} -> $OUT  (short window, $NCOL columns)"
    else
        echo "  ${START}-${END} -> $OUT"
    fi

    # This window already reached the last column: another start would only
    # produce a stub nested inside it.
    [ "$END" -ge "$LEN" ] && break

    START=$((START + STEP))
done

echo

# ------------------------------------------------------------------------------
# Tree per window
# ------------------------------------------------------------------------------

FAILED=()

for WIN in "$OUTDIR"/"${BASE}"_w*.fasta; do

    echo "=== $(basename "$WIN") ==="

    # Don't let one bad window abort the whole run
    if ! iqtree3 -s "$WIN" -redo > "${WIN%.fasta}.log" 2>&1; then
        echo "    FAILED - see ${WIN%.fasta}.log"
        FAILED+=("$(basename "$WIN")")
    fi

done

echo

# Collect the window trees into one file for MAST
NTREE=0
: > "${BASE}_candidates.tre"
for T in "$OUTDIR"/"${BASE}"_w*.treefile; do
    [ -e "$T" ] || continue
    cat "$T" >> "${BASE}_candidates.tre"
    NTREE=$((NTREE + 1))
done

NWIN=$(find "$OUTDIR" -maxdepth 1 -name "${BASE}_w*.fasta" -type f | wc -l | tr -d ' ')

echo "Windows:    $NWIN"
echo "Trees:      $NTREE  -> ${BASE}_candidates.tre"

if [ "${#FAILED[@]}" -gt 0 ]; then
    echo
    echo "WARNING: ${#FAILED[@]} window(s) produced no tree:"
    printf '  %s\n' "${FAILED[@]}"
fi

if [ "$STEP" -lt "$WINDOW" ]; then
    echo
    echo "NOTE: windows overlap, so neighbouring trees share data and many will be"
    echo "      identical in topology. Collapse duplicates before running MAST:"
    echo "        iqtree3 -t ${BASE}_candidates.tre -rf_all"
fi