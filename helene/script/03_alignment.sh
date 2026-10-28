#!/usr/bin/env bash
set -e

# Usage: ./run_phylo.sh S M L

cd ${WD:-/Users/srz223/Documents/courses/MolecularEvolutionPhylogenomics/shamonda_virus/helene/out}

for INPUT in "$@"; do

    echo "=== $INPUT ==="

    # Align
    mafft --auto --adjustdirection ${INPUT}.fasta > ${INPUT}_align.fasta

    # Tree
    iqtree3 -s ${INPUT}_align.fasta -B 1000 -T 1 --prefix ${INPUT}

done

