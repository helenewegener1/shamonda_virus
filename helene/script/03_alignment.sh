#!/usr/bin/env bash
set -e

cd /Users/srz223/Documents/courses/MolecularEvolutionPhylogenomics/shamonda_virus/helene/out

for INPUT in S M L; do

    echo "=== $INPUT ==="

    # Align
    mafft --auto --adjustdirection ${INPUT}.fasta > ${INPUT}_align.fasta

    # Tree
    iqtree3 -s ${INPUT}_align.fasta -B 1000 -T AUTO --prefix ${INPUT}

done

