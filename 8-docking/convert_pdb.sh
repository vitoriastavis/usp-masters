#!/bin/bash

mkdir -p pdbs

find complexes -name "*.pdbqt" -print0 |
    xargs -0 -P 8 -I {} bash -c '
        f="$1"
        obabel "$f" -O "pdbs/$(basename "${f%.pdbqt}.pdb")"
    ' _ {}