#!/bin/bash

for f in complexes/*.pdbqt; do
    obabel "$f" -O "pdbs/$(basename "${f%.pdbqt}.pdb")" &

    while [ "$(jobs -rp | wc -l)" -ge 12 ]; do
        wait -n
    done
done

wait
