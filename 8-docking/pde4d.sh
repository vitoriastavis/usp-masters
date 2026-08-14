for f in complexes/PDE4D*.pdbqt; do
    obabel "$f" -O "pdbs/$(basename "${f%.pdbqt}.pdb")" &

    while [ "$(jobs -rp | wc -l)" -ge 8 ]; do
        wait -n
    done
done

wait