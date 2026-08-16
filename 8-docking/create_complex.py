from pathlib import Path

outdir = Path("complexes")
outdir.mkdir(exist_ok=True)

for receptor in Path(".").glob("PDE4B_dimer*.pdbqt"):

    target = receptor.stem
    target_dir = Path("results_docking") / target

    if not target_dir.is_dir():
        continue

    for ligand_dir in target_dir.iterdir():

        if not ligand_dir.is_dir():
            continue

        ligand = ligand_dir.name
        ligand_file = ligand_dir / f"{ligand}_rep1_out.pdbqt"

        if not ligand_file.exists():
            continue

        outfile = outdir / f"{target}_{ligand}.pdbqt"

        with open(outfile, "w") as out:

            # receptor
            with open(receptor) as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")):
                        out.write(line)

            out.write("TER\n")

            with open(receptor) as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")):
                        resnum = int(line[22:26])
                        new_resnum = resnum - 165

                        line = f"{line[:22]}{new_resnum:4d}{line[26:]}"
                        out.write(line)
            out.write("TER\n")

            # primeira pose do ligante
            inside_model1 = False

            with open(ligand_file) as f:
                for line in f:

                    if line.startswith("MODEL 1"):
                        inside_model1 = True
                        continue

                    if inside_model1 and line.startswith("ENDMDL"):
                        break

                    if inside_model1 and line.startswith(("ATOM", "HETATM")):
                        out.write(line)

            out.write("END\n")
