import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from plip.structure.preparation import PDBComplex

# Protein names (longest first to avoid partial matches)
proteins = [
    "PDE4B_catalytic",
    "PDE4B_monomer",
    "PDE4B_dimer",
    "PDE4D_catalytic",
    "PDE4D_monomer",
    "PDE4D_dimer"
]

interaction_types = [
    "hydrophobic_contacts",
    "hbonds_ldon",
    "hbonds_pdon",
    "water_bridges",
    "saltbridge_lneg",
    "saltbridge_pneg",
    "pistacking",
    "pication_laro",
    "pication_paro",
    "halogen_bonds",
    "metal_complexes",
]

pdb_folder = "../pdbs"
log_file = "protein_interactions.log"


def process_pdb(pdb_file):
    """Analyze one PDB file and return interaction rows."""

    basename = os.path.splitext(pdb_file)[0]

    # Identify protein and ligand from filename
    protein = None
    ligand = None

    for p in proteins:
        if basename.startswith(p + "_"):
            protein = p
            ligand = basename[len(p) + 1:]
            break

    if protein is None:
        return None, []

    pdb_path = os.path.join(pdb_folder, pdb_file)

    rows = []

    mol = PDBComplex()
    mol.load_pdb(pdb_path)
    mol.analyze()

    for bsid, inter in mol.interaction_sets.items():

        for interaction_type in interaction_types:

            if not hasattr(inter, interaction_type):
                continue

            interactions = getattr(inter, interaction_type)

            for i in interactions:

                rows.append({
                    "protein": protein,
                    "ligand": ligand,
                    "binding_site": bsid,
                    "interaction_type": interaction_type,
                    "residue_number": getattr(i, "resnr", None),
                    "residue_name": getattr(i, "restype", None),
                    "residue_chain": getattr(i, "reschain", None),
                })

    return protein, rows


def main():

    pdb_files = sorted(
        f for f in os.listdir(pdb_folder)
        if f.endswith(".pdb")
    )

    # Count how many PDBs belong to each protein
    total_per_protein = {protein: 0 for protein in proteins}

    for pdb_file in pdb_files:
        basename = os.path.splitext(pdb_file)[0]

        for protein in proteins:
            if basename.startswith(protein + "_"):
                total_per_protein[protein] += 1
                break

    # Track completed PDBs
    completed_per_protein = {protein: 0 for protein in proteins}

    # Start a new log
    with open(log_file, "w") as log:
        log.write("PLIP interaction analysis\n")
        log.write("=========================\n")

    all_rows = []

    with ProcessPoolExecutor(max_workers=10) as executor:

        for protein, rows in executor.map(process_pdb, pdb_files):

            if protein is None:
                continue

            all_rows.extend(rows)

            completed_per_protein[protein] += 1

            # Protein is finished when all its PDBs are done
            if completed_per_protein[protein] == total_per_protein[protein]:

                with open(log_file, "a") as log:
                    log.write(
                        f"Finished: {protein} "
                        f"({completed_per_protein[protein]}/"
                        f"{total_per_protein[protein]} PDBs)\n"
                    ) 

    # Convert to DataFrame
    df = pd.DataFrame(all_rows)

    # Create residue ID
    df["residue_id"] = (
        df["residue_name"].astype(str) +
        df["residue_number"].astype("Int64").astype(str) +
        ":" +
        df["residue_chain"].astype(str)
    )

    df.to_csv("protein_interactions.csv", index=False)

    with open(log_file, "a") as log:
        log.write(f"\nDone. {len(df)} interactions written to protein_interactions.csv\n")

  
if __name__ == "__main__":
    main()