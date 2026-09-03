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

# All PLIP interaction attributes to extract
# General attributes
general_attributes = {
    "resnr": "resnr",
    "restype": "restype",
    "reschain": "reschain",
    "resnr_lig": "resnr_l",
    "restype_lig": "restype_l",
    "reschain_lig": "reschain_l",
    "dist": "distance"
}

# Hydrogen bonds
hbond_attributes = {
    "sidechain": "sidechain",
    "dist_h_a": "distance_ah",
    "dist_d_a": "distance_ad",
    "don_angle": "angle",
    "protisdon": "protisdon",
    "donoridx": "d_orig_idx",
    "donortype": "dtype",
    "acceptoridx": "a_orig_idx",
    "acceptortype": "atype",
}

# Water bridges
water_bridge_attributes = {
    "dist_a_w": "distance_aw",
    "dist_d_w": "distance_dw",
    "water_angle": "w_angle",
    "water_idx": "water_orig_idx",
}

# Halogen bonds
halogen_attributes = {
    "acc_angle": "acc_angle",
}

# Salt bridges
saltbridge_attributes = {
    "protispos": "protispos",
}
    
# Pi-stacking
pistacking_attributes = {
    "cent_dist": "distance",
    "angle": "angle",
    "offset": "offset",
    "type": "type",
}

metal_attributes = {
    "metal_idx": "metal_orig_idx",
    "metal_type": "metal_type",
    "target_idx": "target_orig_idx",
    "target_type": "target_type",
    "coordination": "coordination_num",
    "location": "location",
    "rms": "rms",
    "geometry": "geometry",
    "complexnum": "complexnum",
}

# Hydrophobic contacts
hydrophobic_attributes = {
    "dist": "distance",
}

# Pi-cation
pication_attributes = {
    "cent_dist": "distance",
    "offset": "offset",
    "type": "type",
    "protcharged": "protcharged",
}

interaction_attributes = {
    "hydrophobic_contacts": hydrophobic_attributes,
    "hbonds_ldon": hbond_attributes,
    "hbonds_pdon": hbond_attributes,
    "water_bridges": water_bridge_attributes,
    "saltbridge_lneg": saltbridge_attributes,
    "saltbridge_pneg": saltbridge_attributes,
    "pistacking": pistacking_attributes,
    "pication_laro": pication_attributes,
    "pication_paro": pication_attributes,
    "halogen_bonds": halogen_attributes,
    "metal_complexes": metal_attributes,
}

pdb_folder = "../pdbs"
log_file = "protein_interactions.log"

def process_pdb(pdb_file):
    """Analyze one PDB file and return all PLIP interaction rows."""

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

            # Get interaction-specific attributes
            specific_attributes = interaction_attributes.get(
                interaction_type, {}
            )

            # Combine general and interaction-specific attributes
            attributes = {
                **general_attributes,
                **specific_attributes
            }

            for i in interactions:

                row = {
                    "protein": protein,
                    "ligand": ligand,
                    "binding_site": bsid,
                    "interaction_type": interaction_type,
                }

                # Extract general + interaction-specific attributes
                for attr, raw_attr in attributes.items():
                    row[attr] = getattr(i, raw_attr, None)

                # Salt-bridge ligand functional group
                if interaction_type == "saltbridge_lneg":
                    row["lig_group"] = i.negative.fgroup

                elif interaction_type == "saltbridge_pneg":
                    row["lig_group"] = i.positive.fgroup
                # Metal coordinates
                if interaction_type == "metal_complexes":
                    row["metalcoo"] = tuple(float(x) for x in i.metal.coords)
                    row["targetcoo"] = tuple(float(x) for x in i.target.atom.coords)

                rows.append(row)

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

    # Parallel processing
    with ProcessPoolExecutor(max_workers=12) as executor:

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
        df["restype"].astype(str) +
        df["resnr"].astype("Int64").astype(str) +
        ":" +
        df["reschain"].astype(str)
    )

    # Save
    df.to_csv("protein_interactions.csv", index=False)

    with open(log_file, "a") as log:
        log.write(
            f"\nDone. {len(df)} interactions written "
            f"to protein_interactions.csv\n"
        )


if __name__ == "__main__":
    main()