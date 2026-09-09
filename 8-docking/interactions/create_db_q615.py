import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

from plip.structure.preparation import PDBComplex
from plip.basic import config


# ============================================================
# Settings
# ============================================================

proteins = [
    "PDE4B_catalytic",
    "PDE4B_monomer",
    "PDE4B_dimer",
    "PDE4D_catalytic",
    "PDE4D_monomer",
    "PDE4D_dimer"
]

pdb_folder = "../pdbs"


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


general_attributes = {
    "resnr": "resnr",
    "restype": "restype",
    "reschain": "reschain",
    "resnr_lig": "resnr_l",
    "restype_lig": "restype_l",
    "reschain_lig": "reschain_l",
    "dist": "distance"
}


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


water_bridge_attributes = {
    "dist_a_w": "distance_aw",
    "dist_d_w": "distance_dw",
    "water_angle": "w_angle",
    "water_idx": "water_orig_idx",
}


halogen_attributes = {
    "acc_angle": "acc_angle"
}


saltbridge_attributes = {
    "protispos": "protispos"
}


pistacking_attributes = {
    "cent_dist": "distance",
    "angle": "angle",
    "offset": "offset",
    "type": "type",
}


hydrophobic_attributes = {
    "dist": "distance"
}


pication_attributes = {
    "cent_dist": "distance",
    "offset": "offset",
    "type": "type",
    "protcharged": "protcharged",
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


# ============================================================
# Get protein chains from PDB
# ============================================================

def get_protein_chains(pdb_path):
    """
    Get chain IDs from ATOM records.
    """

    chains = set()

    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                chain = line[21].strip()

                if chain:
                    chains.add(chain)

    return sorted(chains)


# ============================================================
# Process one PDB
# ============================================================

def process_q615_pdb(pdb_file):

    basename = os.path.splitext(pdb_file)[0]

    # --------------------------------------------------------
    # Identify protein and ligand from filename
    # --------------------------------------------------------

    protein = None
    ligand = None

    for p in proteins:
        if basename.startswith(p + "_"):
            protein = p
            ligand = basename[len(p) + 1:]
            break

    if protein is None:
        return []

    pdb_path = os.path.join(pdb_folder, pdb_file)

    rows = []

    # --------------------------------------------------------
    # Get protein chains
    # --------------------------------------------------------

    chains = get_protein_chains(pdb_path)

    # --------------------------------------------------------
    # Analyze each chain in INTRA mode
    # --------------------------------------------------------

    for chain in chains:

        print(
            f"Analyzing {pdb_file}, chain {chain}, "
            f"intra-chain interactions..."
        )

        # ----------------------------------------------------
        # Tell PLIP to analyze intra-chain interactions
        # ----------------------------------------------------

        config.INTRA = chain

        # Make sure other special modes are disabled
        config.CHAINS = None
        config.REGIONS = None
        config.PEPTIDES = []

        mol = PDBComplex()
        mol.load_pdb(pdb_path)
        mol.analyze()

        # ----------------------------------------------------
        # Extract interactions
        # ----------------------------------------------------

        for bsid, inter in mol.interaction_sets.items():

            for interaction_type in interaction_types:

                if not hasattr(inter, interaction_type):
                    continue

                interactions = getattr(inter, interaction_type)

                specific_attributes = interaction_attributes.get(
                    interaction_type, {}
                )

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

                    # ----------------------------------------
                    # Extract PLIP attributes
                    # ----------------------------------------

                    for attr, raw_attr in attributes.items():
                        row[attr] = getattr(i, raw_attr, None)

                    # ----------------------------------------
                    # Check both sides for Q615
                    # ----------------------------------------

                    q615_is_protein_side = (
                        row.get("restype") == "GLN"
                        and row.get("resnr") == 615
                    )

                    q615_is_ligand_side = (
                        row.get("restype_lig") == "GLN"
                        and row.get("resnr_lig") == 615
                    )

                    if not (
                        q615_is_protein_side
                        or q615_is_ligand_side
                    ):
                        continue

                    # ----------------------------------------
                    # Determine Q615 and partner
                    # ----------------------------------------

                    if q615_is_protein_side:

                        q615_chain = row["reschain"]

                        partner_type = row["restype_lig"]
                        partner_number = row["resnr_lig"]
                        partner_chain = row["reschain_lig"]

                    else:

                        q615_chain = row["reschain_lig"]

                        partner_type = row["restype"]
                        partner_number = row["resnr"]
                        partner_chain = row["reschain"]

                    # ----------------------------------------
                    # Make absolutely sure this is intra-chain
                    # ----------------------------------------

                    if q615_chain != partner_chain:
                        continue

                    # ----------------------------------------
                    # Don't count Q615 interacting with itself
                    # ----------------------------------------

                    if (
                        partner_type == "GLN"
                        and partner_number == 615
                        and partner_chain == q615_chain
                    ):
                        continue

                    # ----------------------------------------
                    # Store residue IDs
                    # ----------------------------------------

                    row["q615_residue_id"] = (
                        f"GLN615:{q615_chain}"
                    )

                    row["partner_residue_id"] = (
                        f"{partner_type}{partner_number}:{partner_chain}"
                    )

                    rows.append(row)

    return rows


# ============================================================
# Main
# ============================================================

def main_q615():

    pdb_files = sorted(
        f
        for f in os.listdir(pdb_folder)
        if f.endswith(".pdb")
    )

    all_rows = []

    # --------------------------------------------------------
    # Parallel processing
    # --------------------------------------------------------

    with ProcessPoolExecutor(max_workers=12) as executor:

        for rows in executor.map(
            process_q615_pdb,
            pdb_files
        ):

            all_rows.extend(rows)

    # --------------------------------------------------------
    # Create DataFrame
    # --------------------------------------------------------

    df = pd.DataFrame(all_rows)

    # --------------------------------------------------------
    # Remove exact duplicates
    # --------------------------------------------------------

    if len(df) > 0:

        df = df.drop_duplicates()

        preferred_columns = [
            "protein",
            "ligand",
            "q615_residue_id",
            "partner_residue_id",
            "interaction_type",
            "dist",
            "resnr",
            "restype",
            "reschain",
            "resnr_lig",
            "restype_lig",
            "reschain_lig",
        ]

        existing_columns = [
            col
            for col in preferred_columns
            if col in df.columns
        ]

        remaining_columns = [
            col
            for col in df.columns
            if col not in existing_columns
        ]

        df = df[existing_columns + remaining_columns]

    # --------------------------------------------------------
    # Save
    # --------------------------------------------------------

    output_file = "Q615_intra_interactions.csv"

    df.to_csv(
        output_file,
        index=False
    )

    print(
        f"Done. {len(df)} Q615 intra-chain interactions "
        f"written to {output_file}"
    )


if __name__ == "__main__":
    main_q615()