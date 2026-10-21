import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from plip.structure.preparation import PDBComplex
from plip.basic import config

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
log_file = "interactions_ucr2.log"

def process_ucr2_pdb(pdb_file):
    """Analyze UCR2 <-> binding-site interactions for dimer PDBs only."""

    basename = os.path.splitext(pdb_file)[0]

    protein = None
    ligand = None

    for p in proteins:
        if "dimer" in p and basename.startswith(p + "_"):
            protein = p
            ligand = basename[len(p) + 1:]
            break

    if protein is None:
        return None, []

    pdb_path = os.path.join(pdb_folder, pdb_file)

    ucr2_residues = [226, 227, 230, 233, 234, 274]

    binding_site_residues = [
        278, 405, 476, 519, 565, 567, 575,
        579, 582, 583, 586, 614, 615, 618
    ]

    # Detect the two protein chains
    mol = PDBComplex()
    mol.load_pdb(pdb_path)

    chains = sorted(set(
        atom.OBAtom.GetResidue().GetChain()
        for atom in mol.atoms.values()
        if atom.OBAtom.GetResidue() is not None
    ))

    if len(chains) < 2:
        return protein, []

    chain_a, chain_b = chains[:2]

    # UCR2 A -> binding site B
    # UCR2 B -> binding site A
    config.REGIONS = [
        (
            {chain_a: ucr2_residues},
            {chain_b: binding_site_residues}
        ),
        (
            {chain_b: ucr2_residues},
            {chain_a: binding_site_residues}
        )
    ]

    # Re-load after setting REGIONS
    mol = PDBComplex()
    mol.load_pdb(pdb_path)
    mol.analyze()

    rows = []

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
                    "ucr2_residue_id": None,
                    "binding_site_residue_id": None,
                }

                for attr, raw_attr in attributes.items():
                    row[attr] = getattr(i, raw_attr, None)

                # Create IDs for both participating residues
                if (
                    row["restype_lig"] is not None
                    and row["resnr_lig"] is not None
                    and row["reschain_lig"] is not None
                ):
                    row["ucr2_residue_id"] = (
                        f"{row['restype_lig']}{row['resnr_lig']}:{row['reschain_lig']}"
                    )

                if (
                    row["restype"] is not None
                    and row["resnr"] is not None
                    and row["reschain"] is not None
                ):
                    row["binding_site_residue_id"] = (
                        f"{row['restype']}{row['resnr']}:{row['reschain']}"
                    )

                rows.append(row)

    return protein, rows

def main_ucr2():

    pdb_files = sorted(
        f for f in os.listdir(pdb_folder)
        if f.endswith(".pdb") and "dimer" in f
    )

    all_rows = []

    with ProcessPoolExecutor(max_workers=12) as executor:

        for protein, rows in executor.map(
            process_ucr2_pdb,
            pdb_files
        ):
            if protein is None:
                continue

            all_rows.extend(rows)

    df = pd.DataFrame(all_rows)

    df.to_csv("ucr2_interactions.csv", index=False)

    print(
        f"Done. {len(df)} UCR2 interactions written "
        f"to ucr2_interactions.csv"
    )

if __name__ == "__main__":
    main_ucr2()