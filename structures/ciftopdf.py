import sys
from Bio import PDB

def convert_cif_to_pdb(cif_file):
    try:
        # Define o nome do arquivo de saída
        pdb_file = cif_file.replace(".cif", ".pdb")
        
        # Carrega o arquivo .cif
        parser = PDB.MMCIFParser()
        structure = parser.get_structure("structure", cif_file)
        
        # Salva como .pdb
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(pdb_file)
        
        print(f"Conversão concluída: {pdb_file}")
    except Exception as e:
        print(f"Erro ao converter {cif_file}: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python script.py arquivo.cif")
    else:
        convert_cif_to_pdb(sys.argv[1])
