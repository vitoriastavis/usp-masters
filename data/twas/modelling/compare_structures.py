from modeller import *
import os

# Get all .pdb files in the current directory
pdb_files = [f for f in os.listdir('.') if f.endswith('.pdb')]

# Process each filename
pdb_ids = []
for filename in pdb_files:
    # Get the final part after splitting by a delimiter, e.g., underscore or dot
    final_part = filename.rsplit('.', 1)[0].split('_')[-1]
    pdb_ids.append((final_part, 'A'))

# print(pdb_ids)

env = Environ()
aln = Alignment(env)
for (pdb, chain) in pdb_ids:
    m = Model(env, file=pdb, model_segment=('FIRST:'+chain, 'LAST:'+chain))
    aln.append_model(m, atom_files=pdb, align_codes=pdb+chain)
aln.malign()
aln.malign3d()
aln.compare_structures()
aln.id_table(matrix_file='family.mat')
env.dendrogram(matrix_file='family.mat', cluster_cut=-1.0)