from modeller import *
from modeller.automodel import *

env = Environ()
env.io.atom_files_directory = ['/home/gbiuser/Documents/vitoria/usp-masters/6-preparation/target/1-loopmodel']
env.io.hetatm = False
env.io.water = False

aln = Alignment(env)

# Sequence from PDB (ucr2+lr2)
mdl = Model(env, file='4X0F_nohetatm')
aln.append_model(mdl, align_codes='4X0F_nohetatm', atom_files='4X0F_nohetatm.pdb')

# Target sequence (without mutations WT) from .ali 
aln.append(file='structures.ali', align_codes='target')

# 2D alignment
aln.align2d(max_gap_length=50)

# Save
aln.write(file='structures_align2d.ali', alignment_format='PIR')
aln.write(file='structures_align2d.pap', alignment_format='PAP')  

class MyModel(LoopModel):

    def select_loop_atoms(self):
        return Selection(self.residue_range('20:A', '57:A'),
                          self.residue_range('117:A', '160:A'))

    def select_atoms(self):
        # todos os átomos entram na otimização normalmente
        return Selection(self)

    def tether_ca(self, rsr, residuos, stdev):
        for r in residuos:
            ca = self.residues[r].atoms['CA']
            rsr.add(forms.Gaussian(group=physical.xy_distance,
                                    feature=features.XCoordinate(ca),
                                    mean=ca.x, stdev=stdev))
            rsr.add(forms.Gaussian(group=physical.xy_distance,
                                    feature=features.YCoordinate(ca),
                                    mean=ca.y, stdev=stdev))
            rsr.add(forms.Gaussian(group=physical.xy_distance,
                                    feature=features.ZCoordinate(ca),
                                    mean=ca.z, stdev=stdev))

    def special_restraints(self, aln):
        rsr = self.restraints

        residuos_proximos = [
            '109:A', '113:A', '114:A', '354:A'
        ]

        residuos_distantes = [
            '240:A', '241:A',
            '400:A', '402:A', '410:A',
            '413:A', '414:A', '417:A',
            '418:A', '421:A', '438:A',
            '449:A', '450:A', '453:A'
        ]

        # próximos do loop: mola mais frouxa -> mais liberdade pra geometria local
        self.tether_ca(rsr, residuos_proximos, stdev=0.3)

        # distantes: mola mais rígida -> mantém posição quase intacta
        self.tether_ca(rsr, residuos_distantes, stdev=0.05)


a = MyModel(
    env,
    alnfile  = 'structures_align2d.ali',
    knowns   = '4X0F_nohetatm',
    sequence = 'target'
    
)
 
a.starting_model = 1
a.ending_model   = 1

a.loop.starting_model = 1
a.loop.ending_model   = 50             
a.loop.md_level       = refine.very_slow 

a.make()
