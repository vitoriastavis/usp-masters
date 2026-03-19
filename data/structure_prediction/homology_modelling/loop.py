from modeller import *
from modeller.automodel import *

class MyLoop(LoopModel):
    def select_loop_atoms(self):
        return Selection(
            self.residue_range('166:A', '187:A'), #ucr1
            self.residue_range('188:A', '210:A'), #lr1
            self.residue_range('211:A', '282:A'), #ucr2
            self.residue_range('283:A', '327:A'), #lr2
        )

env = Environ()
env.rand_seed = -1

a = MyLoop(env,
           alnfile='pde4d-pde4b.ali',
           knowns=('pde4d_swissmodel', 'pde4b_3g45'),
           sequence='pde4b_target',
           inimodel='pde4b_target.B99990008.pdb',
           assess_methods=assess.DOPE)

a.loop.starting_model = 1
a.loop.ending_model = 5
a.loop.md_level = refine.slow

a.make()