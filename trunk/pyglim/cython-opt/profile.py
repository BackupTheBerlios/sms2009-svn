# A simple profiling wrapper for testpyglim

import testpyglim as tp
import cProfile
import pstats

fn = "ism.profile"

cProfile.run('tp.runISM()',fn)
p=pstats.Stats(fn)
p.sort_stats('cumulative').print_stats(25)

