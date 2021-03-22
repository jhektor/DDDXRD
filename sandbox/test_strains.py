import dddxrd.utils.crystallography as cry
import dddxrd.visualization.polefigures as pf
import dddxrd.visualization.Grainmap as gm
import dddxrd.utils.strain as strain
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import sys

mapfile = "dddxrd/tests/test.map"
#mapfile = "/Users/al8720/Box/projects/castIron/DCI_1_01_merged.map"
grains = gm.Grainmap(mapfile)
grains.plot_3d_map(coloring=grains.I1,cmap=cm.RdBu_r)

d0 = strain.average_cell(grains.grains,make_cubic=True)
E,e = strain.calc_strain(grains.grains[0].ubi,d0)
print(E)
Egm = grains.Green_strain[0]
print(strain.voight_to_matrix(Egm))

plt.figure()
plt.hist(grains.I1)

plt.show()

