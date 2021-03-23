import dddxrd.visualization.Grainmap as gm
import dddxrd.utils.strain as strain
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

mapfile = "dddxrd/tests/test.map"
#mapfile = "/Users/al8720/Box/projects/castIron/DCI_1_01_merged.map"
grains = gm.Grainmap(mapfile)
#grains.plot_3d_map(coloring=grains.I1,cmap=cm.RdBu_r)

d0 = strain.average_cell(grains.grains,make_cubic=True)
fig,axes = plt.subplots(2,1,sharex=True)
ax = axes.ravel()
epsilon = np.linspace(0,0.1)
for eps in epsilon:
    d0[0] = (1-eps)*d0[0]
    d0[1] = (1-0.3*eps)*d0[1]
    d0[2] = (1-0.3*eps)*d0[2]
    E,e = strain.calc_strain(grains.grains[0].ubi,d0,ml=1,mr=-1)
    ee = strain.calc_strain(grains.grains[0].ubi,d0,finite_strain=False)
    ax[0].plot(eps,E[0,0],'or')
    ax[0].plot(eps,E[1,1],'sr')
    ax[0].plot(eps,E[2,2],'+r')
    ax[0].plot(eps,e[0,0],'ob')
    ax[0].plot(eps,e[1,1],'sb')
    ax[0].plot(eps,e[2,2],'+b')

    ax[1].plot(eps,E[0,1],'or')
    ax[1].plot(eps,E[0,2],'sr')
    ax[1].plot(eps,E[1,2],'+r')
    ax[1].plot(eps,e[0,1],'ob')
    ax[1].plot(eps,e[0,2],'sb')
    ax[1].plot(eps,e[1,2],'+b')
labels_normal = ['E_aa','E_bb','E_cc','e_xx','e_yy','e_zz']
labels_shear = ['E_ab','E_ac','E_bc','e_xy','e_xz','e_yz']
ax[0].legend(labels_normal)
ax[1].legend(labels_shear)
ax[1].set_xlabel('Change in lattice d0')
ax[1].set_ylabel('Strain component')
ax[0].set_ylabel('Strain component')
Egm = grains.Green_strain[0]
# print(strain.voight_to_matrix(Egm))

plt.figure()
plt.hist(grains.I1)

plt.show()

