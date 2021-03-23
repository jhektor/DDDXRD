import dddxrd.utils.strain as strain
import dddxrd.utils.crystallography as cry
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

d0 = [1,1,1,90,90,90] # reference lattice parameters
ubi0 = np.array([[d0[0],0,0],[0,d0[1],0],[0,0,d0[2]]]) # ubi aligned with lab frame
R = cry.rodrigues_rotation([1,1,1],30,radians=False) # Rotation matrix

nu = 0.3 # Poisson ratio
eps = np.linspace(0,0.5) # applied strain
fig,axes = plt.subplots(3,1,sharex=True)
ax = axes.ravel()
ax[0].set_title('[0,0]')
ax[1].set_title('[1,1]')
ax[2].set_title('[2,2]')
legend = ('E','e','eps crystal','eps lab')

ubi = np.dot(np.eye(3),ubi0)
for ep in eps:
    # elongate the c-axis
    ubi[2,2] = (1+ep)*ubi0[2,2]
    # compress the a and b
    ubi[0,0] = (1-nu*ep)*ubi0[0,0]
    ubi[1,1] = (1-nu*ep)*ubi0[1,1]
    #rotate grain 45 degrees around 111
    ubi = np.dot(R,ubi.T).T
    # Calculate strain
    E,e = strain.calc_strain(ubi,d0,mr=1,ml=-1,finite_strain=True)
    U,eps_xfab = strain.calc_strain(ubi,d0,finite_strain=False)
    eps_xfab_lab = np.dot(U,np.dot(eps_xfab,U.T))
    # plot
    ax[0].plot(ep,E[0,0],'ro')
    ax[0].plot(ep,e[0,0],'b+')
    ax[0].plot(ep,eps_xfab[0,0],'gs')
    ax[0].plot(ep,eps_xfab_lab[0,0],'kd')
    ax[1].plot(ep,E[1,1],'ro')
    ax[1].plot(ep,e[1,1],'b+')
    ax[1].plot(ep,eps_xfab[1,1],'gs')
    ax[1].plot(ep,eps_xfab_lab[1,1],'kd')
    ax[2].plot(ep,E[2,2],'ro')
    ax[2].plot(ep,e[2,2],'b+')
    ax[2].plot(ep,eps_xfab[2,2],'gs')
    ax[2].plot(ep,eps_xfab_lab[2,2],'kd')
    # rotate back
    ubi = np.dot(np.eye(3),ubi0)
print(ubi)
ax[2].legend(legend)
ax[2].set_xlabel('Applied strain')
ax[1].set_ylabel('Strain component')

plt.show()