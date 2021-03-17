import dddxrd.utils.crystallography as cry
import dddxrd.visualization.polefigures as pf
import dddxrd.visualization.Grainmap as gm
import matplotlib.pyplot as plt
import numpy as np
import sys

# hh = np.arange(-10,11)
# kk = np.arange(-10,11)
# ll = np.arange(0,11)

# cc = []
# pc = []
# for h in hh:
#     for k in kk:
#         for  l in ll:
#             if (h==0) and (k==0) and (l==0):
#                 continue
#             else:
#                 print(h,k,l)
#                 px,py,r,beta = cry.stereographic_projection([h,k,l],carthesian=True,polar=True)
#                 cc.append([px,py])
#                 pc.append([r,beta])
# cc = np.array(cc)
# pc = np.array(pc)
# print(cc)
# print(pc[0,:])
# plt.figure()
# plt.plot(cc[:,0],cc[:,1],'.')

# plt.figure()
# plt.polar(pc[:,1],pc[:,0],'.')

mapfile = "dddxrd/tests/test.map"
grains = gm.Grainmap(mapfile)
x=np.array(grains.xy)[:,0]
y=np.array(grains.xy)[:,1]
fig,ax=pf.ipf_contour(x,y,nnodes=100,bins=20,symmetry='cubic')
plt.show()



