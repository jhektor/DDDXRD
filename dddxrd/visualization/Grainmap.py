from __future__ import print_function, division
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
#matplotlib.use('TkAgg')

import numpy as np
import glob
from ImageD11.grain import read_grain_file, write_grain_file
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import sys
import dddxrd.visualization.Cube as Cube
import dddxrd.utils.plotting as plu
import dddxrd.visualization.polefigures as pf
import dddxrd.utils.strain as strain
import dddxrd.utils.crystallography as cry
import copy

class Grainmap:
    def __init__(self,mapfile,d0=None):
        self.grains = read_grain_file(mapfile)
        #center of mass
        com = []
        npks = []
        for gr in self.grains:
            #com.append(gr.translation)
            npks.append(int(gr.npks))
        #self.com = np.array(com)
        self.ngrains = len(self.grains)
        self._strains(d0)
        self._make_cubes()

    def _strains(self,d0=None):
        self.I1 = []
        self.J2 = []
        self.Green_strain = []
        self.Almansi_strain = []
        if d0 is None:
            print('No reference lattice parameters supplied. Will use the average of all grains')
            d0 = cry.average_cell(self.grains,make_cubic=True)
        self.d0=d0
        for gr in self.grains:
            E,e = strain.calc_strain(gr.ubi,d0)
            i1,j2 = strain.tensor_invariants(e)
            self.I1.append(i1)
            self.J2.append(j2)
            self.Green_strain.append(strain.matrix_to_voigt(E))
            self.Almansi_strain.append(strain.matrix_to_voigt(e))
        return

    def com_axis(self,bins=100):
        # bin z coordinate
        _,bin = np.histogram(self.com[:,2],bins=bins)
        place = np.digitize(self.com[:,2],bin)
        com_axis = np.zeros((bins,3))
        for i in range(bins):
            com = self.com[np.where(place==i+1)[0],:]
            size = self.size[np.where(place==i+1)[0]]
            for j in range(3):
                com_axis[i,j] = np.sum(com[:,j]*size)/np.sum(size)
                #com_axis[i,j] = np.sum(com[:,j])/com.shape[0]
        return com_axis

    def update_com(self):
        for g,com in zip(self.grains,self.com):
            g.translation = com
        self._make_cubes()

    def _make_cubes(self):
        self.cubes = []
        self.size = []
        self.com = []
        self.xy = [] #ipf coordinates in carhesian
        self.rphi = []
        self.rgb = []
        for g in self.grains:
            r = np.sqrt(g.translation[0]**2+g.translation[1]**2+g.translation[2]**2)
            if 1:#r<4000:#1:#np.abs(g.translation[2])<30:
                c = Cube.Cube(g)
                c.name = g.name.strip()
                # if c.size>35: continue
                self.cubes.append(c)
                self.size.append(c.size)
                self.com.append(c.com)
                self.xy.append([c.px,c.py])
                self.rphi.append([c.r,c.phi])
                self.rgb.append(c.ipf_color)

        self.size = np.array(self.size)
        self.com = np.array(self.com)

    def plot_3d_map(self,coloring="ipf",alpha=0.6,linewidths=0.2,edgecolors='k',cmap=cm.viridis,filter=None,**kwargs):
        fig,ax = self._prepare_figure(projection="3d")
        bb = []
        alphal = np.full((self.ngrains),alpha)
        if filter:
            for i,c in enumerate(self.cubes):
                if c.name not in filter:
                    alphal[i]=0 #make grains transparent if filtered out
        if coloring != "ipf":
            if not 'vmin' in kwargs:
                vmin = np.min(coloring)
            else:
                vmin = kwargs['vmin']
            if not 'vmax' in kwargs:
                vmax = np.max(coloring)
            else:
                vmax = kwargs['vmax']
            norm = matplotlib.colors.Normalize(clip=True, vmin=vmin,vmax=vmax)
            mapper = cm.ScalarMappable(norm=norm, cmap=cmap)
            plt.colorbar(mapper)
        for i,c in enumerate(self.cubes):
            bb.append(c.get_bounding_box())
            faces = Poly3DCollection(c.edges, linewidths=linewidths, edgecolors=edgecolors)
            if coloring == "ipf":
                #print(c.ipf_color)
                if len(c.ipf_color) == 3:
                    faces.set_facecolor(tuple(c.ipf_color)+(alphal[i],))
                else:
                    faces.set_facecolor(tuple(c.ipf_color))
            else: #for coloring based on scalar data
                faces.set_facecolor(mapper.to_rgba(coloring[i],alpha=alphal[i]))
            ax.add_collection3d(faces)
        bb = np.array(bb)
        #print(bb.shape)
        ll = np.min(bb[:,:,0],axis=0)
        ul = np.max(bb[:,:,1],axis=0)
        ax.set_xlim(ll[0],ul[0])
        ax.set_ylim(ll[1],ul[1])
        ax.set_zlim(ll[2],ul[2])
        plu.set_axes_equal(ax)
        #ax.set_aspect('equal')
        ax.grid(False)
        return fig,ax

    def _prepare_figure(self,**kwargs):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection=kwargs.get("projection"))
        if kwargs.get("xlabel"):
            ax.set_xlabel(kwargs.get("xlabel"))
        if kwargs.get("ylabel"):
            ax.set_ylabel(kwargs.get("ylabel"))
        if kwargs.get("zlabel"):
            ax.set_zlabel(kwargs.get("zlabel"))
        try:
            ax.view_init(elev=kwargs.get("elev"), azim=kwargs.get("azim"))
        except AttributeError:
            pass
        return fig,ax

    def scatterplot(self,coords,cbar=True,**kwargs):
        #Check if 3D
        if coords.shape[1]==3:
            kwargs["projection"] = "3d"
        fig,ax = self._prepare_figure(**kwargs)
        #remove stuff from kwargs because scatter is not well implemented
        for key in ["projection","xlabel","ylabel","zlabel","elev","azim"]:
            if key in kwargs:
                del kwargs[key]
        if coords.shape[1]==3:
            ax.scatter3D(coords[:,0],coords[:,1],coords[:,2],**kwargs)
            plu.set_axes_equal(ax)

        else:
            im=ax.scatter(coords[:,0],coords[:,1],**kwargs)
            ax.set_aspect(1)
            if cbar:
                fig.colorbar(im)
        return fig,ax

def main(maps):
    gms = []
    stems = []
    for g in maps:
        print('Adding {} to figure list'.format(g))
        gms.append(Grainmap(g))
        stems.append(g.split('.map')[0])

    #gm = gms[0]
    #stem = stems[0]
    # print(gm.ngrains)
    for gm,stem in zip(gms,stems):
        print('Making plots of {}'.format(gm))
        fig,ax = gm.plot_3d_map()
        #fig.savefig('{}_3d.svg'.format(stem),dpi=300,format='svg',bbox_inches='tight')

        view = ('top','side','front')
        cind = ([0,1],[0,2],[1,2])
        lab = (['x','y'],['x','z'],['y','z'])

        for v,c,l in zip(view,cind,lab):
            #print(gm.com[:,c[1]].shape,gm.size.shape)
            #print(np.array([gm.com[:,c[0]],gm.com[:,c[1]]]).T.shape)
            fig,ax = gm.scatterplot(np.array([gm.com[:,c[0]],gm.com[:,c[1]]]).T,c=gm.rgb,s=gm.size/2,alpha=0.6,edgecolor='k',marker='o')
            ax.set_xlabel(l[0])
            ax.set_ylabel(l[1])
            #fig.savefig('{}_{}.svg'.format(stem,v),dpi=300,format='svg',bbox_inches='tight')

        gm.xy = np.array(gm.xy)
        fig,ax = pf.ipf_figure(gm.xy,gm.rgb,gm.size/10)

        fig,ax=pf.ipf_contour(gm.xy[:,0],gm.xy[:,1],nnodes=100,bins=20,symmetry='cubic')
        #fig.savefig('{}_ipf_all.svg'.format(stem),dpi=300,format='svg',bbox_inches='tight')

        plt.figure()
        vals,be,_ = plt.hist(gm.size/gm.size.max(),bins=100)#gm.ngrains//10)
        bc = (be[:-1] + be[1:]) / 2
        plt.xlabel('Relative grain size')
        plt.ylabel('Number of grains')
        #plt.savefig('{}_size_hist.svg'.format(stem),dpi=300,format='svg',bbox_inches='tight')
        header = 'relative grain size, nbr grains'
        #np.savetxt('{}_size_hist.csv'.format(stem),np.column_stack((bc,vals)),fmt=['%.3f','%d'],header=header,delimiter=',')
        #np.savetxt('{}_sizes.csv'.format(stem),gm.size)
        print('done')
        plt.show()

if __name__=='__main__':
    main(sys.argv[1:])
        