# methods to draw a sphere with Matplotlib (from https://stackoverflow.com/questions/44881885/python-draw-parallelepiped)
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import dddxrd.visualization.polefigures as pf
import dddxrd.utils.crystallography as cry
from ImageD11.grain import read_grain_file
import dddxrd.utils.plotting as plu

import sys
import xfab.symmetry 

xa,ya = pf._ipf_triangle_cubic()
v0 =[0,0]
v1 = [np.max(xa),0]
v2 = [xa[np.argmax(ya)],np.max(ya)]


class Sphere:
    def __init__(self,grain,axs=[0,0,1]):
        self.ubi = grain.ubi
        try:
            self.u = grain.U
        except AttributeError:
            self.u = grain.u
        #center of mass
        self.com = grain.translation
        try:
            intensity = float(grain.intensity_info.split(',')[4].split('=')[-1])
        except IndexError:
            intensity = 0
        self.size =  intensity**(1./3.)*5
        self.points,self.edges = self._sphere_points()
        self.ipf_color,self.px,self.py,self.r,self.phi = self._ipf_color(axs=axs)


    def _sphere_points(self):
        x0,y0,z0 = self.com
        V = self.size
        r = 5*(0.75*V/np.pi)**(1./3)

        phi = np.linspace(0,360, 20)/180.0*np.pi
        theta = np.linspace(-90,90,10)/180.0*np.pi

        points = [self.com]
        edges = []
        for i  in range(len(phi)-1):
            for j in range(len(theta)-1):
                cp0= r*np.cos(phi[i])
                cp1= r*np.cos(phi[i+1])
                sp0= r*np.sin(phi[i])
                sp1= r*np.sin(phi[i+1])

                ct0= np.cos(theta[j])
                ct1= np.cos(theta[j+1])

                st0= r*np.sin(theta[j])
                st1= r*np.sin(theta[j+1])

                verts=[]
                verts.append((cp0*ct0+x0, sp0*ct0+y0, st0+z0))
                verts.append((cp1*ct0+x0, sp1*ct0+y0, st0+z0))
                verts.append((cp1*ct1+x0, sp1*ct1+y0, st1+z0))
                verts.append((cp0*ct1+x0, sp0*ct1+y0, st1+z0))
                edges.append(verts)
        
        return points,edges

    def _ipf_color(self,axs=[0,0,1]):
        axis = np.abs(np.dot(self.u.T,axs)) 
        rgb,px,py,r,phi = pf.ipf_xfab(axis)
        return rgb,px,py,r,phi

    def get_bounding_box(self):
        x_min,y_min,z_min = np.min(self.points,axis=0)
        x_max,y_max,z_max = np.max(self.points,axis=0)

        max_range = np.array(
            [x_max-x_min, y_max-y_min, z_max-z_min]).max() / 2.0

        mid_x = (x_max+x_min) * 0.5
        mid_y = (y_max+y_min) * 0.5
        mid_z = (z_max+z_min) * 0.5

        return [
            [mid_x - max_range, mid_x + max_range],
            [mid_y - max_range, mid_y + max_range],
            [mid_z - max_range, mid_z + max_range]
        ]


    def plot_sphere(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        faces = Poly3DCollection(self.edges, linewidths=1, edgecolors='k')
        faces.set_facecolor(tuple(self.ipf_color)+(0.4,))
        # faces.set_alpha(0.2)

        ax.add_collection3d(faces)

        bounding_box = self.get_bounding_box(self.points)

        ax.set_xlim(bounding_box[0])
        ax.set_ylim(bounding_box[1])
        ax.set_zlim(bounding_box[2])

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_aspect('equal')
        plt.show()

if __name__ == '__main__':
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    g = read_grain_file(sys.argv[1])
    bb = []
    for i,gr in enumerate(g):
        c = Sphere(gr)
        bb.append(c.get_bounding_box())
        faces = Poly3DCollection(c.edges, linewidths=0.2, edgecolors='k')
        faces.set_facecolor(tuple(c.ipf_color)+(0.6,))

        ax.add_collection3d(faces)
    bb = np.array(bb)
    ll = np.min(bb[:,:,0],axis=0)
    ul = np.max(bb[:,:,1],axis=0)
    print(ll,ul)
    ax.set_xlim(ll[0],ul[0])
    ax.set_ylim(ll[1],ul[1])
    ax.set_zlim(ll[2],ul[2])
    #ax.set_aspect('equal')
    plu.set_axes_equal(ax)

    ax.grid(False)
    plt.show()
