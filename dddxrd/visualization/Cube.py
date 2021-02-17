# methods to draw a cube with Matplotlib (from https://stackoverflow.com/questions/44881885/python-draw-parallelepiped)
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from ImageD11.grain import read_grain_file
import sys

class Cube:
    def __init__(self,grain):
        self.ubi = grain.ubi
        self.u = grain.U
        #center of mass
        self.com = grain.translation
        try:
            intensity = float(grain.intensity_info.split(',')[4].split('=')[-1])
        except IndexError:
            intensity = 0
        self.size =  intensity**(1./3.)/5
        self.points,self.edges = self._cube_points()
        self.ipf_color,self.r,self.phi = self._ipf_color()


    def _cube_points(self):
        vecs = self.size*self.u #base vectors of the cube
        start = self.com - 0.5*(vecs[0]+vecs[1]+vecs[2]) #this will put the com in the midle of the cube
        # find all corner points
        points = [start,
                  start + vecs[0],
                  start + vecs[0] + vecs[1],
                  start + vecs[1],
                  start + vecs[2],
                  start + vecs[0] + vecs[2],
                  start + vecs[0] + vecs[1] + vecs[2],
                  start + vecs[1] + vecs[2]]

        edges = [[points[0], points[1], points[2], points[3]],
                 [points[0], points[3], points[7], points[4]],
                 [points[0], points[1], points[5], points[4]],
                 [points[6], points[7], points[3], points[2]],
                 [points[6], points[5], points[4], points[7]],
                 [points[6], points[5], points[1], points[2]]
                  ]
        return points,edges

    def _ipf_color(self):
        #from xfab/plot_gff not sure how it works...
        axis = abs(self.u[2,:])
        rgb = np.zeros((3))
        for j in range(3):
            for k in range(j+1,3):
                if (axis[j]>axis[k]):
                    rgb[0]=axis[j]
                    axis[j]=axis[k]
                    axis[k]=rgb[0]
        r=np.sqrt(axis[0]*axis[0]/((axis[2]+1))/((axis[2]+1))+(axis[1]/(axis[2]+1)+1)*(axis[1]/(axis[2]+1)+1))
        if axis[1]==0:
            phi=0
        else:
            phi=np.arctan(axis[0]/axis[1])
        rgb[0]=((np.sqrt(2.0)-r)/(np.sqrt(2.0)-1))**.5
        rgb[1]=((1-4*phi/np.pi)*((r-1)/(np.sqrt(2.0)-1)))**.5
        rgb[2]=(4*phi/np.pi*((r-1)/(np.sqrt(2.0)-1)))**.5
        mx = max(rgb)
        rgb = rgb/mx
        phi = axis[1]/(1+axis[2])
        r = axis[0]/(1+axis[2])
        return rgb,r,phi

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


    def plot_cube(self):
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
    for gr in g:
        c = Cube(gr)
        bb.append(c.get_bounding_box())
        faces = Poly3DCollection(c.edges, linewidths=0.2, edgecolors='k')
        faces.set_facecolor(tuple(c.ipf_color)+(0.6,))

        ax.add_collection3d(faces)
    bb = np.array(bb)
    ll = np.min(bb[:,:,0],axis=0)
    ul = np.max(bb[:,:,1],axis=0)
    ax.set_xlim(ll[0],ul[0])
    ax.set_ylim(ll[1],ul[1])
    ax.set_zlim(ll[2],ul[2])
    ax.set_aspect('equal')
    ax.grid(False)
    plt.show()
