import os
import subprocess
import dddxrd.utils.parser as parser
import dddxrd.utils.plotting as plu
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

lunarc = True
if lunarc:
    from dddxrd.indexing.grid_index_parallel import grid_index_parallel
else:
    from ImageD11.grid_index_parallel import grid_index_parallel

def make_grid(xrange,yrange,zrange,xstep,ystep,zstep,shape='cube',plot=False):
    """ Make a grid. Should perhaps be moved to utils"""
    cube = [(x,y,z) 
        for x in range(xrange[0],xrange[1],xstep)
        for y in range(yrange[0], yrange[1], ystep)
        for z in range(zrange[0], zrange[1], zstep)
        ]
    if shape == 'cube':
        grid = cube
    elif shape == 'cylinder':
        xr = np.abs(xrange[1]-xrange[0])
        yr = np.abs(yrange[1]-yrange[0])
        r = 0.5*np.mean([xr,yr])
        print(r)
        grid = [(x,y,z) for (x,y,z) in cube if x**2+y**2<r**2]
    else:
        print('grid shape not implemented')
        raise ValueError

    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        gg = np.array(grid)
        ax.scatter3D(gg[:,0],gg[:,1],gg[:,2],marker = 'x')
        plu.set_axes_equal(ax)
        plt.show()
    return grid

def index_and_map(par_file=None, pars=None):
    """ Performs indexation and mapping using grid_index_parallel from ImageD11"""
    if (pars is None) and (par_file is not None):
        pars = parser.parse_parameters(par_file)
    elif (pars is None) and (par_file is None):
        print('Must supply either par_file or pars to run_peaksearcher')
        raise ValueError
    grid = make_grid(pars['xrange'],pars['yrange'],pars['zrange'],pars['xstep'],pars['ystep'],pars['zstep'],plot=pars['plotgrid'],shape=pars['gridshape'])
    grid_index_parallel(pars['flt_file'],pars['par_file'],pars['stem'],pars,grid)

def main():
    par_file = 'dddxrd/tests/indexing.yaml'
    index_and_map(par_file=par_file)

if __name__ == "__main__":
    main()