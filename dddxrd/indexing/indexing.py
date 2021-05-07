import os
import subprocess
import dddxrd.utils.parser as parser
import dddxrd.utils.plotting as plu
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ImageD11.grain import read_grain_file, write_grain_file
import shutil

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

def match_grains(mapfiles,outfiles,dist_tol=100,ang_tol=0.5):
    """ Match and relabel grains across map-files.
    mapfiles: list of map files to match. Will edit the "name" key of the grain object so that the name is the same for common grains.
    outfiles: list of filenames to save the matched maps 
    dist_tol: tolerance for distance between CoM position
    ang_tol: tolerance for misorientation angles """
    nmaps = len(mapfiles)
    assert nmaps>1, 'Must provide more than one mapfile to match_grains'
    #rename grains in reference map
    for i in range(nmaps-1):
        nmatch = 0
        #match grains from two consecutive map-files
        match_map = mapfiles[i+1]
        print('Match map: {}'.format(match_map))
        if i ==0:
            ref_map = mapfiles[i]
            print('Reference map: {}'.format(ref_map))
            gm_ref = read_grain_file(ref_map)
            print('Relabelling first reference map')
            for j,g in enumerate(gm_ref):
                g.name = str(j)
            print('Saving reference .map to {}'.format(outfiles[i]))
            write_grain_file(outfiles[i],gm_ref)
            ngrains = len(gm_ref)
        gm_match = read_grain_file(match_map)
        #loop over match map
        print('Starting matching. This might take a while...')
        for k,g in enumerate(gm_match): #Parallelize
            if k%10 == 0:
                print('Matching grain {:d} of {:d}'.format(k,len(gm_match)))
            r = g.Rod
            x,y,z = g.translation
            # print('Matching grain {}'.format(g.name))
            dmin = dist_tol
            amin = ang_tol
            idx = None
            #loop over reference map
            for j,gr in enumerate(gm_ref):
                r0 = gr.Rod
                x0,y0,z0 = gr.translation
                dist=np.sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2) #distance between CoM in ref and match map
                ang = np.arccos(np.dot(r0,r)/(np.linalg.norm(r)*np.linalg.norm(r0)))*180/np.pi #angle between Rodriguez vectors
                #find the best match
                if (dist<dmin) and (ang<amin):
                    dmin=dist
                    amin=ang
                    idx = j
            if idx: #match these grains
                nmatch += 1
                match = gm_ref[idx]
                g.name = match.name
            else: # new grain
                ngrains += 1
                g.name = str(ngrains)


        print('Matched {:d} out of {:d} grains'.format(nmatch,len(gm_match)))
        print('Saving matched .map to {}'.format(outfiles[i+1]))
        write_grain_file(outfiles[i+1], gm_match) 
        gm_ref = gm_match
    print('The total dataset contains {:d} grains.'.format(ngrains))
    return

def main():
    par_file = 'dddxrd/tests/indexing.yaml'
    index_and_map(par_file=par_file)

if __name__ == "__main__":
    main()