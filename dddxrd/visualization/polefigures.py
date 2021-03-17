import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from numpy.random import uniform
import dddxrd.utils.crystallography as cry
import scipy.interpolate as sciint

def _ipf_triangle_cubic(npoints=21):
    """ Generate points for plotting the fundamental zone triangle. Does not include origin"""
    xa = np.zeros((npoints))
    ya = np.zeros((npoints))
    for i in range(npoints):
        ua = np.array([i/(npoints-1.), 1., 1.])
        UA = np.linalg.norm(ua)
        za = ua[2]+UA
        xa[i] = ua[1]/za
        ya[i] = ua[0]/za
    return xa, ya
def ipf_figure(coord,color,size):
    """ Plots and inverse pole figure where the size of the circles are proportional to the grain volumes.
    Only works for cubic materials"""
    #copied from plot_gff
    fig = plt.figure(10,frameon=False,figsize=plt.figaspect(.9))
    ax = plt.Axes(fig,[.2,.2,.7,.7])
    ax.set_axis_off()
    fig.add_axes(ax)
    #plot triangle
    xa,ya = _ipf_triangle_cubic(npoints=21)
    plt.plot(xa,ya,'black') # Curved edge
    plt.plot([0,xa[0]],[0,0.0001],'black',linewidth=2) #lower line
    plt.plot([xa[20],0],[ya[20],0],'black') # upper line
    plt.text(-0.01,-0.02,'[001]')
    plt.text(xa[0]-0.01,-0.02,'[101]')
    plt.text(xa[-1]-0.01,ya[-1]+0.005,'[111]')
    for c,rgb,s in zip(coord,color,size):
        plt.plot(c[1],c[0],'o',color=rgb,alpha=0.6,mec='k',ms=s)
    return fig,ax
    #plt.plot(coord[:,0],coord[:,1],'o',color=color)

def ipf_xfab(axis):
    #from xfab/plot_gff not sure how it works...
    rgb = np.zeros((3))
    for j in range(3):
        for k in range(j+1,3):
            if (axis[j]>axis[k]):
                rgb[0]=axis[j]
                axis[j]=axis[k]
                axis[k]=rgb[0]
    rr=np.sqrt(axis[0]**2/(axis[2]+1)**2+(axis[1]/(axis[2]+1)+1)**2)
    if axis[1]==0:
        beta=0
    else:
        beta=np.arctan(axis[0]/axis[1])
    rgb[0]=((np.sqrt(2.0)-rr)/(np.sqrt(2.0)-1))**.5
    rgb[1]=((1-4*beta/np.pi)*((rr-1)/(np.sqrt(2.0)-1)))**.5
    rgb[2]=(4*beta/np.pi*((rr-1)/(np.sqrt(2.0)-1)))**.5
    mx = max(rgb)
    rgb = rgb/mx
    px,py, r,phi = cry.stereographic_projection(axis,polar=True,carthesian=True)
    return rgb, px, py, r, phi

def pf_contour(x,y,num_bins=96,polar=False,**kwargs):
    """ Makes a contourplot of pole figure data..
    polar: x = theta, y = r 
    """
    if polar:
        fig,ax = _polar_contour(x,y,num_bins=num_bins,**kwargs)
    else:
        xs = np.linspace(x.min(),x.max(),num=num_bins,endpoint=True)
        ys = np.linspace(y.min(),y.max(),num=num_bins,endpoint=True)
        vals,xe,ye = np.histogram2d(x,y,bins=num_bins)
        xgr,ygr = np.meshgrid(xs,ys)
        fig,ax = plt.subplots()
        im = ax.contourf(xgr,ygr,vals,**kwargs)
        fig.colorbar(im)
    return fig,ax

def ipf_contour(x,y,nnodes=200,bins=20,symmetry='cubic',**kwargs):
    """ Makes a contourplot of inverse pole figure data in carthesian coordinates"""
    assert symmetry=='cubic', "ipf_contour only works for cubic materials at the moment"
    #generate fundamental zone
    xi,yi = _ipf_triangle_cubic(npoints=21)
    # add origin points
    xa = np.hstack((0,xi,0))
    ya = np.hstack((0,yi,0))
    # fill with nodes
    xn,yn = _tri_create_nodes([xa[1],ya[1]],[xa[-2],ya[-2]],nno=nnodes)
    xa = np.hstack((xa,xn))
    ya = np.hstack((ya,yn))
    # Create Delaunay triangulation
    triangles = tri.Triangulation(xa, ya)
    # bin the data
    H,xe,ye = np.histogram2d(x,y,bins=bins,density=False)
    xc = (xe[:-1] + xe[1:]) / 2
    yc = (ye[:-1] + ye[1:]) / 2
    # fit a spline
    f = sciint.interp2d(xc,yc,H)
    # Compute values at the nodes
    z = np.zeros(xa.shape)
    for i,(x,y) in enumerate(zip(xa,ya)):
        z[i]=f(x,y)
    #plot
    fig,ax = plt.subplots(frameon=False)
    ax.set_aspect('equal')
    im = ax.tricontourf(triangles, z)
    fig.colorbar(im)
    ax.plot(xi,yi,'black') # Curved edge
    ax.plot([0,xi[0]],[0,0.0001],'black',linewidth=2) #lower line
    ax.plot([xi[20],0],[yi[20],0],'black') # upper line
    ax.text(-0.01,-0.02,'[001]')
    ax.text(xi[0]-0.01,-0.02,'[101]')
    ax.text(xi[-1]-0.01,yi[-1]+0.005,'[111]')
    ax.set_axis_off()
    # plt.title('Contour plot of Delaunay triangulation')
    # plt.figure()
    # plt.plot(xa,ya,'o')
    return fig,ax

def _check_inside(v0,v1,v2,v):
    """ check if point v in inside the triangle v0-v1-v2
    from https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle"""
    area = 0.5*(-v1[1]*v2[0]+v0[1]*(-v1[0]+v2[0])+v0[0]*(v1[1]-v2[1])+v1[0]*v2[1])
    s = 1/(2*area)*(v0[1]*v2[0]-v0[0]*v2[1]+(v2[1]-v0[1])*v[0]+(v0[0]-v2[0])*v[1])
    t = 1/(2*area)*(v0[0]*v1[1]-v0[1]*v1[0]+(v0[1]-v1[1])*v[0]+(v1[0]-v0[0])*v[1])
    if (s>=0) and (t>=0) and (1-s-t >=0):
        return True
    else:
        return False

def _tri_create_nodes(v1,v2,nno=100):
    """ Create random node coordinates inside the triangle [0,0] - [v1] - [v2] """
    rng = np.random.default_rng(seed=12356)
    xa = np.zeros((nno))
    ya = np.zeros((nno))
    hit = 0
    while hit < nno:
        rnd = rng.uniform(0,1,2)
        x = rnd[0]*v1[0]+rnd[1]*v2[0]
        y = rnd[0]*v1[1]+rnd[1]*v2[1]
        #check if inside triangle
        if _check_inside([0,0],v1,v2,[x,y]):
            xa[hit] = x
            ya[hit] = y
            hit += 1
    return xa,ya

def _polar_contour(thetas, rs, num_bins=96, theta_lim=[0,2*np.pi], r_lim=[0,1], radians=True, return_data = False):
    """
    Makes a filled contourplot of pole figure data
    """
    if not radians:
        thetas *= np.pi/180. #makes theta to radians

    r =  np.linspace(*r_lim,num=num_bins,endpoint=True)
    th =  np.linspace(*theta_lim,num=num_bins,endpoint=True)
    vals,xe,ye = np.histogram2d(thetas,rs,bins=num_bins) #bin ipf data in theta - r polar coordinates
    rgr,thgr = np.meshgrid(r,th) #make a grid
    if return_data:
        return thgr, rgr, vals
    else:
        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        im = ax.contourf(thgr,rgr,vals)
        fig.colorbar(im)
        return fig,ax 
    