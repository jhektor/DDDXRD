import numpy as np
import xfab.symmetry
from scipy import stats

def crystal_to_sample(X,U):
    return np.dot(U,np.dot(X,U.T))

def rodrigues_rotation(n,theta,radians=True):
    """ Returns a 3D rotation matrix according to the Rodrigues formula. 
    n: rotation axis
    theta: rotation angles
    """
    if not radians:
        theta *= np.pi/180
    n = n/np.linalg.norm(n)
    R = np.zeros((3,3))
    ct = np.cos(theta)
    st = np.sin(theta)
    R[0,:] = [ct + n[0]**2*(1-ct), n[0]*n[1]*(1-ct)-n[2]*st, n[0]*n[2]*(1-ct)+n[1]*st]
    R[1,:] = [n[0]*n[1]*(1-ct)+n[2]*st, ct + n[1]**2*(1-ct), n[1]*n[2]*(1-ct)-n[0]*st]
    R[2,:] = [n[0]*n[2]*(1-ct)-n[1]*st, n[1]*n[2]*(1-ct)+n[0]*st, ct + n[2]**2*(1-ct)]
    return R

def energy_to_wavelength(energy):
    ''' Convert energy (keV) to wavelength (Å)'''
    lam = 1.2398/energy # in nm
    lam *= 10
    return lam

def braggs_law(wavelength, theta=None, degree=False, d=None):
    assert theta or d, 'Must give either theta or d'
    if degree:
        theta = np.deg2rad(theta)
    if not d:
        return wavelength/(2*np.sin(theta))
    else:
        tth=2*(np.arcsin(wavelength/(2*d))) #this is two theta
        if degree:
            tth = np.rad2deg(tth)
        return tth

def check_reflection(hkl,symmetry='fcc'):
    '''Check if a refection is allowd of forbidden based on the symmetry. Only fcc or bcc'''
    assert symmetry in ['primitive','fcc','bcc'], 'check_reflection only implemented for primitive, fcc or bcc'
    h,k,l=hkl[0], hkl[1],hkl[2]
    if h+k+l==0:
        allowed=False
    else:
        if symmetry=='bcc':
            allowed=(h+k+l)%2==0
        if symmetry=='fcc':
            allowed = ((h%2==1 and k%2==1 and l%2==1) or (h%2==0 and k%2==0 and l%2==0)) # all odd or all even
    return allowed


def angle_to_latpar(angle,energy,symmetry='cubic',hkl=[1,1,1],two_theta=True, degrees=True):
    '''Convert diffraction angle to lattice parameter. Does not check for forbidden reflections. Only cubic for now.'''
    if two_theta:
        angle *= 0.5 #theta
    if degrees:
        angle = np.deg2rad(angle)
    lam = energy_to_wavelength(energy)
    d = braggs_law(lam,angle)
    if symmetry in ['cubic','fcc','bcc']:
        a = d*np.sqrt(hkl[0]**2+hkl[1]**2+hkl[2]**2)
    else:
        print('angle_to_latpar is only implemented for cubic crystals')
        raise ValueError
    return a

def d_spacing(hkl,latticeparameter, symmetry='cubic'):
    assert symmetry in ['cubic','fcc','bcc'], 'd_spacing only works for cubic crystals'
    d2=latticeparameter**2/(hkl[0]**2+hkl[1]**2+hkl[2]**2)
    return np.sqrt(d2)
    

def average_cell(grains,make_cubic=False):
    """ Calculates the average unit cell parameters for a list of grains """
    ls = []
    for g in grains:
        ls.append(g.unitcell)
    ls = np.array(ls)
    for i,r in enumerate(ls.T):
        # z-score to remove outliers
        z = np.abs(stats.zscore(r))
        r = np.where(z>3,np.nan,r)
        ls.T[i] = r
    avg_cell = np.nanmean(ls,axis=0)
    if make_cubic:
        a = (avg_cell[0]+avg_cell[1]+avg_cell[2])/3
        avg_cell = [a,a,a,90,90,90]
    print('Average parameters: ',avg_cell)
    return avg_cell

def stereographic_projection(g,carthesian=True,polar=False):
    """ Stereographic projection of a coordinate vector"""
    gn = g/np.linalg.norm(g) #normalize the vector
    alpha = np.arccos(gn[2])
    beta = np.arccos(gn[0]/(gn[0]**2+gn[1]**2))
    phi = np.arctan2(gn[1],gn[0])
    r = np.tan(alpha/2.)
    px = gn[0]/(1+gn[2])#r*np.cos(phi)
    py = gn[1]/(1+gn[2])#r*np.sin(phi)    
    # px1 = r*np.cos(phi)
    # py1 = r*np.sin(phi)
    #print(px1,px,py1,py)
    if carthesian and polar:
        return px,py, r, phi
    elif carthesian:
        return px,py
    elif polar:
        return r,phi
    else:
        print('Must choose carhesian or/and polar coordinates in stereographic projection')
        raise AttributeError
    
