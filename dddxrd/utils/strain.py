import numpy as np
import xfab.tools
from ImageD11 import finite_strain as fs

def matrix_to_voigt(e):
    """ Returns the matrix in Voigt notation [xx,yy,zz,yz,xz,xy]"""
    return [e[0,0],e[1,1],e[2,2],e[1,2],e[0,2],e[0,1]]
def voight_to_matrix(e):
    return np.array([[e[0],e[5],e[4]],[e[5],e[1],e[3]],[e[4],e[3],e[2]]])
def average_cell(grains,make_cubic=False):
    """ Calculates the average unit cell parameters for a list of grains """
    avg_cell = [0,0,0,0,0,0]
    for g in grains:
        avg_cell += g.unitcell
    avg_cell /= len(grains)
    if make_cubic:
        a = (avg_cell[0]+avg_cell[1]+avg_cell[2])/3
        avg_cell = [a,a,a,90,90,90]
    return avg_cell

def calc_strain(ubi,d0,mr=1,ml=-1,finite_strain=True):
    """ Calculates the strain tensor for one grain in crystal and lab frame
    The strain tensor in crystal fram is the Green-Lagrange: E = 1/2*(C-I)=1/2*(F^T*F-I)
    In sample frame it is the Euler-Almansi strain: e = 1/2*(I-c)=1/2*(I-F^(-T)F^(-1))
    ubi: (UB)^-1 matrix of the grain
    d0: reference lattice parameters [a,b,c,alpha,beta,gamma]
     """
    if not finite_strain:
        U,eps = xfab.tools.ubi_to_u_and_eps(ubi,d0)
        return finite_strain.e6_to_symm(eps)
    # make reference UBI
    ubi0 = np.diag([d0[0],d0[1],d0[2]]) #only for 90 degree cells
    assert (d0[3]==90) and (d0[3]==90) and (d0[3]==90), "Finite strains are only for cubic refence cells for now"
    ub0 = np.linalg.inv(ubi0)
    F = fs.DeformationGradientTensor( ubi,ub0 )
    E = F.finite_strain_ref(m=mr) # strain in crystal
    e = F.finite_strain_lab(m=ml) # strain in lab
    return E, e

def tensor_invariants(e):
    """ Calculates invariants of a tensor. 
    If the tensor is strain/stress:
    I1 is the volumetric strain/stress
    J2 is the magnitude of shear strain/stress"""
    I1 = np.trace(e)
    I2 = 0.5*(I1**2-np.trace(np.linalg.matrix_power(e,2)))
    J2 = I1**2-2*I2
    return I1,J2