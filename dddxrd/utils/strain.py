import numpy as np
import xfab.tools
from ImageD11 import finite_strain as fs
import matplotlib.pyplot as plt
import dddxrd.utils.crystallography as cry

def matrix_to_voigt(e):
    """ Returns the matrix in Voigt notation [xx,yy,zz,yz,xz,xy]"""
    return [e[0,0],e[1,1],e[2,2],e[1,2],e[0,2],e[0,1]]
def voight_to_matrix(e):
    return np.array([[e[0],e[5],e[4]],[e[5],e[1],e[3]],[e[4],e[3],e[2]]])

def calc_strain(ubi,d0,mr=1,ml=-1,finite_strain=True):
    """ Calculates the strain tensor for one grain in crystal and lab frame
    By default, the strain tensor in crystal fram is the Green-Lagrange: E = 1/2*(C-I)=1/2*(F^T*F-I)
    and in sample frame it is the Euler-Almansi strain: e = 1/2*(I-c)=1/2*(I-F^(-T)F^(-1))
    ubi: (UB)^-1 matrix of the grain
    d0: reference lattice parameters [a,b,c,alpha,beta,gamma]
     """
    if not finite_strain:
        U,eps = xfab.tools.ubi_to_u_and_eps(ubi,d0)
        return U,fs.e6_to_symm(eps)
    # make reference UBI
    ubi0 = np.diag([d0[0],d0[1],d0[2]]) #only for 90 degree cells
    assert (d0[3]==90) and (d0[3]==90) and (d0[3]==90), "Finite strains are only for cubic refence cells for now"
    ub0 = np.linalg.inv(ubi0)
    dg = fs.DeformationGradientTensor( ubi,ub0 )
    #F=dg.F
    #E = 0.5*(np.dot(F.T,F)-np.eye(3))
    #e = 0.5*(np.eye(3)-np.dot(np.linalg.inv(F).T,np.linalg.inv(F)))
    E = dg.finite_strain_ref(m=mr) # strain in crystal
    e = dg.finite_strain_lab(m=ml) # strain in sample
    return E, e

def tensor_invariants(e):
    """ Calculates invariants of a tensor. 
    If the tensor is strain/stress:
    I1 is the volumetric strain/stress
    J2 is the magnitude of shear strain/stress"""
    e = np.array(e)
    if e.ndim==1:
        e = voight_to_matrix(e)
    I1 = np.trace(e)
    I2 = 0.5*(I1**2-np.trace(np.linalg.matrix_power(e,2)))
    J2 = I1**2-2*I2
    return I1,J2

def _fill_C3(cs):
    """ helper function for filling the C matrix for isotropy and cubic anisotropy"""
    C = np.zeros((6,6))
    C[0,:] = [cs[0],cs[2],cs[2],0,0,0]
    C[1,:] = [cs[2],cs[0],cs[2],0,0,0]
    C[2,:] = [cs[2],cs[2],cs[0],0,0,0]
    C[3,:] = [0,0,0,cs[1],0,0]
    C[4,:] = [0,0,0,0,cs[1],0]
    C[5,:] = [0,0,0,0,0,cs[1]]
    return C
def _fill_C9(cs):
    """ helper function for filling the C matrix for transverse isotropy and orthotropy"""
    C = np.zeros((6,6))
    C[0,:] = [cs[0],cs[8],cs[7],0,0,0]
    C[1,:] = [cs[8],cs[1],cs[6],0,0,0]
    C[2,:] = [cs[7],cs[6],cs[2],0,0,0]
    C[3,:] = [0,0,0,cs[3],0,0]
    C[4,:] = [0,0,0,0,cs[4],0]
    C[5,:] = [0,0,0,0,0,cs[5]]
    return C

def elasticty_matrix(cs):
    """ Returns the 6x6 elasticity matrix from the coefficients cs. The lengts of coeff determines how the matrix looks:
    Isotropic: cs = [K,mu]
    Cubic anisotropy: cs = [C11,C44,C12]
    Transverse isotropy: cs = [C11,C33,C44,C66,C13]
    Orthotropic: cs = [C11,C22,C33,C44,C55,C66,C12,C13,C22,C23]
    """
    if len(cs) == 2:
        print('Setting up isotropic elasticity matrix')
        K = cs[0]
        mu = cs[1]
        d = K+4.*mu/3.
        od = K - 2.*mu/3.
        C = _fill_C3([d,mu,od])
    elif len(cs) == 3:
        print('Setting up cubic anisotropy')
        C = _fill_C3(cs)
    elif len(cs) == 5:
        print('Setting up transverse anisotropy')
        cc = [cs[0],cs[0],cs[1],cs[2],cs[2],cs[3],cs[4],cs[4],cs[0]-2*cs[3]]
        C = _fill_C9(cc)
    elif len(cs) == 9:
        C = _fill_C9(cs)
    else:
        raise ValueError("C matrix with {:d} independent parameters is not implemented.".format(len(cs)))
    return C
 
def calc_stress(C,E,PK2=True, cauchy=False, F=None):
    """Calculates the 2nd Piola-Kirchoff and/or Cauchy stress.
    Second Piola Kirchoff is:
    S = C:E and is defined in the crystal frame.
    Cauchy stress is:
    sigma = (1/det(F))*F*S*F^T) and is defined in the sample frame"""
    S = np.dot(C,E) # 2nd Piola-Kirchoff
    if cauchy:
        assert F is not None, "Need to supply the deformation gradient for computing Cauchy stress"
        S = voight_to_matrix(S)
        J = np.linalg.det(F)
        sigma = (1./J)*np.dot(F,np.dot(S,F.T))
        #sigma =  matrix_to_voigt(sigma)
    sigma = matrix_to_voigt(sigma)
    S = matrix_to_voigt(S)
    if PK2:
        if cauchy:
            return S, sigma
        else:
            return S
    else:
        return sigma


def strain_energy_density(S,E):
    """ Calculates the strain energy density as W=0.5*(S:E)."""
    S = np.array(S)
    E = np.array(E)
    if S.ndim == 1:
        S = voight_to_matrix(S)
    if E.ndim == 1:
        E = voight_to_matrix(E)
    return 0.5*np.tensordot(S,E,axes=2)
