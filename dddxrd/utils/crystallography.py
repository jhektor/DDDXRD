import numpy as np

def energy_to_wavelength(energy):
    ''' Convert energy (keV) to wavelength (Å)'''
    lam = 1.2398/energy # in nm
    lam *= 10
    return lam

def braggs_law(wavelength, theta, degree=False):
    if degree:
        theta = np.deg2rad(theta)
    return wavelength/(2*np.sin(theta))

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

        