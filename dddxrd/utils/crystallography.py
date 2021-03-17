import numpy as np
import xfab.symmetry

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

#### CODE FOR DOING POLEFIGURES
def stereographic_projection(g,carthesian=True,polar=False):
    """ Stereographic projection of a coordinate vector"""
    gn = g/np.linalg.norm(g) #normalize the vector
    alpha = np.arccos(gn[2])
    beta = np.arccos(gn[0]/(gn[0]**2+gn[1]**2))
    phi = np.arctan2(gn[1],gn[0])
    r = np.tan(alpha/2.)
    px = r*np.cos(phi)
    py = r*np.sin(phi)
    if carthesian and polar:
        return px,py, r, phi
    elif carthesian:
        return px,py
    elif polar:
        return r,phi
    else:
        print('Must choose carhesian or/and polar coordinates in stereographic projection')
        raise AttributeError
    

def old_pf_projection(ub,h,rotation=7):
    """
    Stereographic projection for pole figure.
    """
    # h = h/np.linalg.norm(h)
    # hp = np.dot(ub,h)
    # ub = np.linalg.inv(ub)
    Rsym = xfab.symmetry.rotations(rotation)
    rs = []
    ths = []
    R=  np.pi/2 #radius of sphere
    for rot in Rsym:
        hp = np.dot(rot,h) #apply symmetry
        # hp = hp/np.linalg.norm(hp)
        hpr = np.dot(ub,hp) #translate from crystal to sample
        hpr = hpr/np.linalg.norm(hpr)

        #Stereographic projection from hkl to x,y

        # th = np.arccos(hpr[2])
        # ph = np.arctan2(hpr[1],hpr[0]) #this chooses the correct quadrant
        # xo = np.tan(th/2.)*np.cos(ph)
        # yo = np.tan(th/2.)*np.sin(ph)

        alpha = np.arccos(hpr[2])#-np.pi/2
        if alpha<-np.pi/2 or alpha > np.pi/2: continue #skip poles from south hemisphere?
            # alpha = np.pi-alpha
        # print alpha*180/np.pi #must be in (0,90)
        phi = np.arctan2(hpr[1],hpr[0]) #this chooses the correct quadrant
        # r = R*np.tan(alpha/2)
        # alpha = np.pi-alpha
        # phi = -phi
        r=R*np.tan(alpha/2)
        # xo = r*np.cos(phi)
        # yo = r*np.sin(phi)


        # xo = hpr[0]/(hpr[2]+np.sqrt(hpr[0]**2+hpr[1]**2+hpr[2]**2))
        # yo = hpr[1]/(hpr[2]+np.sqrt(hpr[0]**2+hpr[1]**2+hpr[2]**2))
        # xo = hpr[0]/(1-hpr[2])
        # yo = hpr[1]/(1-hpr[2])
        # xo = hpr[0]/(-hpr[2]+np.sqrt(hpr[0]**2+hpr[1]**2+hpr[2]**2))
        # yo = hpr[1]/(-hpr[2]+np.sqrt(hpr[0]**2+hpr[1]**2+hpr[2]**2))

        # r = np.sqrt(xo**2+yo**2)
        # theta = np.arctan2(yo,xo)
        rs.append(r)
        ths.append(phi)
    return ths,rs

def old_ipf_projection(ubi,h0, symmetry=7):
    """
    Stereographic projection for inverse pole figure. DOES NOT WORK AT THE MOMENT
    symmetry:
        4   tetragonal
        7   cubic
    """
    # h0 = h0/np.linalg.norm(h0)
    hp = np.dot((ubi),h0) #Translates h0 from sample to crystal
    # hp = hp/np.linalg.norm(hp)

    Rsym = xfab.symmetry.rotations(symmetry)

    R = 1#np.pi/2.
    rgb = np.zeros((3))

    for rot in Rsym:
        hpr = np.dot(rot,hp) #applies symmetry
        hpr = hpr/np.linalg.norm(hpr)

        alpha = np.arccos(hpr[2])
        #if alpha<-np.pi/2 or alpha > np.pi/2:
            #alpha = np.pi-alpha
        phi = np.arctan2(hpr[1],hpr[0]) #this chooses the correct quadrant
        r = R*np.tan(alpha/2)

        #xo = r*np.cos(phi)
        #yo = r*np.sin(phi)
        xo = np.tan(alpha/2)*np.cos(phi)
        yo = np.tan(alpha/2)*np.sin(phi)
        ro = np.sqrt((xo+1)**2+yo**2) #radius from [-1,0]
        ## alpha = np.pi-alpha
        ## phi = -phi
        ## r=R/np.tan(alpha/2)
        #xo = R*np.cos(phi)*np.tan(alpha/2)
        #yo = R*np.sin(phi)*np.tan(alpha/2)

        ## #
        ## th = np.arccos(hpr[2])
        ## ph = np.arctan2(hpr[1],hpr[0])
        ## xo = np.tan(th/2.)*np.cos(ph)
        ## yo = np.tan(th/2.)*np.sin(ph)

        ##
        ## # xo = hpr[0]/(hpr[2]+np.sqrt(hpr[0]**2+hpr[1]**2+hpr[2]**2))
        ## # yo = hpr[1]/(hpr[2]+np.sqrt(hpr[0]**2+hpr[1]**2+hpr[2]**2))
        ## # xo = hpr[0]/(1-hpr[2])
        ## # yo = hpr[1]/(1-hpr[2])
        ##
        ## r = np.sqrt(xo**2+yo**2)
        ## theta = np.arctan2(yo,xo)
        ## if xo < 0: #2 or 3 quadrant
        ##     theta = theta+np.pi
        ## elif yo<0: #4 quadrant
        ##     theta = theta + 2*np.pi
        #if phi >= 0 and phi <= 45*np.pi/180. and ro <= np.sqrt(2): #this should be the correct orientation
        if ro <= np.sqrt(2): #this should be the correct orientation
            #rgb = ipf_coloring(xo,yo, symmetry = symmetry)
            rgb[0]=((np.sqrt(2.0)-r)/(np.sqrt(2.0)-1))**.5
            rgb[1]=((1-4*phi/np.pi)*((r-1)/(np.sqrt(2.0)-1)))**.5
            rgb[2]=(4*phi/np.pi*((r-1)/(np.sqrt(2.0)-1)))**.5
            mx = max(rgb)
            rgb = rgb/mx
            return rgb,r,phi
    phi = np.nan
    r = np.nan
    return rgb,r,phi

        