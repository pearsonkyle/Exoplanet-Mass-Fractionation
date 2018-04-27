import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq as findzero
from scipy.interpolate import interp1d
from scipy.integrate import trapz
from astropy.io import fits

kb = 1.38e-16  # g cm2 s-2 K-1 [erg/K]
G = 6.674e-8 # cm3 g-1 s-2
h = 6.626e-27 # [erg/s]
c = 2.9979e10 # cm s-1
mp = 1.67e-24 # g
# barye = g cm-1 s-2
rsun = 6.957e10 # cm 
rjup = 6.9911e9 # cm 
mjup = 1.898e30 # g
rearth = 637.1e6 # cm
mearth = 5.97e28 # g
AU = 1.496e+13 # cm 
pc = 3.086e+18 # cm 

# load in atomic data for kinetic radius  
datatype = [
    ('atomic number','i8'),
    ('atom','U3'),
    ('mass','f8'),
    ('electronegativity','f8'),
    ('radius','f8'),
    ('radius error','f8')
]
atom_data = np.genfromtxt('atomic van der waal radii.csv',delimiter=',',dtype=datatype,skip_header=1)

class AttributeDict(dict):
    def __getattr__(self, name):
        return self[name]
    def __setattr__(self,name,val):
        self[name] = val

class Interpolator(object):
    def __init__(self,data):
        self.data = np.array(data)
        self.fn = interp1d(self.data[:,0],self.data[:,1],assume_sorted=False,kind='linear',fill_value='extrapolate')
    def __call__(self,P):
        return self.fn(P)
    @property
    def x(self):
        return self.data[:,0]
    @property
    def y(self):
        return self.data[:,1]

def Q(d1,d2):
    # d1,d2 collision diameters
    return (np.pi/16)*(d1+d2)**2

def binarydiffusion(T,m1,m2,d1,d2):
    '''
    Hard-sphere approximation for binary diffusion coefficient [1/cm/s]

    input parameters: 
        T - Temperature [k]
        m1/2 - mass of constituent 1/2 [g]
        d1/2 - diameter of constituent 1/2 [cm]
    '''
    q = Q(d1,d2)
    b = 3./(64*q) * ( 2*np.pi*kb*T*(m1+m2)/(m1*m2) )**0.5
    return b

kinetic_radius = Interpolator( np.array([atom_data['mass'],atom_data['radius']]).T )

def crossover_zero(m2,*args):
    '''
        Transcendtal form of cross over mass with 
        variable diffusion coefficient that depends on
        the cross over mass

        input parameters:
            m2 - crossover mass [amu]    
    '''
    Feuv,R,M,T,m1,x1 = args 

    # get the van der waals radius for mass of atom
    d1 = 2*kinetic_radius(m1)*1e-8 # convert Angstrom to cm
    d2 = 2*kinetic_radius(m2*mp)*1e-8 
    b1 = binarydiffusion(T,m1,m2,d1,d2)

    return crossover(*args,b1)/mp - m2

def massloss(Feuv,m1,M,R):
    '''
    input parameters:
        Feuv - Stellar Flux [erg/cm2/s]
        m1 - mass of constituent 1 being carried away [g]
        M - mass of planet [g]
        R - radius of planet [cm]
        returns [g/s]?
    '''
    gp = G*M/R # grav potential
    # pi r^2 being absorbed but 4*pi^2 being emitted 
    return (1./4) * Feuv * 1./(gp*m1)

def crossover(Feuv,R,M,T,m1,x1,b1):
    '''
    Crossover mass at which a heavier second constituent will not be dragged along with constituent 1
    See Hunten 1987 

    input parameters:
        Feuv - energy flux loss from planet (ratio of incident flux to grav potential) [erg/cm2/s]
        m1 - mass of constituent 1 being carried away [g]
        M - mass of planet [g]
        R - radius of planet [cm]
        T - Temperature [K]
        b1 - binary diffusion coefficient [1/cm/s]
        x1 - mole fraction of constituent 1
    '''
    mass = m1 + 0.25*Feuv*(R**3)*kb*T / ( (G*M)**2 *m1 *b1 *x1 )
    return mass

if __name__ == "__main__":
    
    # parameters for planet 
    pars = AttributeDict({
        'R':0.3767*rjup, # Radius of Planet [cm]
        'M':0.0727*mjup, # Mass of Planet [g]
        'a':0.0386*AU, # Semi-major Axis [AU]
        'm1':1*mp, # Mass of Constituent 1 [g]
        'starspec':'Spectra/hlsp_muscles_multi_multi_gj436_broadband_v22_const-res-sed.fits',
        'b1':lambda T:  1.04*10**18 * T**0.732, # binary diffusion coefficient [1/cm/s]
        'x1':0.9, # Fraction of Constituent 1 
        'T':600, # Temperature of Planet [K]
        'distance':10.14 # Distance to planet [parsec]
    })

    # H-He diffusion coefficients are different by a factor of 3
    # the paramterized function is larger

    # load data from star into grid 
    spec = fits.getdata(pars.starspec,1)
    mask = spec['WAVELENGTH'] < 912 # EUV wavelength
    Feuv = spec['FLUX'][mask] * (pars.distance*pc)**2/pars.a**2  # * spec['WAVELENGTH'][mask] * 1e-8 / (h*c)
    pars.Feuv = trapz(Feuv,spec['WAVELENGTH'][mask])
    pars.FLUX = spec['FLUX'][mask]
    pars.WAVELENGTH = spec['WAVELENGTH'][mask]

    cmass = findzero(crossover_zero, 1, 100, args=(pars.Feuv,
                                                    pars.R,
                                                    pars.M,
                                                    pars.T,
                                                    pars.m1,
                                                    pars.x1) )
    print('crossover mass: {:.1f} amu'.format(cmass) )    
