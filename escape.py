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

class AttributeDict(dict):
    def __getattr__(self, name):
        return self[name]
    def __setattr__(self,name,val):
        self[name] = val

class grid(object):
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

def binaryparam(T,m1,m2,d1,d2):
    '''
    Hard-sphere approximation for binary diffusion coefficient 

    input parameters: 
        T - Temperature [k]
        m1/2 - mass of constituent 1/2 [g]
        d1/2 - diameter of constituent 1/2 [cm]
    '''
    q = Q(d1,d2)
    b = 3./(64*q) * ( 2*np.pi*kb*T*(m1+m2)/(m1*m2) )**0.5
    return b

def crossover_eval(m2,*args):
    Feuv,R,M,T,x1 = args 
    b1 = binaryparam(T,m1,m2,d1,d2)
    value = m1 + 0.25*Feuv*(R**3)*kb*T / ( (G*M)**2 *m1 *b1 *x1 ) - m2
    return value 

# 2 grids, one more crossover mass, the other for collision diameters
# def photorate(opt_ni, *args) 
#findzero( crossover_eval, 0, 40, args=(pars,) )

def massloss(Feuv,m1,M,R):
    '''
    input parameters:
        Feuv - Stellar Flux [erg/cm2/s]
        m1 - mass of constituent 1 being carried away [g]
        M - mass of planet [g]
        R - radius of planet [cm]
    '''
    gp = G*M/R # grav potential
    # pi r^2 being absorbed but 4*pi^2 being emitted 
    return (1./4) * Feuv * 1./(gp*m1)

def crossover(Feuv,R,M,T,m1,b1,x1):
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
        'radius':0.3767*rjup, # cm
        'distance':0.02872, # AU
        'mass':0.0727*mjup, # g
        'm1':1*mp, # g
        'starspec':'hlsp_muscles_multi_multi_gj436_broadband_v22_const-res-sed.fits',
        'b1':lambda T:  1.04*10**18 * T**0.732, 
        'X1':0.9,
        'T':600, 
    })

    # load data from star into grid 
    spec = fits.getdata(pars.starspec,1)
    mask = spec['WAVELENGTH'] < 912
    Feuv = spec['FLUX'][mask] * spec['WAVELENGTH'][mask] * 1e-8 / (h*c) / pars.distance**2
    pars.Feuv = trapz(Feuv,spec['WAVELENGTH'][mask])
    pars.FLUX = spec['FLUX'][mask]
    pars.WAVELENGTH = spec['WAVELENGTH'][mask]

    pars.g = G*pars.mass/pars.radius**2
    
    pars.escape = escape(pars.Feuv,
                         pars.m1,
                         pars.mass,
                         pars.radius)

    pars.crossover = crossover(pars.m1,
                                pars.escape,
                                pars.T,
                                pars.b1(pars.T),
                                pars.g,
                                pars.X1)

    f,ax = plt.subplots(1)
    ax.plot(spec['WAVELENGTH'][mask],spec['FLUX'][mask])
    ax.set_xlabel('Wavelength (Angstroms)')
    ax.set_ylabel('Flux Density (erg/cm2/s/Ang)')
    plt.show()
    
