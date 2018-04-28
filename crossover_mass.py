from scipy.optimize import brentq as findzero
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy as np
import pickle

from get_parameters import get_pars 

kb = 1.38e-16  # g cm2 s-2 K-1 [erg/K]
G = 6.674e-8 # cm3 g-1 s-2
h = 6.626e-27 # [erg/s]
c = 2.9979e10 # cm s-1
mp = 1.67e-24 # g
# barye = g cm-1 s-2
rsun = 6.957e10 # cm 
rjup = 6.9911e9 # cm 
mjup = 1.898e30 # g
rearth = 6.371e8 # cm
mearth = 5.97e27 # g
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
    b1 = 2*binarydiffusion(T,m1,m2,d1,d2)

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
    check FEUV unit
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
    #pars = get_pars('GJ 436 b',T=750)
    #pars = get_pars('GJ 1214 b',T=400)
    pars = get_pars('HD 97658 b',T=750)
    
    cmass = findzero(crossover_zero, 1, 1000, args=(pars.Feuv,
                                                    pars.R,
                                                    pars.MASS,
                                                    pars.T,
                                                    pars.m1,
                                                    pars.x1) )
    
    upper = findzero(crossover_zero, 1, 1000, args=(pars.Feuv,
                                                    pars.R+pars.UR,
                                                    pars.MASS-pars.UMASS,
                                                    pars.T,
                                                    pars.m1,
                                                    pars.x1) )

    lower = findzero(crossover_zero, 1, 1000, args=(pars.Feuv,
                                                    pars.R-pars.UR,
                                                    pars.MASS+pars.UMASS,
                                                    pars.T,
                                                    pars.m1,
                                                    pars.x1) )

    print('crossover mass: {:.1f} - {:.1f} [amu]'.format(lower,upper) )    

    print(pars.Leuv)

    # atomic radii plot 
    #f,ax = plt.subplots(1)
    #ax.plot(atom_data['mass'],atom_data['radius'],'k-')
    #ax.set_xlabel("Mass (amu)")
    #ax.set_ylabel("Radius (A)")
    #ax.set_title("Atomic van der Waal Radii")
    #plt.show()