from scipy.optimize import brentq as findzero
from scipy.interpolate import interp1d
from scipy.integrate import trapz
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import pickle

from util import *
kinetic_radius = Interpolator( np.array([atom_data['mass'],atom_data['radius']]).T )

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

def crossover_zero(m2,*args):
    '''
        Transcendtal form of cross over mass with 
        variable diffusion coefficient that depends on
        the cross over mass

        input parameters:
            m2 - crossover mass [amu]    
    '''
    T,Mdot,M,R,x1,m1, = args 

    # get the van der waals radius for mass of atom
    d1 = 2*kinetic_radius(m1)*1e-8 # convert Angstrom to cm
    d2 = 2*kinetic_radius(m2*mp)*1e-8 
    b1 = 2*binarydiffusion(T,m1,m2,d1,d2)

    return crossover(*args,b1)/mp - m2

def massloss(Mdot,R):
    '''
    input parameters:
        Mdot - Mass loss flux from planet [g/cm2/s]
        R - radius of planet [cm]

    returns: [g/s] 
    '''
    return Mdot*4*np.pi*R**2

def massflux(Leuv,A,MASS,R, ruvrp2,K,eta, **kwargs):
    '''
    Mass loss flux from the planet due to incident stellar radiation

    input parameters:
        Leuv - Luminosity of star [erg/s]
        A - semi major axis of planet [cm]
        MASS - mass of planet [g]
        R - radius of planet [cm]

        # correction factors 
        ruvrp2 - scaling of absorption radius to emitted (R_uv/R_p)^2
        K - Roche lobe effect 
        eta - efficiency of heating 
    '''
    Feuv = eta*Leuv/(4*np.pi*A**2) # portion of flux used for escape [erg/cm2/s]
    grav = K*G*MASS/R # graviational potential 
    Mdot = Feuv * 0.25 * ruvrp2 / grav 
    return Mdot

def crossover(T,Mdot,M,R,x1,m1, b1):
    '''
    Crossover mass at which a heavier second constituent will not be dragged along with constituent 1
    See Hunten 1987 

    input parameters:
    check FEUV unit
        Feuv - mass loss flux from energy ratio of incident stellar flux to grav potential [g/cm2/s]
        m1 - mass of constituent 1 being carried away [g]
        M - mass of planet [g]
        R - radius of planet [cm]
        T - Temperature [K]
        b1 - binary diffusion coefficient [1/cm/s]
        x1 - mole fraction of constituent 1

        fudge factors
        eta - efficiency of heating 
        K - Roche lobe effect (Erkaev+2007)
        a - scaling of absorption radius to emitted 
    '''
    m2 = m1 + kb*T*Mdot*R**2/(b1*G*M*x1*m1)
    return m2

def equilibrium_temp(TEFF,AR,F,Ab,**kwargs):
    '''
    Equilbirum temperature (Southworth 2010)
    input parameters:
        Teff - effective temperature of star [K]
        AR - scaled semi-major axis of planet [a/R*]
        Ab - Bond Albedo
        F - Heat Redistribution factor
            0.25 = heating averaged over planet, 0.5 just heating on dayside
    '''
    Teq = TEFF * ( (1-Ab)/(4*F) )**0.25 * (1./AR)**0.5 * (1/2)**0.5
    return Teq 



def get_pars(name='GJ 436 b', m1=1*mp, x1=0.9, T=0, heatefficiency=0.2, spectrafile='', Ab=0, F=0.25):
    '''
    get parameters to perform the calculation for crossover mass 
    will retrieve exoplanet parameters from exoplanets.org csv file 

    input parameters:
        m1 - mass of constituent 1 being carried away [g]
        T - temperature, if 0 the equilibrium temperature will be calculated [K] 
        Ab - bond albedo 
    '''

    # get star name     
    sname = ''.join(name.split(' ')[:-1]).lower()
    if spectrafile == '':
        spectrafile = 'Spectra/hlsp_muscles_multi_multi_{}_broadband_v22_const-res-sed.fits'.format(sname)

    # get entry from exoplanets.org database 
    planet = Exoplanet(name)

    # only get keys of interest from exoplanets.org
    keys = ['MASS','UMASS','R','UR','DIST','UDIST','A','UA','TEFF','AR','MSTAR','UMSTAR']
    mod = [mjup,mjup,rjup,rjup,pc,pc,AU,AU,1,1,msun,msun] # conversion to cgs

    pars = AttributeDict({
            'name':planet['EANAME'],
            'm1':m1, # Mass of Constituent 1 [g]
            'starspec':spectrafile,
            'x1':x1, # Fraction of Constituent 1 
            'T':T # Temperature of Planet [K]
         })
    
    for i in range(len(keys)):
        k = keys[i]
        pars[k] = planet[k] * mod[i]
    
    # compute Feuv and Leuv from spectrum 
    spec = fits.getdata(pars.starspec,1)
    mask = (spec['WAVELENGTH'] < 912) #& (spec['WAVELENGTH'] > 100) # EUV wavelength
    Feuv = spec['FLUX'][mask] * pars.DIST**2/pars.A**2 # erg/cm2/s/A
    pars.Feuv = trapz(Feuv,spec['WAVELENGTH'][mask]) # erg/cm2/s
    pars.Leuv = 4*np.pi*pars.DIST**2 * trapz(spec['FLUX'][mask], spec['WAVELENGTH'][mask]) # erg/s
    pars.FLUX = spec['FLUX'][mask]
    pars.WAVELENGTH = spec['WAVELENGTH'][mask]

    # heating efficiency/fraction of absorbed EUV energy that drives escape 
    pars.eta = heatefficiency

    # equilibrium temperature (A=0,F=0.25)
    pars.Teq = equilibrium_temp(**pars,F=F,Ab=Ab)
    if T!=0: pars.Teq = T # override equilibrium temp calc

    # compute radius scale factor 
    # HI cross section between 13-20eV: ~2-7e-18 cm^2
    # tau=1 pressure of X-ray/EUV absorption between 13-20eV Verner1996
    pars.gravity = G*pars.MASS/pars.R**2
    P = mp*pars.gravity / 5e-18
    pars.scalefactor = -1*np.log(P/10) # scale heights above 10 bar
    pars.scaleheight = kb*pars.Teq / (2.3*mp*pars.gravity)
    
    # ratio between the planet's EUV absorbing radius and planetary radius
    pars.ruvrp2 = ((pars.scaleheight*pars.scalefactor + pars.R) / pars.R)**2
    
    # potential energy reduction factor due to roche lobe effect
    delta = pars.MASS / pars.MSTAR
    lam = pars.A / pars.R
    eta = lam * (delta/3)**(1./3)
    pars.K = 1-3/(2*eta) + 1./(2*eta**3)

    # TODO 
    # solve for sonic point height
    # compute transonic heating rate from Johnson 2013

    # compute mass loss rate 
    pars.Mdot = massflux(**pars)
    pars.massloss = massloss(pars.Mdot, pars.R)

    return pars 


if __name__ == "__main__":
    
    # parameters for planet 
    pars = get_pars('GJ 436 b')
    #pars = get_pars('GJ 1214 b')
    #pars = get_pars('HD 97658 b')

    cmass = findzero(crossover_zero, 1, 100, args=(pars.Teq,
                                                    pars.Mdot,
                                                    pars.MASS,
                                                    pars.R,
                                                    pars.x1,
                                                    pars.m1) )

    print('cross over mass: {:.1f}'.format(cmass) )

    # Atmospheric hydrodynamic escape calculation for the crossover mass at which a heavier second constituent will not be dragged along with a lighter constituent (i.e. H).