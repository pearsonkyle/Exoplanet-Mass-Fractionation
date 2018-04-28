from scipy.integrate import trapz
from astropy.io import fits
import pandas as pd
import numpy as np

mp = 1.67e-24 # g
rjup = 6.9911e9 # cm 
mjup = 1.898e30 # g
rearth = 6.371e8 # cm
mearth = 5.97e27 # g
AU = 1.496e+13 # cm 
pc = 3.086e+18 # cm 

# database from exoplanets.org 
exoplanets = pd.read_csv('exoplanets.csv',low_memory=False)

class AttributeDict(dict):
    def __getattr__(self, name):
        return self[name]
    def __setattr__(self,name,val):
        self[name] = val

def get_pars(name='GJ 436 b', m1=1*mp, x1=0.9, T=750, starspec=''):

    # get star name     
    sname = ''.join(name.split(' ')[:-1]).lower()
    spectrafile = 'Spectra/hlsp_muscles_multi_multi_{}_broadband_v22_const-res-sed.fits'.format(sname)

    # get entry from exoplanets.org database 
    entry = exoplanets[ exoplanets.EANAME == name]
    assert(entry.shape[0]>=1)

    # TODO equilibrium temperature calculation 
    keys = ['MASS','UMASS','R','UR','DIST','UDIST','A','UA']
    mod = [mjup,mjup,rjup,rjup,pc,pc,AU,AU] # conversion to cgs

    pars = AttributeDict({
            'name':entry['EANAME'].item(),
            'm1':m1, # Mass of Constituent 1 [g]
            'starspec':spectrafile,
            'x1':x1, # Fraction of Constituent 1 
            'T':T # Temperature of Planet [K]
         })
    
    for i in range(len(keys)):
        k = keys[i]
        pars[k] = entry[k].item() * mod[i]
    

    # compute Feuv and Leuv from spectrum 
    spec = fits.getdata(pars.starspec,1)
    mask = (spec['WAVELENGTH'] < 912) #& (spec['WAVELENGTH'] > 100) # EUV wavelength
    Feuv = spec['FLUX'][mask] * pars.DIST**2/pars.A**2 
    pars.Feuv = trapz(Feuv,spec['WAVELENGTH'][mask]) # erg/cm2/s
    
    Leuv = spec['FLUX'][mask] * 4 * np.pi * pars.DIST**2
    pars.Leuv = trapz(Leuv,spec['WAVELENGTH'][mask]) # erg/s

    pars.FLUX = spec['FLUX'][mask]
    pars.WAVELENGTH = spec['WAVELENGTH'][mask]

    return pars 


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    transit_list = [
        'GJ 436 b',
        'GJ 1214 b',
        'HD 97658 b'
    ]

    # CREATE LATEX TABLE 
    # TODO equilibrium temperature calculation 
    keys = ['MASS','UMASS','R','UR','DIST','UDIST','A','UA']
    mod = [mjup,mjup,rjup,rjup,pc,pc,AU,AU] # conversion to cgs

    for i in range(len(transit_list)):
        pars = get_pars(transit_list[i])
        print("{} & {:.2f}$\pm${:.2f} & {:.2f}$\pm${:.2f} & {:.4f}$\pm${:.4f}& {:.2f}$\pm${:.2f} & {:.1f} ".format(pars.name,
                                            pars.R/rearth,pars.UR/rearth,
                                            pars.MASS/mearth,pars.UMASS/mearth,
                                            pars.A/AU,pars.UA/AU,
                                            pars.DIST/pc,pars.UDIST/pc,
                                            pars.T) )

    rv_list = [
        'GJ 176 b',

        'GJ 667 C b',
        'GJ 667 C c',

        'GJ 832 b',
        'GJ 832 c',

        'GJ 876 b',
        'GJ 876 c',
        'GJ 876 d',
        'GJ 876 e',

        'HD 40307 b',
        'HD 40307 c',
        'HD 40307 d',
        ]
