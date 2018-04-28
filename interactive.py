import matplotlib.pyplot as plt
import matplotlib.colors as colors
from astropy.io import fits
import numpy as np


from get_parameters import get_pars
from crossover_mass import *

# parameters for planet 
#pars = get_pars('GJ 436 b',T=750)
pars = get_pars('GJ 1214 b',T=400)
pars = get_pars('')


# load spectra and compute gravity
spec = fits.getdata(pars.starspec,1)
mask = spec['WAVELENGTH'] < 912 # EUV wavelength
Feuv = spec['FLUX'][mask] * (pars.DIST*pc)**2/pars.A**2

f,axf = plt.subplots(1)
axf.plot(spec['WAVELENGTH'], spec['FLUX']* (pars.DIST*pc)**2/pars.A**2,'k-')
axf.plot(spec['WAVELENGTH'][mask],Feuv,'r-')
axf.set_xlabel('Wavelength (Angstroms)')
axf.set_ylabel('Flux Density (erg/cm2/s/Ang)')
axf.set_title(pars.name)
a = plt.axes([.65, .6, .2, .2])

a.plot(spec['WAVELENGTH'], spec['FLUX']* (pars.DIST*pc)**2/pars.A**2,'k-')
a.plot(spec['WAVELENGTH'][mask],Feuv,'r-')
a.set_xlim([-20,1200])
a.set_ylim([0,3])
a.set_xlabel('Wavelength (A)')
a.set_ylabel('Flux Density')


M = np.linspace( 1, 25, 90)*mearth
R = np.linspace( 1, 4.5, 40)*rearth
mg, rg = np.meshgrid(M,R)

# TODO TEST THIS
for i in range(R.shape[0]):
    for j in range(M.shape[0]):
                
        try:
                cmass = findzero(crossover_zero, 1, 10000, args=(pars.Feuv,
                                                        R[i],
                                                        M[j],
                                                        pars.T,
                                                        pars.m1,
                                                        pars.x1) )
                mg[i,j] = cmass
        except:
                mg[i,j] = 101
                #import pdb; pdb.set_trace()

                
mask = mg > 100
mg[mask] = 101

f,axg = plt.subplots(1)
cp = axg.contourf( M/mearth, R/rearth, mg, 20,cmap='viridis',vmin=1,vmax=100 )
levels = [2,4,8,18,28,44]
CS = plt.contour(mg, levels,
                 origin='lower',
                 linewidths=2,
                 vmin=1,vmax=100,
                 extent=(min(M/mearth), max(M/mearth), min(R/rearth), max(R/rearth)),
                 colors='w' )
plt.clabel(CS,inline=1, fontsize=10, fmt='%1.1f')

#axg.plot( pars.M/mearth, pars.R/rearth, 'ko')
axg.errorbar(pars.MASS/mearth, pars.R/rearth,yerr=pars.UR/rearth,xerr=pars.UMASS/mearth,fmt='o',zorder=10,color='black')
axg.set_title(pars.name)
axg.set_ylabel('Radius (Earth)')
axg.set_xlabel("Mass (Earth)")
cbar = f.colorbar(cp, ax=axg, extend='both')
cbar.ax.set_ylabel("Crossover mass (amu)")
plt.show() 
