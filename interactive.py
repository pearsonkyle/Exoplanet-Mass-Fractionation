from crossover_mass import *
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
import matplotlib.colors as colors

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

# load spectra and compute gravity
spec = fits.getdata(pars.starspec,1)
mask = spec['WAVELENGTH'] < 912 # EUV wavelength
Feuv = spec['FLUX'][mask] * (pars.distance*pc)**2/pars.a**2  # * spec['WAVELENGTH'][mask] * 1e-8 / (h*c)
pars.Feuv = trapz(Feuv,spec['WAVELENGTH'][mask])
pars.FLUX = spec['FLUX'][mask]
pars.WAVELENGTH = spec['WAVELENGTH'][mask]

f,axf = plt.subplots(1)
axf.plot(spec['WAVELENGTH'], spec['FLUX'],'k-')
axf.plot(spec['WAVELENGTH'][mask],Feuv,'r-')
axf.set_xlabel('Wavelength (Angstroms)')
axf.set_ylabel('Flux Density (erg/cm2/s/Ang)')


M = np.linspace( 0.5, 4, 90)*mearth
R = np.linspace( 0.5, 6, 40)*rearth
mg, rg = np.meshgrid(M,R)

T=500

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
                import pdb; pdb.set_trace()

                
mask = mg > 100
mg[mask] = 101

f,axg = plt.subplots(1)
cp = axg.contourf( M/mearth, R/rearth, mg, 40,cmap='nipy_spectral',vmin=1,vmax=100 )
axg.plot( pars.M/mearth, pars.R/rearth, 'ko')
axg.set_ylabel('Radius (Earth)')
axg.set_xlabel("Mass (Earth)")
cbar = f.colorbar(cp, ax=axg, extend='both')
cbar.ax.set_ylabel("Crossover mass (amu)")
plt.show() 
