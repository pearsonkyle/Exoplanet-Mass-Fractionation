# Exoplanet-Mass-Fractionation
Atmospheric hydrodynamic escape calculation for the crossover mass at which a heavier second constituent will not be dragged along with a lighter constituent (i.e. H).

The summary of this project is as follows: 
- Use Navier-Stokes Equation including Drag forces to model atmospheric escape from an exoplanet atmosphere
- Derive the energy limited escape rate from the planet using EUV spectra from MUSCLES Treasury Survey
- Derive diffusion coefficients using kinetic theory and van der Waal radius of atoms
- Compute the cossover mass (see Hunten 1987 for derivation) 

This project was completed for class PTYS 598B - Special Topics in Planetary Science (Atmospheric Escape) taught by Roger Yelle. 


## File Guide
[mass-fractionation-exoplanet.pdf](https://github.com/pearsonkyle/Exoplanet-Mass-Fractionation/blob/master/mass-fractionation-exoplanet.pdf) - White paper for the project.

[crossover_mass.py](https://github.com/pearsonkyle/Exoplanet-Mass-Fractionation/blob/master/crossover_mass.py) - calculations for energy limited escape, mass los rates and crossover mass

[atomic van der waal radii.csv](https://github.com/pearsonkyle/Exoplanet-Mass-Fractionation/blob/master/atomic%20van%20der%20waal%20radii.csv) - Van der Waal radii of atoms used when computing diffusion coefficients

[exoplanets.csv](https://github.com/pearsonkyle/Exoplanet-Mass-Fractionation/blob/master/exoplanets.csv) - Exoplanet parameters from exoplanets.org 

[util.py](https://github.com/pearsonkyle/Exoplanet-Mass-Fractionation/blob/master/util.py) - Useful constants and data structures 
