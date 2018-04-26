# Exoplanet-Mass-Fractionation
Atmospheric hydrodynamic escape calculation for the crossover mass at which a heavier second constituent will not be dragged along with a lighter constituent (i.e. H).


The summary of this project is as follows: 
- Use Navier-Stokes Equation including Drag forces to model atmospheric escape from an exoplanet atmosphere
- Derive the energy limited escape rate from the planet using EUV spectra from MUSCLES Treasury Survey
- Derive diffusion coefficients using kinetic theory and van der Waal radius of atoms
- Compute the cossover mass (see Hunten 1987 for derivation) 

This project was completed for class PTYS 598B - Special Topics in Planetary Science (Atmospheric Escape) taught by Roger Yelle

## File Guide
[escape.py](https://github.com/pearsonkyle/Exoplanet-Mass-Fractionation/blob/master/escape.py) - Basic functions to compute crossover mass, energy limited escape and diffusion coefficients 

[Exoplanet Mass Fractionation.ipynb](https://github.com/pearsonkyle/Exoplanet-Mass-Fractionation/blob/master/Exoplanet%20Mass%20Fractionation.ipynb) - Interactive plots for crossover mass that let you adjust the planetary parameters. 

[atomic van der waal radii.csv](https://github.com/pearsonkyle/Exoplanet-Mass-Fractionation/blob/master/atomic%20van%20der%20waal%20radii.csv) - Van der Waal radii of atoms used when computing diffusion coefficients

[interactive.py](https://github.com/pearsonkyle/Exoplanet-Mass-Fractionation/blob/master/interactive.py) - same program as the jupyter notebook 
