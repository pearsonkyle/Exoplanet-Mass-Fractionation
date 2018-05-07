# Exoplanet-Mass-Fractionation
Exoplanets with orbits close to their host star are subject to intense stellar radiation. High energy radiation from x-ray to ultraviolet drives photochemistry and atmospheric escape in the upper atmosphere of planets. If the escape rate of a lighter constituent (i.e. H) is high enough, it can drag heavier species along with it thus fractionating the atmosphere. Here I propose to use hydrodynamic escape in an energy limited regime to explore the effects of drag on carrying away heavier species (Z>1).

The summary of this project is as follows: 
- Use Navier-Stokes equation including drag forces to model atmospheric escape from an exoplanet atmosphere
- Derive the energy limited escape rate from the planet using EUV spectra from MUSCLES Treasury Survey
- Derive diffusion coefficients using kinetic theory and van der Waal radii of atoms
- Compute the crossover mass (see Hunten 1987 for derivation) 

This project was completed for class PTYS 598B - Special Topics in Planetary Science (Atmospheric Escape) taught by Roger Yelle (Spring2018). 


## File Guide
[mass-fractionation-exoplanet.pdf](https://github.com/pearsonkyle/Exoplanet-Mass-Fractionation/blob/master/mass-fractionation-exoplanet.pdf) - White paper for the project.

[crossover_mass.py](https://github.com/pearsonkyle/Exoplanet-Mass-Fractionation/blob/master/crossover_mass.py) - calculations for energy limited escape, mass los rates and crossover mass

[atomic van der waal radii.csv](https://github.com/pearsonkyle/Exoplanet-Mass-Fractionation/blob/master/atomic%20van%20der%20waal%20radii.csv) - Van der Waal radii of atoms used when computing diffusion coefficients

[exoplanets.csv](https://github.com/pearsonkyle/Exoplanet-Mass-Fractionation/blob/master/exoplanets.csv) - Exoplanet parameters from exoplanets.org 

[util.py](https://github.com/pearsonkyle/Exoplanet-Mass-Fractionation/blob/master/util.py) - Useful constants and data structures 
