from util import *
from crossover_mass import get_pars

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    transit_list = [
        'HD 97658 b',
        'GJ 436 b',
        'GJ 1214 b',
    ]

    for i in range(len(transit_list)):
        pars = get_pars(transit_list[i])
        print("{} & {:.2f}$\pm${:.2f} & {:.2f}$\pm${:.2f} & {:.4f}$\pm${:.4f}& {:.2f}$\pm${:.2f} & {:.1f} ".format(pars.name,
                                            pars.R/rearth,pars.UR/rearth,
                                            pars.MASS/mearth,pars.UMASS/mearth,
                                            pars.A/AU,pars.UA/AU,
                                            pars.DIST/pc,pars.UDIST/pc,
                                            pars.T) )
        print("{} - {:.2e} g/s".format(pars.name,pars.massloss))

    # TODO: generate spectra plot 
    
    '''
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
    '''