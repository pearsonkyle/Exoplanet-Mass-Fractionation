import matplotlib.pyplot as plt

from util import *
from crossover_mass import get_pars

if __name__ == "__main__":

    transit_list = [
        'HD 97658 b',
        'GJ 436 b',
        'GJ 1214 b',
    ]

    f,ax = plt.subplots(1)

    for i in range(len(transit_list)):
        pars = get_pars(transit_list[i])
        print("{} & {:.2f}$\pm${:.2f} & {:.2f}$\pm${:.2f} & {:.4f}$\pm${:.4f}& {:.2f}$\pm${:.2f} & {:.1f} & {:.1f} ".format(pars.name,
                                            pars.R/rearth,pars.UR/rearth,
                                            pars.MASS/mearth,pars.UMASS/mearth,
                                            pars.A/AU,pars.UA/AU,
                                            pars.DIST/pc,pars.UDIST/pc,
                                            pars.Teq,
                                            pars.TEFF) )
        print("{} - {:.2e} g/s".format(pars.name,pars.massloss))

        ax.semilogx(pars.WAVELENGTH,pars.FLUX* pars.DIST**2/pars.A**2 ,label=pars.name)

    ax.set_ylabel('Stellar Flux at Planet [erg/cm2/s/A]')
    ax.set_xlabel('Wavelength (nm)')
    ax.legend(loc='best')
    plt.show()

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