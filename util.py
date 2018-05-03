from scipy.interpolate import interp1d
import pandas as pd

kb = 1.38e-16  # g cm2 s-2 K-1 [erg/K]
G = 6.674e-8 # cm3 g-1 s-2
h = 6.626e-27 # [erg/s]
c = 2.9979e10 # cm s-1
mp = 1.67e-24 # g
# barye = g cm-1 s-2
rsun = 6.957e10 # cm 
msun = 1.989e33 # g
rjup = 6.9911e9 # cm 
mjup = 1.898e30 # g
rearth = 6.371e8 # cm
mearth = 5.97e27 # g
AU = 1.496e+13 # cm 
pc = 3.086e+18 # cm 

# database from exoplanets.org 
exoplanets = pd.read_csv('exoplanets.csv',low_memory=False)

# load in atomic data for kinetic radius  
atom_data = pd.read_csv('atomic van der waal radii.csv')

class AttributeDict(dict):
    def __getattr__(self, name):
        return self[name]
    def __setattr__(self,name,val):
        self[name] = val

class Exoplanet(AttributeDict):
    def __init__(self,name):
        entry = exoplanets[ exoplanets.EANAME == name]
        assert(entry.shape[0]>=1)

        # load all keys/columns as a class property
        for k in exoplanets.keys():
            self[k] = entry[k].item()

class Interpolator(object):
    def __init__(self,data):
        self.data = data
        self.fn = interp1d(self.data[:,0],self.data[:,1],assume_sorted=False,kind='linear',fill_value='extrapolate')
    def __call__(self,P):
        return self.fn(P)
    @property
    def x(self):
        return self.data[:,0]
    @property
    def y(self):
        return self.data[:,1]