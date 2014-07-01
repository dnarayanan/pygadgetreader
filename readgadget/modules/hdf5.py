import h5py as h5py
import numpy as np
from common import METALFACTOR

HDF5_NAMES = {'pos':'Coordinates',
              'vel':'Velocities',
              'pid':'ParticleIDs',
              'mass':'Masses',
              'u':'InternalEnergy',
              'rho':'Density',
              'hsml':'SmoothingLength',
              'ne':'ElectronAbundance',
              'nh':'NeutralHydrogenAbundance',
              'sfr':'StarForamtionRate',
              'metallicity':'Metallicity',
              'metalarray':'Metallicity',
              'age':'StellarFormationTime',
              'pot':'Potential',
              'fh2':'FractionH2',
              'sigma':'Sigma'}

def hdf5_readposvel(f,h,ptype): 
    arr = f['PartType%d/%s' % (ptype,HDF5_NAMES[h.reading])]
    return np.asarray(arr)*h.convert
    

def hdf5_readpid(f,h,ptype):    
    arr = f['PartType%d/%s' % (ptype,HDF5_NAMES[h.reading])]
    return np.asarray(arr)


def hdf5_readmass(f,h,ptype):
    ## check to see if mass block exists
    if HDF5_NAMES[h.reading] in f['PartType%d' % ptype]:
        arr = f['PartType%d/%s' % (ptype,HDF5_NAMES[h.reading])]
    else:
        arr = np.zeros(h.npart[ptype])
        arr.fill(h.mass[ptype])
    return np.asarray(arr)*h.convert


def hdf5_readgasprop(f,h):
    if HDF5_NAMES[h.reading] in f['PartType0']:
        arr = f['PartType0/%s' % (HDF5_NAMES[h.reading])]
    else:
        arr = np.zeros(h.npart[0])

    if h.units and h.reading == 'u':
        ne = f['PartType0/%s' % (HDF5_NAMES['ne'])]
        import common as common
        h.convert = common.getTfactor(np.asarray(ne),h)

    return np.asarray(arr)*h.convert


def hdf5_readmetals(f,h,ptype,single=1):
    if HDF5_NAMES[h.reading] in f['PartType%d' % ptype]:
        metals = f['PartType%d/%s' % (ptype,HDF5_NAMES[h.reading])]
        if single:
            arr = np.sum(metals,axis=1) * METALFACTOR
        else:
            arr = np.asarray(metals)
    else:
        if single:
            arr = np.zeros(h.npart[ptype])
        else:
            arr = np.zeros((h.npart[ptype],h.flag_metals))
    return arr


def hdf5_readpotentials(f,h,ptype):
    arr = f['PartType%d/%s' % (ptype,HDF5_NAMES[h.reading])]
    return np.asarray(arr)


def hdf5_readage(f,h):
    arr = f['PartType4/%s' % HDF5_NAMES[h.reading]]
    return np.asarray(arr)
