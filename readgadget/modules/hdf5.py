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


def hdf5_general(f,h,ptype):

    if ('PartType%d' % ptype) in f:
        if HDF5_NAMES[h.reading] in f['PartType%d' % ptype]:
            arr = f['PartType%d/%s' % (ptype,HDF5_NAMES[h.reading])]
        else:
            arr = np.zeros(h.npart[ptype])

        if h.units and h.reading == 'u':
            ne = f['PartType0/%s' % (HDF5_NAMES['ne'])]
            import common as common
            h.convert = common.getTfactor(np.asarray(ne),h)
    else:
        arr = np.zeros(0)

    return np.asarray(arr)*h.convert


def hdf5_readmass(f,h,ptype):

    if ('PartType%d' % ptype) in f:
        ## check to see if mass block exists
        if HDF5_NAMES[h.reading] in f['PartType%d' % ptype]:
            arr = f['PartType%d/%s' % (ptype,HDF5_NAMES[h.reading])]
        else:
            arr = np.zeros(h.npart[ptype])
            arr.fill(h.mass[ptype])
    else:
        arr = np.zeros(0)

    return np.asarray(arr)*h.convert


def hdf5_readmetals(f,h,ptype,single=1):
    
    if ('PartType%d' % ptype) in f:
        if HDF5_NAMES[h.reading] in f['PartType%d' % ptype]:
            metals = f['PartType%d/%s' % (ptype,HDF5_NAMES[h.reading])]
            if single and h.flag_metals > 1:
                arr = np.sum(metals,axis=1) * METALFACTOR
            else:
                arr = np.asarray(metals)
        else:
            if single:
                arr = np.zeros(h.npart[ptype])
            else:
                arr = np.zeros((h.npart[ptype],h.flag_metals))
    else:
        if single:
            arr = np.zeros(0)
        else:
            arr = np.zeros((0,h.flag_metals))

    return arr


def hdf5_read(f,h,p):
    """Main driver for reading HDF5 files"""
    if h.reading == 'metallicity':
        arr = hdf5_readmetals(f,h,p)
    elif h.reading == 'metalarray':
        arr = hdf5_readmetals(f,h,p,single=0)
    elif h.reading == 'mass':
        arr = hdf5_readmass(f,h,p)
    else:
        arr = hdf5_general(f,h,p)

    if h.reading != 'pid':
        arr = arr.astype(np.float64)

    return arr
