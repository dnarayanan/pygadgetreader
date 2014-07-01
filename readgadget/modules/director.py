
import hdf5 as hdf5
import tipsy as t
import gadget as g
import numpy as np

def readposvel(f,h,p):
    if h.hdf5_file:
        arr = hdf5.hdf5_readposvel(f,h,p)
    elif h.tipsy_file:
        arr = t.tipsy_read(f,h,p)
    elif h.gadget_file:
        arr = g.gadget_readposvel(f,h,p)
    return arr.astype(np.float64)

def readpid(f,h,p):
    if h.hdf5_file:
        arr = hdf5.hdf5_readpid(f,h,p)
    elif h.tipsy_file:
        arr = t.tipsy_read(f,h,p)
    elif h.gadget_file:
        arr = g.gadget_readpid(f,h,p)
    return arr

def readmass(f,h,p):
    if h.hdf5_file:
        arr = hdf5.hdf5_readmass(f,h,p)
    if h.tipsy_file:
        arr = t.tipsy_read(f,h,p)
    elif h.gadget_file:
        arr = g.gadget_readmass(f,h,p)
    return arr.astype(np.float64)

def readgasprop(f,h):
    if h.hdf5_file:
        arr = hdf5.hdf5_readgasprop(f,h)
    elif h.tipsy_file:
        arr = t.tipsy_read(f,h,0)
    elif h.gadget_file:
        arr = g.gadget_readgasprop(f,h)
    return arr.astype(np.float64)

def readmetals(f,h,p,single=1):
    if h.hdf5_file:
        arr = hdf5.hdf5_readmetals(f,h,p,single=single)
    elif h.tipsy_file:
        arr = t.tipsy_read(f,h,p)
    elif h.gadget_file:
        arr = g.gadget_readmetals(f,h,p,single=single)
    return arr.astype(np.float64)

def readpotentials(f,h,p):
    if h.hdf5_file:
        arr = hdf5.hdf5_readpotentials(f,h,p)
    elif h.tipsy_file:
        arr = t.tipsy_read(f,h,p)
    elif h.gadget_file:
        arr = g.gadget_readpotentials(f,h,p)
    return arr.astype(np.float64)

def readage(f,h):
    if h.hdf5_file:
        arr = hdf5.hdf5_readage(f,h)
    elif h.tipsy_file:
        arr = t.tipsy_read(f,h,4)
    elif h.gadget_file:
        arr = g.gadget_readage(f,h)
    return arr.astype(np.float64)
