from modules.common import *
import modules.header as HEAD
import modules.gadget as gadget
import modules.tipsy as tipsy
import modules.hdf5 as hdf5
import numpy as np
import sys

def readhead(snap,data,**kwargs):
    """Read and return desired header info"""
    h = HEAD.Header(snap,data,0,kwargs)
    h.f.close()
    return h.vals[headerTypes[data]]


def readsnap(snap,data,ptype,**kwargs):
    """Read and return desired snapshot data"""
    pollOptions(kwargs,data,ptype)

    d = dataTypes[data]
    p = pTypes[ptype]

    h = HEAD.Header(snap,d,0,kwargs)
    f = h.f
    initUnits(h)
    
    #print 'reading %d files...' % (h.nfiles)

    for i in range(0,h.nfiles):
        if i > 0:
            h = HEAD.Header(snap,d,i,kwargs)
            f = h.f
            initUnits(h)

        if h.hdf5_file:
            arr = hdf5.hdf5_read(f,h,p)
        elif h.tipsy_file:
            arr = tipsy.tipsy_read(f,h,p)
        elif h.gadget_file:
            gadget.skipblocks(f,h,d)
            arr = gadget.gadget_read(f,h,p)

        f.close()

        if i > 0:
            if len(arr) > 0:
                return_arr = np.concatenate((return_arr,arr))
        else:
            return_arr = arr

            if not h.supress:
                ## print statement
                printer = 'Returning %s' % dataNames[d]
                if h.units:
                    if d in dataUnits:
                        if d == 'u':
                            printer = 'Returning Temperature %s' % (dataUnits[d])
                        else:
                            printer = '%s %s' % (printer,dataUnits[d])
                else:
                    if d in dataDefaultUnits:
                        printer = '%s %s' % (printer,dataDefaultUnits[d])
                    else:
                        printer = '%s in code units' % printer
            
                print printer

    return return_arr


