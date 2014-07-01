from modules.common import *
import modules.header as HEAD
from modules.director import *
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

        if h.gadget_file:
            g.skipblocks(f,h,d)
        if d == 'pos' or d == 'vel':
            arr = readposvel(f,h,p)
        elif d == 'pid':
            arr = readpid(f,h,p)
        elif d == 'mass':
            arr = readmass(f,h,p)
        elif d in GasProps:
            if p != 0: 
                print 'WARNING!! you requested ParticleType%d for %s, returning GAS instead' % (p,d)
            arr = readgasprop(f,h)
        elif d == 'metallicity':
            arr = readmetals(f,h,p)
        elif d == 'metalarray':
            arr = readmetals(f,h,p,single=0)
        elif d == 'pot':
            arr = readpotentials(f,h,p)
        elif d == 'age':
            arr = readage(f,h)
        else:
            print 'no %s block!!! (requested %s)' % (d,data)
            sys.exit()

        f.close()

        if i > 0:
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
                        printer = '%s %s' % (printer,defaultDataUnits[d])
                    else:
                        printer = '%s in code units' % printer
            
                print printer

    return return_arr


