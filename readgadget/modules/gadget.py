import sys
import numpy as np
from common import METALFACTOR,headerTypes,GasProps,GasStarProps



def skip(f):
    skipval = np.fromfile(f,dtype=np.uint32,count=1)
    #print skipval
    return skipval[0]

def errorcheck(s1,s2,block):
    if s1!=s2:
        print 'wtf issue here - block %s >> %d vs %d' % (block,s1,s2)
        sys.exit()

def skipblocks(f,h,val):
    """Skip to desired block!"""
    def skiptypes(ptypes,block,multiplier=1):
        """Skip an entire block of given type"""
        if ptypes == -1:
            ptypes = [0,1,2,3,4,5]
        if isinstance(ptypes,int):
            ptypes = [ptypes]
        nparts = np.sum(h.npart[ptypes])
        if nparts == 0:
            return
            
        s1 = skip(f)
        for i in ptypes:
            f.seek(4 * multiplier * h.npart[i],1)
        s2 = skip(f)
        errorcheck(s1,s2,block)
            
    def skipmasses():
        """Skip mass block if it exists"""
        mb_exists = 0
        for i in range(0,len(h.npart)):
            if h.mass[i] == 0 and h.npart[i] > 0:
                mb_exists = 1
                break
        if mb_exists:
            s1 = skip(f)
            for i in range(0,len(h.npart)):
                if h.mass[i] == 0 and h.npart[i] > 0:
                    f.seek(4 * h.npart[i],1)
            s2 = skip(f)
            errorcheck(s1,s2,'mass')

    
    for key,items in h.BLOCKORDER.iteritems():
        if val == key: 
            if h.debug:
                print 'returning for key %s' % key
            return

        if val == 'metalarray' and key == 'metallicity': return

        multi = 1
        if key == 'pos' or key == 'vel':
            multi = 3
        if key == 'metallicity' or key == 'metalarray':
            multi = h.flag_metals

        if h.debug: print 'skipping %s' % key
        if key == 'mass':
            skipmasses()
        elif len(items) == 1:
            skiptypes(items[0],key,multiplier=multi)
        elif len(items) > 1:
            if h.vals[headerTypes[items[1]]]:
                skiptypes(items[0],key,multiplier=multi)
    return

    """
    if val == 'pos': return
    skiptypes(-1,'pos',multiplier=3)

    if val == 'vel': return
    skiptypes(-1,'vel',multiplier=3)

    if val == 'pid': return
    skiptypes(-1,'pid')

    if val == 'mass': return
    mb_exists = 0
    for i in range(0,len(h.npart)):
        if h.mass[i] == 0 and h.npart[i] > 0:
            mb_exists = 1
            break
    if mb_exists:
        s1 = skip(f)
        for i in range(0,len(h.npart)):
            if h.mass[i] == 0 and h.npart[i] > 0:
                f.seek(4 * h.npart[i],1)
        s2 = skip(f)
        errorcheck(s1,s2,'mass')
                
    if val == 'u': return
    skiptypes(0,'u')

    if val == 'rho':return
    skiptypes(0,'rho')

    if val == 'ne': return
    skiptypes(0,'ne')

    if val == 'nh': return
    skiptypes(0,'nh')
    
    if val == 'hsml': return
    skiptypes(0,'hsml')

    if val == 'sfr': return
    skiptypes(0,'sfr')
    
    if val == 'delaytime': return
    skiptypes(0,'delaytime')

    if h.flag_fH2:
        if val == 'fh2': return
        skiptypes(0,'fh2')
        if val == 'sigma': return
        skiptypes(0,'sigma')
    
    if h.flag_age:
        if val == 'age': return
        skiptypes(4,'age')
        
    if h.flag_metals:
        if val == 'metallicity' or val == 'metalarray': return
        skiptypes([0,4],'metal',multiplier=h.flag_metals)

    if h.flag_tmax:
        if val == 'tmax': return
        skiptypes([0,4],'tmax')

    if val == 'nspawn': return
    skiptypes([0,4],'nspawn')
    
    if val == 'pot': return
    skiptypes(-1,'pot')
    """
    


def gadget_readposvel(f,h,ptype):
    skip1 = skip(f)
    for i in range(0,ptype):
        f.seek(4*3*h.npart[i],1)
    posvel = np.fromfile(f,dtype=np.float32,count=h.npart[ptype]*3)
    for i in range(ptype+1,len(h.npart)):
        f.seek(4*3*h.npart[i],1)
    skip2 = skip(f)
    errorcheck(skip1,skip2,'posvel')
    
    posvel = posvel.reshape(h.npart[ptype],3)
    return posvel*h.convert

def gadget_readpid(f,h,ptype):
    skip1 = skip(f)
    for i in range(0,ptype):
        f.seek(4 * h.npart[i],1)
    pid = np.fromfile(f,dtype=np.uint32,count=h.npart[ptype])
    for i in range(ptype+1,len(h.npart)):
        f.seek(4 * h.npart[i],1)
    skip2 = skip(f)
    errorcheck(skip1,skip2,'pids')
    return pid

def gadget_readmass(f,h,ptype):
    mb_exists = 0
    for i in range(0,len(h.npart)):
        if h.mass[i] == 0 and h.npart[i] > 0:
            mb_exists = 1
            break

    if mb_exists:
        skip1 = skip(f)
        for i in range(0,ptype):
            if h.mass[i] == 0 and h.npart[i] > 0:
                f.seek(4 * h.npart[i],1)
                
        if h.mass[ptype] > 0:
            mass = np.zeros(h.npart[ptype],dtype=np.float32)
            mass.fill(h.mass[ptype])
        elif h.mass[ptype] == 0 and h.npart[ptype] > 0:
            mass = np.fromfile(f,dtype=np.float32,count=h.npart[ptype])

        for i in range(ptype+1,len(h.npart)):
            if h.mass[i] == 0 and h.npart[i] > 0:
                f.seek(4 * h.npart[i],1)

        skip2 = skip(f)
        errorcheck(skip1,skip2,'mass')
        return mass*h.convert
        #return mass.astype(np.float64)*h.convert
    else:
        mass = np.zeros(h.npart[ptype],dtype=np.float32)
        mass.fill(h.mass[ptype])
        return mass

def gadget_readgasprop(f,h):
    skip1   = skip(f)
    gasprop = np.fromfile(f,dtype=np.float32,count=h.npart[0])
    skip2   = skip(f)
    errorcheck(skip1,skip2,'gasprop')

    ## read in Ne for temp calc
    if h.units and h.reading == 'u':
        skip1 = skip(f)
        f.seek(4 * h.npart[0],1)
        skip2 = skip(f)
        errorcheck(skip1,skip2,'rho for Temp')
        
        skip1 = skip(f)
        ne    = np.fromfile(f,dtype=np.float32,count=h.npart[0])
        skip2 = skip(f)
        errorcheck(skip1,skip2,'ne for Temp')

        import common as common
        h.convert = common.getTfactor(ne,h)

    return gasprop*h.convert

def gadget_readgasstarprop(f,h,ptype):
    skip1 = skip(f)
    if ptype == 4:
        f.seek(4 * h.npart[0],1)
    gasstarprop = np.fromfile(f,dtype=np.float32,count=h.npart[ptype])
    if ptype == 0:
        f.seek(4 * h.npart[4],1)
    skip2 = skip(f)
    errorcheck(skip1,skip2,'gas-star prop')
    return gasstarprop*h.convert

def gadget_readmetals(f,h,ptype,single=1):
    skip1 = skip(f)
    if ptype == 4:
        f.seek(4 * h.flag_metals * h.npart[0],1)
    metals = np.fromfile(f,dtype=np.float32,count=h.flag_metals * h.npart[ptype])
    if ptype == 0:
        f.seek(4 * h.flag_metals * h.npart[4],1)
    skip2 = skip(f)
    errorcheck(skip1,skip2,'metals')

    if single and h.flag_metals > 1:
        newZ   = np.zeros(h.npart[ptype],dtype=np.float32)
        metals = metals.reshape(h.npart[ptype],h.flag_metals)
        for i in range(0,h.npart[ptype]):
            tmpz = 0.0
            for k in range(0,h.flag_metals):
                tmpz += metals[i,k]
            newZ[i] = tmpz * METALFACTOR
        metals = newZ
        newZ   = None
    elif h.flag_metals > 1:
        metals = metals.reshape(h.npart[ptype],h.flag_metals)

    return metals

def gadget_readpotentials(f,h,ptype):
    skip1 = skip(f)
    for i in range(0,ptype):
        f.seek(4 * h.npart[i],1)
    potentials = np.fromfile(f,dtype=np.float32,count=h.npart[ptype])
    for i in range(ptype+1,len(h.npart)):
        f.seek(4 * h.npart[i],1)
    skip2 = skip(f)
    errorcheck(skip1,skip2,'potentials')

    return potentials

def gadget_readage(f,h):
    skip1 = skip(f)
    age   = np.fromfile(f,dtype=np.float32,count=h.npart[4])
    skip2 = skip(f)
    errorcheck(skip1,skip2,'age')

    return age



def gadget_read(f,h,p,d):
    """Main driver for reading gadget binaries"""

    skipblocks(f,h,d)

    if h.reading == 'pos' or h.reading == 'vel':
        arr = gadget_readposvel(f,h,p)
    elif h.reading == 'pid':
        arr = gadget_readpid(f,h,p)
    elif h.reading == 'mass':
        arr = gadget_readmass(f,h,p)
    elif h.reading in GasProps:
        if p != 0: 
            print('WARNING!! you requested ParticleType%d for %s, returning GAS instead' 
                  % (p,h.reading))
        arr = gadget_readgasprop(f,h)
    elif h.reading in GasStarProps:
        arr = gadget_readgasstarprop(f,h,p)
    elif h.reading == 'metallicity':
        arr = gadget_readmetals(f,h,p)
    elif h.reading == 'metalarray':
        arr = gadget_readmetals(f,h,p,single=0)
    elif h.reading == 'pot':
        arr = gadget_readpotentials(f,h,p)
    elif h.reading == 'age':
        arr = gadget_readage(f,h)
    else:
        print 'no clue what to read =('
        arr = np.zeros(0)
    
    return arr
