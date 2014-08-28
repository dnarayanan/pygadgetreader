import sys
import struct
import numpy as np
import gadget as g

TYPE2NAMES = {
    'pos':'POS ',
    'vel':'VEL ',
    'pid':'ID  ',
    'mass':'MASS',
    'u':'U   ',
    'rho':'RHO ',
    'ne':'NE  ',
    'nh':'NH  ',
    'hsml':'HSML',
    'sfr':'SFR',
    'metallicity':'Z   ',
    'pot':'POT ',
    'fh2':'fH2 ',
    'sigma':'Sigm',
    'age':'AGE '
}

def findBlock(f,h):
    RETURN = False

    blk1 = np.fromfile(f,dtype=np.uint32,count=1)
    if not blk1:
        return 2
    blk1 = blk1[0]
    NAME = struct.unpack('4s',f.read(4))[0]
    sl   = np.fromfile(f,dtype=np.uint32,count=1)
    blk2 = np.fromfile(f,dtype=np.uint32,count=1)[0]
    g.errorcheck(blk1,blk2,'blk %s' % NAME)

    ## custom return
    if h.reading not in TYPE2NAMES:
        if NAME == h.reading:
            RETURN = True
    elif NAME == TYPE2NAMES[h.reading]:
        RETURN = True

    if RETURN:
        if h.debug: print 'returning for block %s' % NAME
        return 1
    else:
        if h.debug: print 'skipping block %s' % NAME

    s1 = np.fromfile(f,dtype=np.uint32,count=1)[0]
    f.seek(s1, 1)
    s2 = np.fromfile(f,dtype=np.uint32,count=1)[0]
    g.errorcheck(s1,s2,NAME)
    return 0

"""
def gadget_general(f,h,ptype):
    skip1 = np.fromfile(f,dtype=np.uint32,count=1)[0]
    for i in range(0,ptype):
        f.seek(4 * h.npart[i],1)
    vals = np.fromfile(f,dtype=np.float32,count=h.npart[ptype])
    for i in range(ptype+1,len(h.npart)):
        f.seek(4 * h.npart[i],1)
    skip2 = np.fromfile(f,dtype=np.uint32,count=1)[0]
    g.errorcheck(skip1,skip2,"generic")
    return vals
"""

def gadget_general(f,h,ptype):
    skip1 = np.fromfile(f,dtype=np.uint32,count=1)[0]
    
    gOnly   = 0.
    sOnly   = 0.
    gsOnly  = 0.
    allPart = 0.

    if h.npart[0] > 0:
        gOnly   = skip1 / (4. * h.npart[0])
    if h.npart[4] > 0:
        sOnly   = skip1 / (4. * h.npart[4])
    if h.npart[0] > 0 or h.npart[4] > 0:
        gsOnly  = skip1 / (4. * h.npart[0] + 4. * h.npart[1])
    allPart = skip1 / (4. * np.sum(h.npart)) 

    if gOnly == 1.0:
        if ptype != 0:
            print 'block is only present for gas!'
            return
        vals = np.fromfile(f,dtype=np.float32,count=h.npart[ptype])
    elif sOnly == 1.0:
        if ptype != 4:
            print 'block is only present for stars!'
            return
        vals = np.fromfile(f,dtype=np.float32,count=h.npart[ptype])
    elif gsOnly == 1.0:
        if ptype != 0 and ptype != 4:
            print 'block is only present for gas & stars!'
            return
        if ptype == 4:
            f.seek(4 * h.npart[0],1)
        vals = np.fromfile(f,dtype=np.float32,count=h.npart[ptype])
        if ptype == 0:
            f.seek(4 * h.npart[4],1)
    elif allPart == 1.0:
        for i in range(0,ptype):
            f.seek(4 * h.npart[i],1)
        vals = np.fromfile(f,dtype=np.float32,count=h.npart[ptype])
        for i in range(ptype+1,len(h.npart)):
            f.seek(4 * h.npart[i],1)
        
    skip2 = np.fromfile(f,dtype=np.uint32,count=1)[0]
    g.errorcheck(skip1,skip2,"generic read")
    return vals
    
def gadget_type2_read(f,h,p):
    """Main driver for reading gadget type-2 binaries"""

    stat = 0
    while stat == 0:
        stat = findBlock(f,h)
        
    if stat == 2:
        print 'end of file =/'
        print 'scanning did not find block %s' % (h.reading)
        sys.exit()
    
    if h.reading == 'pos' or h.reading == 'vel':
        arr = g.gadget_readposvel(f,h,p)
    elif h.reading == 'pid':
        arr = g.gadget_readpid(f,h,p)
    elif h.reading == 'mass':
        arr = g.gadget_readmass(f,h,p)
    elif h.reading == 'metallicity':
        arr = g.gadget_readmetals(f,h,p,single=0)
    else:
        arr = gadget_general(f,h,p)
        
    return arr
