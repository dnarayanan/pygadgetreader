import numpy as np
from collections import OrderedDict
import sys

METALFACTOR = 0.0189/0.0147
H_MASSFRAC  = 0.76
BOLTZMANN   = 1.3806e-16
PROTONMASS  = 1.6726e-24
GAMMA       = 5.0 / 3.0

UnitLength_in_cm         = 3.085678e21
UnitMass_in_g            = 1.989e43
UnitVelocity_in_cm_per_s = 1.0e5

#############################################
## for gadget type-1 binary block ordering ##
#############################################
## you can easily add your own block ordering!  
# duplicate one of the below, then modify it to
# match your sims.  Don't forget to add yours
# to the dict BLOCKORDERING.
#
# Format of block ordering looks like this:
#
# (BlockName, [ParticleTypes,FlagsToCheck])
#
# -> BlockName : string
#    must be equal to those in dataTypes value
#    (value!!, not key, see below)
# -> ParticleTypes : integer list
#    defines what types of particles have this block
#    -1 means ALL particles
# -> FlagsToCheck : string
#    which flags (if any) to check that determine if block
#    is present
#############################################

## sets the default ordering
DEFAULT_BLOCKORDERING = 'romeel'

## Romeel's block ordering
BLOCKORDERING0 = OrderedDict([
    ('pos' ,[-1]),
    ('vel' ,[-1]),
    ('pid' ,[-1]),
    ('mass',[-1]),
    ('u'   ,[0]),
    ('rho' ,[0]),
    ('ne'  ,[0]),
    ('nh'  ,[0]),
    ('hsml',[0]),
    ('sfr' ,[0]),
    ('delaytime',  [0,'flag_delaytime']),
    ('fh2'  ,      [ 0,'flag_fh2']),
    ('sigma',      [ 0,'flag_fh2']),
    ('age' ,       [ 4,'flag_age']),
    ('metallicity',[[0,4],'flag_metals']),
    ('tmax',       [[0,4],'flag_tmax']),
    ('nspawn',     [[0,4]]),
    ('pot',        [-1,'flag_potential'])
])

## Ken's block ordering
BLOCKORDERING1 = OrderedDict([
    ('pos' ,[-1]),
    ('vel' ,[-1]),
    ('pid' ,[-1]),
    ('mass',[-1]),
    ('u'   ,[0]),
    ('rho' ,[0]),
    ('ne'  ,[0]),
    ('nh'  ,[0]),
    ('hsml',[0]),
    ('sfr' ,[0]),
    ('age' ,       [ 4,'flag_age']),
    ('metallicity',[[0,4],'flag_metals']),
    ('fh2'  ,      [ 0,'flag_fh2']),
    ('sigma',      [ 0,'flag_fh2']),
    ('pot',        [-1,'flag_potential'])
])

## NAME THE BLOCK ORDERINGS ##
BLOCKORDERING = {'romeel':BLOCKORDERING0,
                 'ken'   :BLOCKORDERING1}

## default particle types
pTypes = {0:0,1:1,2:2,3:3,4:4,5:5,
          'gas':0,'dm':1,'disk':2,'bulge':3,'stars':4,'bndry':5}
# account for different names
pTypes['star'] = 4

pNames = {0:'GAS',1:'DM',2:'DISK',3:'BULGE',4:'STAR',5:'BNDRY'}

## default header returns (corresponds to h.vals[KEY])
headerTypes = {'npart':'npart',
               'ngas':'ngas',
               'ndm':'ndm',
               'ndisk':'ndisk',
               'nbulge':'nbulge',
               'nstar':'nstar',
               'nbndry':'nbndry',
               'mass':'mass',
               'time':'time',
               'nfiles':'nfiles',
               'redshift':'redshift',
               'boxsize':'boxsize',
               'O0':'O0',
               'Ol':'Ol',
               'h':'h',
               'flag_sfr':'flag_sfr',
               'flag_cooling':'flag_cooling',
               'flag_sfr':'flag_sfr',
               'flag_fb':'flag_fb',
               'flag_fh2':'flag_fh2',
               'flag_age':'flag_age',
               'flag_metals':'flag_metals',
               'flag_delaytime':'flag_delaytime',
               'flag_tmax':'flag_tmax',
               'flag_potential':'flag_potential'}
# account for different names
headerTypes['num_files']     = 'nfiles'
headerTypes['nstars']        = 'nstar'
headerTypes['gascount']      = 'ngas'
headerTypes['dmcount']       = 'ndm'
headerTypes['diskcount']     = 'ndisk'
headerTypes['bulgecount']    = 'nbulge'
headerTypes['bndrycount']    = 'nbndry'
headerTypes['starcount']     = 'nstar'
headerTypes['a']             = 'time'
headerTypes['z']             = 'redshift'
headerTypes['box']           = 'boxsize'
headerTypes['Omega0']        = 'O0'
headerTypes['OmegaLambda']   = 'Ol'
headerTypes['hubble']        = 'h'
headerTypes['flag_feedback'] = 'flag_fb'
headerTypes['f_sfr']         = 'flag_sfr'
headerTypes['f_fb']          = 'flag_fb'
headerTypes['f_cooling']     = 'flag_cooling'
headerTypes['f_age']         = 'flag_age'
headerTypes['f_fh2']         = 'flag_fh2'
headerTypes['f_metals']      = 'flag_metals'


## default data types
# the VALUE here is the important part
dataTypes = {'pos':'pos',
             'vel':'vel',
             'pid':'pid',
             'mass':'mass',
             'u':'u',
             'rho':'rho',
             'ne':'ne',
             'nh':'nh',
             'hsml':'hsml',
             'sfr':'sfr',
             'delaytime':'delaytime',
             'fh2':'fh2',
             'sigma':'sigma',
             'age':'age',
             'z':'metallicity',
             'zarray':'metalarray',
             'tmax':'tmax',
             'nspawn':'nspawn',
             'pot':'pot'}
# account for different names
dataTypes['positions']  = 'pos'
dataTypes['vels']       = 'vel'
dataTypes['velocity']   = 'vel'
dataTypes['velocities'] = 'vel'
dataTypes['fH2']        = 'fh2'
dataTypes['FH2']        = 'fh2'
dataTypes['metals']     = 'metalarray'

## values used for output logging
dataNames = {'pos':'Positions',
             'vel':'Velocities',
             'pid':'Particle IDs',
             'mass':'Mass',
             'u':'Internal Energy',
             'rho':'Density',
             'ne':'Electron Abundance',
             'nh':'Neutral Hydrogen Density',
             'hsml':'Smoothing Length',
             'sfr':'Star Formation Rate',
             'delaytime':'Delay Time',
             'fh2':'Fractional H2 abundance',
             'sigma':'Surface Density',
             'age':'Stellar Age',
             'metallicity':'Metallicity',
             'metalarray':'Metal Array',
             'tmax':'Maximum Temperature',
             'nspawn':'Number of Stars Spawned',
             'pot':'Potential'}
dataDefaultUnits = {'sfr':'[Msun/yr]'}
dataUnits = {'vel':'[km/s, peculiar]',
             'mass':'[Msun/h]',
             'u':'[Kelvin]',
             'rho':'[h^2 g/cm^3, physical]',
             'sfr':'[Msun/yr]',
             'sigma':'[h g/cm^2, physical]'}


## properties that redirect to readgasprops()
GasProps     = ['u','rho','ne','nh','hsml','sfr',
                'delaytime','fh2','sigma']
GasStarProps = ['tmax','nspawn']


RecognizedOptions = ['units','hdf5','tipsy','supress_output','blockordering','debug','double']
def pollOptions(KWARGS,data,ptype):
    """warn user if option is unrecognized"""
    for key,items in KWARGS.iteritems():
        if key not in RecognizedOptions:
            print 'WARNING!! option not recognized: %s' % key

    kill = 0
    if data not in dataTypes:
        print 'ERROR! %s not a recognized data request' % data
        kill = 1
    if ptype not in pTypes:
        print 'ERROR! %s not a recognized particle type' % ptype
        kill = 1
    if kill:
        sys.exit()


def initUnits(h):
    """initialize conversion factors"""
    convert = 1.0

    if h.units and not h.tipsy_file:
        if h.reading == 'rho':
            if h.boxsize > 0. and h.OmegaLambda > 0:
                convert = ((1.0 + h.redshift)**3 * 
                           h.UnitMass_in_g / h.UnitLength_in_cm**3)
            else:
                convert = h.UnitMass_in_g / h.UnitLength_in_cm**3
        elif h.reading == 'vel':
            if h.boxsize > 0. and h.Ol > 0:
                convert = np.sqrt(h.time)
        elif h.reading == 'sigma':
            convert = h.UnitMass_in_g / h.UnitLength_in_cm**2
        elif h.reading == 'u':
            h.UnitTime_in_s = h.UnitLength_in_cm / h.UnitVelocity_in_cm_per_s
            h.UnitEnergy_in_cgs = (h.UnitMass_in_g * h.UnitLength_in_cm**2 / 
                                   h.UnitTime_in_s**2)
        elif h.reading == 'mass':
            convert = h.UnitMass_in_g
    h.convert = convert


def getTfactor(Ne,h):
    """calculate temperature conversion factor including Ne"""
    MeanWeight = (4.0 / (3.0 * H_MASSFRAC + 1.0 + 4.0 * H_MASSFRAC * Ne) * 
                  PROTONMASS)
    conversion = (MeanWeight / BOLTZMANN * (GAMMA - 1.0) * 
                  h.UnitEnergy_in_cgs / h.UnitMass_in_g)
    return conversion
