import numpy as np
import os,sys
import common as c
import struct

class Header(object):
    def __init__(self,snap,filenum, *args):
        self.snap_passed = snap
        self.args        = args[0]
        self.fileType    = ''

        #self.reading     = reading

        ## make these modifiable at some point
        self.UnitMass_in_g = c.UnitMass_in_g
        self.UnitLength_in_cm = c.UnitLength_in_cm
        self.UnitVelocity_in_cm_per_s = c.UnitVelocity_in_cm_per_s
        
        ## return double array? (instead of float)
        self.double = False
        if 'double' in args[0] and args[0]['double'] == 1:
            self.double = True

        ## debug?
        self.debug = False
        if 'debug' in args[0] and args[0]['debug'] == 1:
            self.debug = True

        ## supress output?
        self.supress = False
        if 'supress_output' in args[0] and args[0]['supress_output'] == 1:
            self.supress = True

        ## allow for different block orderings on the fly
        self.BLOCKORDER = c.BLOCKORDERING[c.DEFAULT_BLOCKORDERING]
        if 'blockordering' in args[0]:
            self.BLOCKORDER = c.BLOCKORDERING[args[0]['blockordering']]

        ## unit conversions?
        self.units = False
        if 'units' in args[0] and args[0]['units'] == 1:
            self.units = True

        ## make sure user has not selected tipsy & hdf5
        if ('hdf5' in args[0]  and args[0]['hdf5']  == 1 and
            'tipsy' in args[0] and args[0]['tipsy'] == 1):
            print 'cannot select tipsy and hdf5!!'
            sys.exit()

        ## detect file type
        self.snap = self.detectFileType(args,filenum)

        ## read header
        if self.fileType == 'gadget':
            s1,s2 = self.read_gadget_header()
            if s1 == 256 and s2 == 256:
                self.extension = ''
            else:
                print 'found GADGET file, but header is off...%d vs %d' % (s1,s2)
                sys.exit()
        elif self.fileType == 'hdf5':
            self.read_hdf5_header()
            self.extension = 'hdf5'
        elif self.fileType == 'tipsy':
            self.read_tipsy_header()
            self.extension = 'bin/aux'
        else:
            print 'should not be here...'
            sys.exit()

        if not self.supress or self.debug:
            tmptxt = 'reading %s' % self.snap
            if self.nfiles > 1:
                tmptxt = '%s (%d/%d files)' % (tmptxt,filenum+1,self.nfiles)
                print tmptxt
            elif self.debug:
                print tmptxt

        if not hasattr(self,'npartThis'):
            self.npartThis = []

        ## assign dictionary
        self.nparticles = np.sum(self.npart)
        self.vals = {'npart':self.npartTotal,
                     'npartThis':self.npartThis,
                     'ngas':self.npartTotal[0],
                     'ndm':self.npartTotal[1],
                     'ndisk':self.npartTotal[2],
                     'nbulge':self.npartTotal[3],
                     'nstar':self.npartTotal[4],
                     'nbndry':self.npartTotal[5],
                     'mass':self.mass,
                     'time':self.time,
                     'nfiles':self.nfiles,
                     'redshift':self.redshift,
                     'boxsize':self.boxsize,
                     'O0':self.Omega0,
                     'Ol':self.OmegaLambda,
                     'h':self.HubbleParam,
                     'flag_cooling':self.flag_cool,
                     'flag_sfr':self.flag_sfr,
                     'flag_fb':self.flag_fb,
                     'flag_fh2':self.flag_fH2,
                     'flag_age':self.flag_age,
                     'flag_metals':self.flag_metals,
                     'flag_potential':self.flag_potential,
                     'flag_delaytime':self.flag_delaytime,
                     'flag_tmax':self.flag_tmax}

    def read_gadget_header(self):
        """read gadget type-1 & 2 binary header"""
        import gadget as g
        f = open(self.snap,'rb')
        self.f = f

        ## test for type 2 binaries
        if np.fromfile(f,dtype=np.uint32,count=1)[0] == 8:
            NAME = struct.unpack('4s',f.read(4))[0]
            np.fromfile(f,dtype=np.uint32,count=1)
            np.fromfile(f,dtype=np.uint32,count=1)
            self.fileType = 'gadget2'
            if self.debug: print 'detected TYPE 2 gadget binary'
        else:
            f.seek(0)

        skip1 = g.skip(f)


        self.npart       = np.fromfile(f,dtype=np.uint32,count=6)
        self.mass        = np.fromfile(f,dtype=np.float64,count=6)
        self.time        = np.fromfile(f,dtype=np.float64,count=1)[0]
        self.redshift    = np.fromfile(f,dtype=np.float64,count=1)[0]
        self.flag_sfr    = np.fromfile(f,dtype=np.int32,count=1)[0]
        self.flag_fb     = np.fromfile(f,dtype=np.int32,count=1)[0]
        self.npartTotal  = np.fromfile(f,dtype=np.uint32,count=6)
        self.flag_cool   = np.fromfile(f,dtype=np.int32,count=1)[0]
        self.nfiles      = np.fromfile(f,dtype=np.int32,count=1)[0]
        self.boxsize     = np.fromfile(f,dtype=np.float64,count=1)[0]
        self.Omega0      = np.fromfile(f,dtype=np.float64,count=1)[0]
        self.OmegaLambda = np.fromfile(f,dtype=np.float64,count=1)[0]
        self.HubbleParam = np.fromfile(f,dtype=np.float64,count=1)[0]
        self.flag_age    = np.fromfile(f,dtype=np.int32,count=1)[0]
        self.flag_metals = np.fromfile(f,dtype=np.int32,count=1)[0]
        self.npartTotalHighWord   = np.fromfile(f,dtype=np.uint32,count=6)
        self.flag_entropy         = np.fromfile(f,dtype=np.int32,count=1)[0]
        self.flag_doubleprecision = np.fromfile(f,dtype=np.int32,count=1)[0]
        self.flag_potential       = np.fromfile(f,dtype=np.int32,count=1)[0]
        self.flag_fH2             = np.fromfile(f,dtype=np.int32,count=1)[0]
        self.flag_tmax            = np.fromfile(f,dtype=np.int32,count=1)[0]
        self.flag_delaytime       = np.fromfile(f,dtype=np.int32,count=1)[0]
        
        self.npartThis.append(self.npart)

        bl = 256
        if self.fileType == 'gadget2':
            bl += 4 + 4 + 4 + 4
        bytes_left = bl + 4 - f.tell()
        f.seek(bytes_left,1)
        
        skip2 = g.skip(f)
        g.errorcheck(skip1,skip2,'header')
        #if skip1 != 256 or skip2 != 256:
        #    print 'not a gadget header!! %d vs %d' % (skip1,skip2)
        return skip1,skip2
        
    def read_hdf5_header(self):
        """read HDF5 header"""
        import h5py as h5py
        f = h5py.File(self.snap,'r')
        self.f = f
        hd = f['Header']
        ha = hd.attrs

        if not hasattr(self,'npartThis'):
            self.npartThis = []


        self.npart       = ha['NumPart_ThisFile']
        self.mass        = ha['MassTable']
        self.time        = ha['Time']
        self.redshift    = ha['Redshift']
        self.flag_sfr    = ha['Flag_Sfr']
        self.flag_fb     = ha['Flag_Feedback']
        self.npartTotal  = ha['NumPart_Total']
        self.flag_cool   = ha['Flag_Cooling']
        self.nfiles      = ha['NumFilesPerSnapshot']
        self.boxsize     = ha['BoxSize']
        self.Omega0      = ha['Omega0']
        self.OmegaLambda = ha['OmegaLambda']
        self.HubbleParam = ha['HubbleParam']
        self.flag_age    = ha['Flag_StellarAge']
        self.flag_metals = ha['Flag_Metals']
        self.npartTotalHighWord   = ha['NumPart_Total_HighWord']
        self.flag_entropy         = 1
        self.flag_doubleprecision = ha['Flag_DoublePrecision']

        self.npartThis.append(self.npart)

        if 'PartType0/Potential' in f:
            self.flag_potential   = 1
        else:
            self.flag_potential   = 0
        if 'PartType0/FractionH2' in f:
            self.flag_fH2         = 1
        else:
            self.flag_fH2         = 0
        if 'PartType0/TemperatureMax' in f:
            self.flag_tmax        = 1
        else:
            self.flag_tmax        = 0
        if 'PartType0/DelayTime' in f:
            self.flag_delaytime   = 1
        else:
            self.flag_delaytime   = 0
            

    def read_tipsy_header(self):
        """read TISPY header"""
        f = open(self.snap,'rb')
        self.f    = f
        self.time = np.fromfile(f,dtype=np.float64,count=1)[0]
        ntotal    = np.fromfile(f,dtype=np.int32,count=1)[0]
        ndim      = np.fromfile(f,dtype=np.int32,count=1)[0]
        ngas      = np.fromfile(f,dtype=np.int32,count=1)[0]
        ndark     = np.fromfile(f,dtype=np.int32,count=1)[0]
        nstar     = np.fromfile(f,dtype=np.int32,count=1)[0]
        alignment = np.fromfile(f,dtype=np.float32,count=1)
        
        self.npart       = [ngas,ndark,0,0,nstar,0]
        self.npartTotal  = [ngas,ndark,0,0,nstar,0]
        self.mass        = [0.,0.,0.,0.,0.,0.]
        self.redshift    = 1.0 / self.time - 1.0
        self.flag_sfr    = 1
        self.flag_fb     = 1
        self.flag_cool   = 1
        self.nfiles      = 1
        self.boxsize     = 0.0
        self.Omega0      = 0.0
        self.OmegaLambda = 0.0
        self.HubbleParam = 0.0
        self.flag_age    = 1
        self.flag_metals = 1
        self.npartTotalHighWord = [0,0,0,0,0,0]
        self.flag_entropy = 0
        self.flag_doubleprecision = 0
        self.flag_potential = 1
        self.flag_fH2    = 0
        self.flag_tmax   = 1
        self.flag_delaytime = 1


    def detectFileType(self,args,filenum):
        
        ## first check to see if the user specified which type
        if 'hdf5' in args[0] and args[0]['hdf5'] == 1:
            snap = '%s.hdf5' % self.snap_passed
            if os.path.isfile(snap):
                self.fileType = 'hdf5'
                return snap
            else:
                snap = '%s.%d.hdf5' % (self.snap_passed,filenum)
                if os.path.isfile(snap):
                    self.fileType = 'hdf5'
                    return snap
                else:
                    print 'could not locate HDF5 file!'
                    sys.exit()
        elif 'tipsy' in args[0] and args[0]['tipsy'] == 1:
            snap = '%s.bin' % self.snap_passed
            if os.path.isfile(snap):
                self.fileType = 'tipsy'
                return snap
            else:
                print 'could not locate TIPSY file!'
                sys.exit()


        ## auto detection routine

        # GADGET BINARY
        snap = self.snap_passed
        if os.path.isfile(snap):
            if self.debug: print 'detected GADGET file'
            self.fileType = 'gadget'
            return snap
        else:
            snap = '%s.%d' % (self.snap_passed,filenum)
            if os.path.isfile(snap):
                if self.debug: print 'detected GADGET multi-part file'
                self.fileType = 'gadget'
                return snap

        # HDF5 
        snap = '%s.hdf5' % self.snap_passed
        if os.path.isfile(snap):
            if self.debug: print 'detected HDF5 file'
            self.fileType = 'hdf5'
            return snap
        else:
            snap = '%s.%d.hdf5' % (self.snap_passed,filenum)
            if os.path.isfile(snap):
                if self.debug: print 'detected HDF5 multi-part file'
                self.fileType = 'hdf5'
                return snap
            
        # TIPSY
        snap = '%s.bin' % self.snap_passed
        if os.path.isfile(snap):
            if self.debug: print 'detected TIPSY file'
            self.fileType = 'tipsy'
            return snap

        print 'COULD NOT FIND FILE %s' % self.snap_passed
        sys.exit()
            
