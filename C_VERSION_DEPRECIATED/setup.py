#!/usr/bin/env python

## Default data type to read, select 1 ##
GADGET_DEFAULT = 1
TIPSY_DEFAULT  = 0
HDF5_DEFAULT   = 0

## HDF5 DIRECTORIES ##
HDF5INCL = "/Users/bob/local/hdf5/include"
HDF5LIB  = "/Users/bob/local/hdf5/lib"

## ALTERNATE BLOCK STRUCTURE FOR GADGET BINARIES ##
ALTBLOCK = 0


##################
## NO TOUCHING! ##
##################
from numpy.distutils.core import setup, Extension
import os,sys

NAME = "readgadget"
FILE = "readgadget.c"

if GADGET_DEFAULT + TIPSY_DEFAULT + HDF5_DEFAULT > 1:
    print 'cannot select more than one file format!'
    sys.exit()
if GADGET_DEFAULT + TIPSY_DEFAULT + HDF5_DEFAULT == 0:
    print 'no file selected, defaulting to GADGET'

INCLDIRS,LIBDIRS,LIBS,MACROS = [],[],[],[]

## check for HDF5
HAVE_HDF5 = 0
if os.path.isdir(HDF5INCL) and os.path.isdir(HDF5LIB):
    HAVE_HDF5 = 1
    INCLDIRS.append(HDF5INCL)
    LIBDIRS.append(HDF5LIB)
    LIBS.append('hdf5')
    MACROS.append(('HAVE_HDF5',HAVE_HDF5))
else:
    if HDF5_DEFAULT:
        print 'could not find HDF5 dirs!'
        sys.exit()

if HDF5_DEFAULT:
    MACROS.append(('HDF5_DEFAULT',HDF5_DEFAULT))
elif TIPSY_DEFAULT:
    MACROS.append(('TIPSY_DEFAULT',TIPSY_DEFAULT))

if ALTBLOCK:
    MACROS.append(('ALTBLOCK',1))


setup(name = "readgadget",
      version="0.5",
      ext_modules=[Extension("readgadget",["readgadget.c"],
                             define_macros = MACROS,
                             include_dirs = INCLDIRS,
                             libraries    = LIBS,
                             library_dirs = LIBDIRS)])
