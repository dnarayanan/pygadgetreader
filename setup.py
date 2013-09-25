#!/usr/bin/env python

from numpy.distutils.core import setup, Extension

setup(name = "readgadget",version="0.5",ext_modules=[Extension("readgadget",["readgadget.c"])])
