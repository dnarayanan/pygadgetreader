#!/usr/bin/env python

from numpy.distutils.core import setup, Extension

setup(name = "readgadget",version="0.01",ext_modules=[Extension("readgadget",["readgadget.c"])])
