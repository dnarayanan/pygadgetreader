try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(name='PyGadgetReader',
      version='1.0',
      description='module to read all sorts of gadget files',
      author='Robert Thompson',
      author_email='rthompsonj@gmail.com',
      url='something.com',
      packages=['readgadget','readgadget.modules'],
      #scripts=['bin/test.py'],
      install_requires=['h5py>=2.2.1'],)
