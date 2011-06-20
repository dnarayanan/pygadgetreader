This module will read in POS,VEL,PID,MASS,U,RHO,NE,NH,HSML,SFR,STELLARAGE,Z from gadget type 1 snapshot files.  It will read in single or multiple files.  So far it's all in code units so the user must convert them after an import.

This code is heavily inspired by gadgetpyio written by Matthew Becker, Matthew Turk, & Peter Teuben - thanks for putting the source to your code online.


REQUIREMENTS:
python
numpy
c compiler

INSTALL:
python setup.py build     --> this builds the module
python setup.py install   --> this installs the module, may require sudo

USAGE:
This module contains many different functions.  Each is called with this format -
function_name('data_file',num_of_source_files,'type')
Where type can be one of these: gas,dm,disk,bulge,stars,bndry
Function names are: readpos,readvel,readpid,readmass,readu,readrho,readNE,readNH,readHSML,readSFR,readage,readZ.  You can get a better description of these functions in ipython by typing the function name followed by a question mark (ie: readpid?)

USAGE FROM WITHIN PYTHON:
example:

	from readgadget import *
	vel=readvel('/home/user/snap_HALO_000',2,'gas')
	vx,vy,vz=hsplit(vel,3)

This 3 line code will read gas velocities from a 2 part snapshot (snap_HALO_000.0,snap_HALO_000.1) and return an NgasX3 array containing vx,vy,vz.  The call to hsplit just splits the array vel into 3 columns.




Any comments or suggestsions feel free to contact me:
Robert Thompson
UNLV Physics & Astronomy
rthompson@physics.unlv.edu

