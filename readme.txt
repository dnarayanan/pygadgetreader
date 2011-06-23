This module will read in POS,VEL,PID,MASS,U,RHO,NE,NH,HSML,SFR,STELLARAGE,Z from gadget type 1 snapshot files.  It will read in single or multiple files.  There is also a nifty toggle that allows you to return data in code units, or physical units.  There are two main modules described below.

REQUIREMENTS:
 -python
 -numpy
 -c compiler

INSTALL:
python setup.py build     --> this builds the module
python setup.py install   --> this installs the module, may require sudo

MODULES:  
(in order to use these modules one must put 'from readgadget import *' at the top of your .py script)
##################################################

readhead - This module simply reads in the gadget header file and returns the value of interest. 
           The values it can read are:

	 time	       - scale factor of the snapshot
	 redshift      - redshift of the snapshot
	 boxsize       - boxsize if present in units of kpc/h
	 O0	       - Omega_0 (Omega_dm+Omega_m)
	 Ol	       - Omega_Lambda
	 h	       - hubble parameter
	 f_sfr	       - Star Formation Rate flag	 0=off 1=on
	 f_fb	       - Feedback flag	     		 0=off 1=on
	 f_cooling     - Cooling flag			 0=off 1=on 
	 f_age	       - Stellar Age tracking flag	 0=off 1=on
	 f_metals      - Metal tracking flag  		 0=off 1=on


Definition:   readhead('a','b',numfiles=0)

	      Parameters
	      ----------
	      a : Input file.  
	      	  Must be input as a string and enclosed in ' ' - see examples.
	      b : Value of interest from the above list.  
	      	  Must be input as a string and enclosed in ' ' - see examples.
	      
	      numfiles: Number of files the snapshot is broken up into. Assumed to be 1 if it is not included.


Example:
	      z=readhead('snap_001','redshift')  
	      		- reads redshift value and assigns it to the z variable
	      
	      h=readhead('snap_005','h',numfiles=2) 
	      boxsize=(readhead('snap_005','boxsize',numfiles=2)/1000/h)**3  
	      		- reads hubble param and assigns it to the h variable, then reads in the boxsize, 
			  converts it to Mpc^3.

##################################################

readsnap - This module does the bulk of the work.  It reads data blocks from the snapshot 
	   file and returns the requested data for specified particle type.  
	   Supported data blocks are:

	   pos	       - Position data
	   vel	       - Velocity data in km/s
	   pid	       - Particle ids
	   mass	       - Particle masses
	   u	       - Gas internal energy
	   rho	       - Gas density
	   ne	       - Number density of free electrons
	   nh	       - Number density of neutral hydrogen
	   hsml	       - Smoothing length of SPH particles
	   sfr	       - Gas star formation rate in Msun/year
	   age	       - Formation time of stellar particles (in terms of the scale factor)
	   z	       - Metallicty of gas & star particles


	   Supported particle types:
	   gas	       - Gas
	   dm	       - Dark Matter
	   disk	       - Disk particles
	   bulge       - Bulge particles
	   star/stars  - Star particles
	   bndry       - Boundary particles
	   

Definition:	readsnap('a','b','c',numfiles=0,units=0)

		Parameters
		----------
		a: Input file.
		   Must be input as a string and enclosed in ' ' - see examples.
		b: Data block you are interested in (see above list)
		   Must be input as a string and enclosed in ' ' - see examples.
		c: Particle type you are interested in (see above list)
		   Must be input as a string and enclosed in ' ' - see examples.
		
		numfiles: Number of files the snapshot is broken up into. Assumed to be 1 if it is not included.
		   units: Can either be 0 or 1.  Assumed to be 0 if not included.  
		          This parameter allows for the data to be returned in cgs units rather than code units.
			  Currently only active for density (rho) and internal energy 
			  (u - returns temperature in K).

Example:
		DMpos=readsnap('snap_001','pos','dm')
		DMx,DMy,DMz=hsplit(DMpos,3)
			- reads in dark matter data and returns an Nx3 array containing positions.  
			  hsplit is then used to split the array into x,y,z positions.

		 grho=readsnap('snap_005','rho','gas',numfiles=2,units=1)
		gtemp=readsnap('snap_005','u','gas',numfiles=2,units=1)
			- reads a multi-file snapshot (2) and returns density and temperature in cgs units.

##################################################

This code is heavily inspired by gadgetpyio written by Matthew Becker, Matthew Turk, & Peter Teuben - thanks for putting the source to your code online.

Any comments or suggestsions feel free to contact me:
Robert Thompson
UNLV Physics & Astronomy
rthompson@physics.unlv.edu

