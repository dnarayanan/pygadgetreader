This module contains a function for reading in data from Gadget type 1 snapshot files, reading in Gadget type 1 header information, and one for reading data from the PStarGroupFinder properties file (galprop).  These functions are described below:

REQUIREMENTS:
 -python
 -numpy
 -c compiler

INSTALL:
python setup.py build     --> this builds the module
python setup.py install   --> this installs the module, may require sudo

FUNCTIONS:  
(in order to use these functions one must put 'from readgadget import *' at the top of your .py script to import the module)
##################################################

readhead - This function simply reads in the gadget header file and returns the value of interest. 
           The values it can read are:

	 time	       - scale factor of the snapshot
	 redshift      - redshift of the snapshot
	 boxsize       - boxsize if present in units of kpc/h
	 O0	       - Omega_0 (Omega_dm+Omega_m)
	 Ol	       - Omega_Lambda
	 h	       - hubble parameter
	 gascount      - total # of gas   particles [type 0]
	 dmcount       - total # of DM    particles [type 1]
	 diskcount     - total # of disk  particles [type 2]
	 bulgecount    - total # of bulge particles [type 3]
	 starcount     - total # of star  particles [type 4]
	 bndrycount    - total # of bndry particles [type 5]
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
	      
	      Optional
	      --------
	      numfiles: DEPRECIATED!  No longer needed, but still here for backwards compatability with previous scripts:
	      Number of files the snapshot is broken up into. Assumed to be 1 if it is not included.


Example:
	      z=readhead('snap_001','redshift')  
	      		- reads redshift value and assigns it to the z variable
	      
	      h=readhead('snap_005','h',numfiles=2) 
	      boxsize=(readhead('snap_005','boxsize',numfiles=2)/1000/h)**3  
	      		- reads hubble param and assigns it to the h variable, then reads in the boxsize, 
			  converts it to Mpc^3.

##################################################

readsnap - This function does the bulk of the work.  It reads data blocks from the snapshot 
	   file and returns the requested data for specified particle type.  
	   Supported data blocks are:

	   --GADGET & TIPSY--
	   pos	       - (all)         Position data
	   vel	       - (all)         Velocity data in km/s
	   pid	       - (all)         Particle ids
	   mass	       - (all)         Particle masses
	   u	       - (gas)         Internal energy
	   rho	       - (gas)         Density
	   ne	       - (gas)         Number density of free electrons
	   nh	       - (gas)         Number density of neutral hydrogen
	   hsml	       - (gas)         Smoothing length of SPH particles
	   sfr	       - (gas)         Star formation rate in Msun/year
	   delaytime   - (gas)         DelayTime (>0 member of wind)
	   fH2	       - (gas)         Fractional Abundance of molecular hydrogen
	   Sigma       - (gas)         Approximate surface density @ each particle (HIdensity * scale_height)
	   age	       - (stars)       Formation time of stellar particles (in terms of the scale factor)
	   z	       - (gas & stars) Metallicty of gas & star particles (returns total Z)
	   tmax        - (gas & stars) Maximum temp
	   nspawn      - (gas & stars) Number of star particles spawned
	   potential   - (all)         Potential of particles

	   metals      - (gas & stars) NMETALS array [C,O,Si,Fe]

	   --TIPSY--
	   s_age       - (aux)    stellar age

	   mhalo       - (envira) mass of parent SKID halo
	   windage     - (envira) time since last launched in a wind for gas
	   rvir	       - (envira) ??
	   vvir	       - (envira) ??

	   starfrac    - (future) ??
	   relaunch    - (future) ??
	   rho	       - (future) gad density in the future
	   u	       - (future) gas particle temp in the future
	   z	       - (future) metal value of particle in the future


	   Supported particle types (note tipsy only returns gas/dm/star particles):
	   gas	       - Gas
	   dm	       - Dark Matter
	   disk	       - Disk particles
	   bulge       - Bulge particles
	   star/stars  - Star particles
	   bndry       - Boundary particles
	   

Definition:	readsnap('a','b','c',numfiles=0,units=0,tipsy=0,future=0)

		Parameters
		----------
		a: Input file.
		   Must be input as a string and enclosed in ' ' - see examples.
		b: Data block you are interested in (see above list)
		   Must be input as a string and enclosed in ' ' - see examples.
		c: Particle type you are interested in (see above list)
		   Must be input as a string and enclosed in ' ' - see examples.
		
		Optional
		--------
		numfiles: Number of files the snapshot is broken up into. Assumed to be 1 if it is not included.
		   units: Can either be 0 for code units or 1 for real units.  Assumed to be 0 if not included.  
		          This parameter allows for the data to be returned in real units(1) rather than code units(0).
			  Currently only active for density (rho), internal energy 
			  (u - returns temperature in K), Mass (returns Msun), and Sigma (returns g/cm^2).
		   tipsy: Must be set to 1 if reading in tipsy binary/aux/envira files.
		  future: integer value for the future file.
			  

Example:
		DMpos=readsnap('snap_001','pos','dm')
		DMx,DMy,DMz=hsplit(DMpos,3)
			- reads in dark matter data and returns an Nx3 array containing positions.  
			  hsplit is then used to split the array into x,y,z positions.

		 grho=readsnap('snap_005','rho','gas',numfiles=2,units=1)
		gtemp=readsnap('snap_005','u','gas',numfiles=2,units=1)
			- reads a multi-file snapshot (2) and returns density and temperature in cgs units.

##################################################

galprop	 - This function will read in the properties file outputted by P-Star group finder.  It then returns one of the following:

	 mstar	       - Mass of the stars within a group
	 bmag	       - B-magnitude of each group
	 imag	       - I-magnitude of each group
	 vmag	       - V-magnitude of each group
	 kmag	       - K-magnitude of each group
	 cm	       - Center of mass positions of each group
	 sfr	       - Instantanious star formation rate of each group
	 mgas	       - Mass of the gas within a group
	 zstar	       - Stellar metallicity of each group
	 zgas	       - Gas metallicity of each group


Definitions:	galprop('a',b,'c',units=0)

		Parameters
		----------
		a: Input directory (location of the property file)
		   Must be input as a string and enclosed in ' ' - see examples.
		b: Snapshot number
		c: Value of interest (see above list)
		   Must be input as a string and enclosed in ' ' - see examples.

		Optional
		--------
		untits: Can either be 0 for code units or 1 for real units.  Assumed to be 0 if not included.
		        This is only active for Gas & Star Masses (1 returns units of Msun)

Example:	group_mstar=galprop('../', 7, 'mstar',units=1)
			- returns the total stellar mass of each group in units of Msun

		group_cm=galprop('../', 7, 'cm')
		gx,gy,gz=hsplit(group_cm,3)
			- reads in group center of mass positions and returns an NgroupX3 array.
			  hsplit is then used to split the array into gx,gy,gz positions.

##################################################

galdata	 - This function will return data from all baryonic particle contained within a specified galprop group.  The data returned is an array with the following particle properties at each index:

	 xpos	       - x-position of the particle (index 0)
	 ypos	       - y-position of the particle (index 1)
	 zpos	       - z-position of the particle (index 2)
	 INDEX	       - index of the particle (in relation to ALL particles from the simulation) (index 3)
	 TYPE	       - TYPE of particle: 0=gas, 4=star (index 4)


Definitions:	galdata('a','b',c,d)

		Parameters
		----------
		a: Snapshot file.
		   Must be input as a string and enclosed in ' ' - see examples.
		b: Input directory (location of the property files)
		   Must be input as a string and enclosed in ' ' - see examples.
		c: Snapshot number
		d: Galaxy number of interest - see examples.

Example:	galaxy=galdata('snap_001', '../', 7, 5)
			- returns an array of data with the above information for galaxy 5 in the galprop list

To match data to the snapshot:
		indexes = galaxy[:,3]  			#unfortunately this returns an array of doubles
		indexes = indexes.astype(int)  	#convert that array to integers
		
		- now you can access the data for this particlular galaxy using these indexes on the full dataset
		  for example:  if you read in the density value via:	rho = readsnap(snap,'rho','gas')
		  				then rho[indexes] will return the density values of only the particles within
		  				the target galaxy.

##################################################

This code is heavily inspired by gadgetpyio written by Matthew Becker, Matthew Turk, & Peter Teuben - thanks for putting the source to your code online.

Any comments or suggestsions feel free to contact me:
Robert Thompson
UNLV Physics & Astronomy
rthompson@physics.unlv.edu

