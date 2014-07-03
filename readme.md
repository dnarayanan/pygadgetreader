# pyGadgetReader

Author: Robert Thompson

E-Mail: rthompsonj@gmail.com

# Contents
* [Summary](#sum)
* [Requirements](#req)
* [Obtaining](#obt)
* [Customization](#cust)
* [Installation](#inst)
* [Usage](#usage)
 * [readhead()](#readhead)
 * [readsnap()](#readsnap) 
 * [readrockstar()](#readrockstar) 

## <a id="sum"></a>Summary
Do you *love* running simulations but *hate* fighting with the different flavors of output?  If so, you've come to the right place!  `pyGadgetReader` is designed to take the headache out of reading your `GADGET` simulation data; it plays the role of interpreter between the binary snapshot & `python`.  The module currently supports the following data types:

- **Gadget** type1-binary (single/multi-part)
- **HDF5** outputs (single/multi-part)
- **TIPSY** binaries (bin/aux)

`pyGadgetReader` attempts to detect which file format is is dealing with so you don't have to!  It does however, assume that the different formats have unique file extensions (*none*,`.hdf5`,& `.bin`).  

`pyGadgetReader` also supports the following non-GADGET files:

- **Rockstar** binary outputs

## <a id="req"></a>Requirements
* python 2.7.x
* numpy
* h5py

## <a id="obt"></a>Obtaining
The easiest way to download the code and stay up to date is to clone a version from [bitbucket](https://bitbucket.org/rthompson/pygadgetreader) to your local computer via Mercurial (hg):

~~~Bash
> hg clone https://bitbucket.org/rthompson/pygadgetreader
~~~


## <a id="cust"></a>Customization
Before building the module, there are a few customizations you may want to tinker with.  These can all be found in *`readgadget/modules/common.py`*.

1. **UNITS:**  The code currently assumes that your `Gadget` length units are `Kpc/h`, mass units are 10^(10)`Msun/h`, and velocity units are `km/s`.  This can be changed by modifying the *`UnitLength_in_cm`*, *`UnitMass_in_g`*, and *`UnitVelocity_in_cm_per_s`* respectively.  I plan to make this a run-time option in the near future. (*note: has NO impact on TIPSY files as of now*)

2. **METALFACTOR:**  If your metal field contains multiple species `(flag_metals > 1)`, this factor is multiplied to their sum to determine the overall particle's `metallicity`.

3. **BLOCKORDERING:** When dealing with `Gadget` type-1 binaries, the block-ordering can be a pain in the neck.  Custom fields can be added the snapshot causing all of your previous readers to start returning bad data; not anymore!  I've tried to design a system where the user can easily customize their block ordering.  `pyGadgetReader` currently has 2 defaults - `BLOCKORDERING0` and `BLOCKORDERING1` which pertain to two different groups who have different block orderings.  The default is set via the `DEFAULT_BLOCKORDERING` variable, which is used in conjunction with the `BLOCKORDERING` dictionary defined farther down in the file.

    If you require a custom block ordering, use one of the already present block orderings as a template.  The first step is creating a new `OrderedDict` named something along the lines of `BLOCKORDERING3`.  Each entry in the `OrderedDict` represents a data block, and must contain both a KEY and a VALUE.  Once your new block ordering is defined, you should add it to the `BLOCKORDERING` dictionary with an appropriate name.  You can also set it as default via the `DEFAULT_BLOCKORDERING` variable if you so wish.

    **KEY:** This is the name of the block - it MUST be present as a value in the `dataTypes` dictionary defined farther down in `common.py`.

    **VALUE:** This is a list that contains 2 important pieces of information: 1) what particles this block is present for, and 2) should the code first check for a specific header flag?

	Below are two examples.  The first represents the 'pos' block (which is present in the `dataTypes` dictionary).  The second entry is a list telling the code what particle types have this block, where -1 meaning ALL particle types.  The second represents the 'metallicity' data block; here we have a list of [0,4] telling the code that this block is present for particle types 0 (gas) & 4 (stars), and to check the `flag_metals` flag before attempting a read.  You are more than welcome to omit the flag checker if you know for certain the data block exists.

~~~python
...
('pos',[-1]),
('metallicity',[[0,4],'flag_metals']),
...
~~~

## <a id="inst"></a>Installation
Once the code is downloaded there are two methods of installation depending on your access rights.  If you have write access to your python distribution, then the preferred method is to execute the following commands:

~~~Bash
	> python setup.py build     ## this builds the module
	> python setup.py install   ## this installs the module, may require sudo
~~~

If you do *not* have write access to your python install, we need to modify the environment variable `PYTHONPATH` to point to the pyGadgetReader directory.  Add these two lines to your `.bashrc/.bash_profile` (or respective shell file):

~~~Bash
PYTHONPATH=/path/to/pyGadgetReader:$PYTHONPATH
export PYTHONPATH
~~~

****
#### NOTE for uninstalling previous versions:
If you had previously installed my `C` version of `pyGadgetReader` you should *remove* it before trying to use the code as there may be some naming conflicts.  First you need to find out *where* python, a point in the general direction is typing `which python` in your terminal, this will return your installation directory.  Next you need to locate your `site-packages` directory which is usually under python's `lib` directory.  Once there you are looking for anything in the form of `readgadget.so`, once this is found remove it.

## <a id="usage"></a>Usage
**IMPORTANT:** When using `pyGadgetReader`, **do NOT include** the snapshot extension or number prefix (for multi-part).  As an example, if your snapshot is named `'snap_N128L16_005.0.hdf5'`, you would only pass `'snap_N128L16_005'` to the functions below.

To gain access to the following functions, place this at the top of your python script:

~~~python
from readgadget import *
~~~


### <a id="readhead"></a>readhead()
This function reads in the header and returns values of interest.  The values it can read in are as follows:

	 time	       - scale factor of the snapshot
	 redshift      - redshift of the snapshot
	 boxsize       - boxsize if present in units of kpc/h
	 O0	       	   - Omega_0 (Omega_dm+Omega_m)
	 Ol	       	   - Omega_Lambda
	 h	       	   - hubble parameter
	 gascount      - total # of gas   particles [type 0]
	 dmcount       - total # of DM    particles [type 1]
	 diskcount     - total # of disk  particles [type 2]
	 bulgecount    - total # of bulge particles [type 3]
	 starcount     - total # of star  particles [type 4]
	 bndrycount    - total # of bndry particles [type 5]
	 f_sfr	       - Star Formation Rate flag	 0=off 1=on
	 f_fb	       - Feedback flag	     		 0=off 1=on
	 f_cooling     - Cooling flag			 	 0=off 1=on 
	 f_age	       - Stellar Age tracking flag	 0=off 1=on
	 f_metals      - Metal tracking flag  		 0=off 1=on


    Definition:   readhead('a','b',hdf5=0,tipsy=0,debug=0)

	      Parameters
	      ----------
	      a : Input file.  
	      	  Must be input as a string and enclosed in ' ' - see examples.
	      b : Value of interest from the above list.  
	      	  Must be input as a string and enclosed in ' ' - see examples.
	      
	      Optional
	      --------
		   hdf5:   hdf5 file?
		   tipsy:  tipsy file?

    Example:
	      z = readhead('snap_001','redshift')  
	      		- reads redshift value and assigns it to the z variable	      
	      h = readhead('snap_005','h') 
				- reads in the HubbleParam from the snapshot header



### <a id="readsnap"></a>readsnap()
This function does the heavy lifting.  It reads data blocks from the snapshot and returns the requested data for a a specified particle type.

	   Supported data blocks are:

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
	   Sigma       - (gas)         Approximate surface density
	   age	       - (stars)       Formation time of stellar particles
 	   z	       - (gas & stars) Metallicty of gas & star particles (returns total Z)
	   tmax        - (gas & stars) Maximum temp
	   nspawn      - (gas & stars) Number of star particles spawned
	   potential   - (all)         Potential of particles
	   zarray      - (gas & stars) NMETALS array [C,O,Si,Fe]

	   Supported particle types (note tipsy only returns gas/dm/star particles):
	   gas	       - Gas
	   dm	       - Dark Matter
	   disk	       - Disk particles
	   bulge       - Bulge particles
	   star/stars  - Star particles
	   bndry       - Boundary particles
	   

    Definition:	readsnap('a','b','c',units=0,hdf5=0,tipsy=0,debug=0,
									 supress_output=0,blockordering='romeel')

		Parameters
		----------
		a: Input file.
		   Must be input as a string and enclosed in ' '
		b: Data block you are interested in (see above list)
		   Must be input as a string and enclosed in ' '
		c: Particle type you are interested in (see above list)
		   Must be input as a string and enclosed in ' ', or an int (0-6)
		
		Optional
		--------
		     units: Can either be 0 for code units or 1 for CGS
			  hdf5: hdf5 file?
		     tipsy: tipsy file?
		     debug: Shows debug information
	supress_output: if set to 1 no output is printed to the command line
	 blockordering: allows for the user to specify which block ordering to use
					(only valid for Gadget type-1 binaries)
			  
    Example:
		DMpos=readsnap('snap_001','pos','dm')
		 grho=readsnap('snap_005','rho','gas',units=1)
		gtemp=readsnap('snap_005','u','gas',units=1)


### <a id="readrockstar"></a>readrockstar()
This function reads rockstar binary data.  

	   Current supported return data types:

	   	halos		- All halo data
	   	particles 	- particle IDs & respective halo membership

		pos
		corevel
		bulkvel
		m
		r
		child_r
		vmax_r
		mgrav
		vrmax
		rvmax
		rs
		klypin_rs
		vrms
		J
		energy
		spin
		alt_m
		Xoff
		Voff
		b_to_a
		c_to_a
		A
		b_to_a2
		c_to_a2
		A2
		bullock_spin
		kin_to_pot
		m_pe_b
		m_pe_d
		num_p
		num_child_particles
		p_start
		desc
		flags
		n_core
		min_pos_err
		min_vel_err
		min_bulkvel_err

    Definition:	readrockstar('a','b')

		Parameters
		----------
		a: Input binary file.
		   Must be input as a string and enclosed in ' '
			do NOT include the number or .bin extensions!!
		b: Data block you are interested in (see above list)
		   Must be input as a string and enclosed in ' '
					  
    Example:
	   halodata = readrockstar('halos_0037','halos')
			      - returns array with ALL of the above halo data
		halo_cm = readrockstar('halos_0037','pos')
			      - returns center of mass positions
   		   PIDs = readrockstar('halos_0037','particles')
        	      - returns an (N,2) array with:
                    [i,0] = particle ID,
				    [i,1] = halo ID

****

Any comments or suggestions feel free to contact me:

    Robert Thompson
    rthompsonj@gmail.com

