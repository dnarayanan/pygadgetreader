#ifdef ENABLE_HDF5

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>
#include <hdf5.h>

void read_gadget_HDF5(){
  int i;
  int ndim;
  unsigned int n;
  unsigned int pc = 0;
  unsigned int nread;
  unsigned int cnt = 0;
  
  char TYPE[500];
  sprintf(TYPE, "/PartType%d", type);

  nread = Nth(nth_Particle,header.npart[type]);

  // UNIT CONVERSION STUFF
  double convert = 1.;
  double UnitTime_in_s, UnitEnergy_in_cgs;
  double MeanWeight = 0.;
  if(Units==1){
    //u conversion to K
    if(values==4){
      UnitTime_in_s      = UnitLength_in_cm / UnitVelocity_in_cm_per_s;      
      UnitEnergy_in_cgs  = UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2);
    }
    
    //rho conversion to CGS
    if(values==5){
      if(header.OmegaLambda > 0. && header.BoxSize > 0.) 
	convert = pow(1. + header.redshift,3.) * UnitMass_in_g/pow(UnitLength_in_cm,3);
      else
	convert = UnitMass_in_g/pow(UnitLength_in_cm,3);
    }

    if(values==13)
      convert = UnitMass_in_g / pow(UnitLength_in_cm,2);

    //mass
    if(values==3)
      convert = UnitMass_in_g / SOLARMASS;
  }

  // ARRAY DECLARATIONS
  unsigned int *piddata;
  float posvel[header.npart[type]][3];
  float metalarray[header.npart[type]][header.flag_metals];
  float *simdata;
  float *Ne;

  // HDF5 STUFF
  hid_t hdf5_file, hdf5_grp, hdf5_dataset;
  hdf5_file    = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  hdf5_grp     = H5Gopen1(hdf5_file, TYPE);

  // CHECK IF DATASET EXISTS
  int status = 0;
  if(values == 0)
    if(H5Lexists(hdf5_grp,"Coordinates",H5P_DEFAULT)==0)
      status = 1;
  if(values == 1)
    if(H5Lexists(hdf5_grp,"Velocities",H5P_DEFAULT)==0)
      status = 1;
  if(values == 2)
    if(H5Lexists(hdf5_grp,"ParticleIDs",H5P_DEFAULT)==0)
      status = 1;
  //if(H5Lexists(hdf5_grp,"Masses",H5P_DEFAULT)==0)
  //  status = 1;
  if(values == 6)
    if(H5Lexists(hdf5_grp,"ElectronAbundance",H5P_DEFAULT)==0)
      status = 1;
  if(values == 4)
    if(H5Lexists(hdf5_grp,"InternalEnergy",H5P_DEFAULT)==0)
      status = 1;
  if(values == 5)
    if(H5Lexists(hdf5_grp,"Density",H5P_DEFAULT)==0)
      status = 1;
  if(values == 7)
    if(H5Lexists(hdf5_grp,"NeutralHydrogenAbundance",H5P_DEFAULT)==0)
      status = 1;
  if(values == 8)
    if(H5Lexists(hdf5_grp,"SmoothingLength",H5P_DEFAULT)==0)
      status = 1;
  if(values == 9)
    if(H5Lexists(hdf5_grp,"StarFormationRate",H5P_DEFAULT)==0)
      status = 1;
  if(values == 10)
    if(H5Lexists(hdf5_grp,"StellarFormationTime",H5P_DEFAULT)==0)
      status = 1;
  if(values == 12)
    if(H5Lexists(hdf5_grp,"FractionH2",H5P_DEFAULT)==0)
      status = 1;
  if(values == 13)
    if(H5Lexists(hdf5_grp,"Sigma",H5P_DEFAULT)==0)
      status = 1;
  if(values == 15)
    if(H5Lexists(hdf5_grp,"TemperatureMax",H5P_DEFAULT)==0)
      status = 1;
  if(values == 16)
    if(H5Lexists(hdf5_grp,"DelayTime",H5P_DEFAULT)==0)
      status = 1;
  if(values == 18)
    if(H5Lexists(hdf5_grp,"Potential",H5P_DEFAULT)==0)
      status = 1;
  if(values == 11 || values == 14)
    if(H5Lexists(hdf5_grp,"Metallicity",H5P_DEFAULT)==0)
      status = 1;
  if(values == 17)
    if(H5Lexists(hdf5_grp,"NstarsSpawn",H5P_DEFAULT)==0)
      status = 1;
  
  int bypass = 0;

  if(status == 1){
    PyErr_Format(PyExc_IndexError, "BLOCK NOT PRESENT!! %d, returning",values);
    return;
  }


  //pos or vel
  if(values == 0 || values == 1){
    for(int i=0; i<header.npart[type];i++)
      for(int j=0; j<3; j++)
	posvel[i][j] = 0.0;
    
    ndim = 2;
    npy_intp dims[2]={nread_total,3};
    array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    
    if(values==0)      hdf5_dataset = H5Dopen1(hdf5_grp, "Coordinates");
    else if(values==1) hdf5_dataset = H5Dopen1(hdf5_grp, "Velocities");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posvel); 
  }
  //PID
  else if(values == 2){
    ndim = 1;
    npy_intp dims[1]={nread_total};
    array   = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_UINT32);
    piddata = (unsigned int*)malloc(header.npart[type]*sizeof(int));
    
    hdf5_dataset = H5Dopen1(hdf5_grp, "ParticleIDs");
    H5Dread(hdf5_dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, piddata);
  }
  //MASS
  else if(values == 3){
    ndim = 1;
    npy_intp dims[1]={nread_total};
    array   = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    simdata = (float*)malloc(header.npart[type]*sizeof(float));
    
    if(header.mass[type] == 0 && header.npart[type] > 0){
      hdf5_dataset = H5Dopen1(hdf5_grp, "Masses");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, simdata);    
    }
    else if(header.mass[type]>0){
      for(int i=0;i<header.npart[type];i++)
	simdata[i] = header.mass[type];
      bypass = 1;
    }
  }
  //standard float values
  //else if((values > 2 && values < 10) || (values > 11 && values < 17) || values == 18){
  //else if((values > 3 && values < 11) || (values > 11 && values < 14) || values == 18){
  else if(values!=0 && values!=1 && values!=2 && values!=3 && values!=11 && values!=14 && values!=17){
    ndim = 1;
    npy_intp dims[1]={nread_total};
    array   = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    simdata = (float*)malloc(header.npart[type]*sizeof(float));
    
    if(values == 4){  
      if(Units==1){  // T!
	Ne = (float*)malloc(header.npart[type]*sizeof(float));
	hdf5_dataset = H5Dopen1(hdf5_grp, "ElectronAbundance");
	H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Ne);
	H5Dclose(hdf5_dataset);
      }

      hdf5_dataset = H5Dopen1(hdf5_grp, "InternalEnergy");
    }
    if(values == 5) hdf5_dataset = H5Dopen1(hdf5_grp, "Density");
    if(values == 6) hdf5_dataset = H5Dopen1(hdf5_grp, "ElectronAbundance");
    if(values == 7) hdf5_dataset = H5Dopen1(hdf5_grp, "NeutralHydrogenAbundance");
    if(values == 8) hdf5_dataset = H5Dopen1(hdf5_grp, "SmoothingLength");
    if(values == 9) hdf5_dataset = H5Dopen1(hdf5_grp, "StarFormationRate");

    if(values == 10) hdf5_dataset = H5Dopen1(hdf5_grp, "StellarFormationTime");

    if(values == 12) hdf5_dataset = H5Dopen1(hdf5_grp, "FractionH2");
    if(values == 13) hdf5_dataset = H5Dopen1(hdf5_grp, "Sigma");

    if(values == 15) hdf5_dataset = H5Dopen1(hdf5_grp, "TemperatureMax");
    if(values == 16) hdf5_dataset = H5Dopen1(hdf5_grp, "DelayTime");
    if(values == 18) hdf5_dataset = H5Dopen1(hdf5_grp, "Potential");

    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, simdata);    
  }
  //METAL
  // reading one float array
  else if((values == 11 || values == 14) && (header.flag_metals == 1)){
    ndim = 1;
    npy_intp dims[1]={nread_total};
    array   = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    simdata = (float*)malloc(header.npart[type]*sizeof(float));
    hdf5_dataset = H5Dopen1(hdf5_grp, "Metallicity");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, simdata);
  }
  // want one float array
  else if(values == 11 && header.flag_metals > 1){
    ndim = 1;
    npy_intp dims[1]={nread_total};
    array   = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    hdf5_dataset = H5Dopen1(hdf5_grp, "Metallicity");
    for(int i=0; i<header.npart[type]; i++)
      for(int j=0; j<header.flag_metals; j++)
	metalarray[i][j] = 0.0;
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, metalarray);
  }
  // metal array
  else if(values == 14 && header.flag_metals > 1){
    ndim = 2;
    npy_intp dims[2]={nread_total,header.flag_metals};
    array   = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    //simdata = (float*)malloc(header.npart[type]*sizeof(float));
    hdf5_dataset = H5Dopen1(hdf5_grp, "Metallicity");
    for(int i=0; i<header.npart[type]; i++)
      for(int j=0; j<header.flag_metals; j++)
	metalarray[i][j] = 0.0;
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, metalarray);
  }
  else if(values == 17){
    ndim = 1;
    npy_intp dims[1]={nread_total};
    array   = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_INT32);
    simdata = (int*)malloc(header.npart[type]*sizeof(int));

    hdf5_dataset = H5Dopen1(hdf5_grp, "NstarsSpawn");
    H5Dread(hdf5_dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, simdata);
  }

  if(bypass == 0)
    H5Dclose(hdf5_dataset);
  H5Gclose(hdf5_grp);
  H5Fclose(hdf5_file);

  float Z_tmp;

  for(n=0;n<header.npart[type];n++)
    {
      if(n % nth_Particle == 0){
	if(values == 0 || values == 1){
	  DATA(array,pc,0) = posvel[n][0];
	  DATA(array,pc,1) = posvel[n][1];
	  DATA(array,pc,2) = posvel[n][2];
	}
	else if(values == 2)
	  PIDDATA(array,pc) = piddata[n];
	
	else if((values > 2 && values < 10) || (values > 11 && values < 17) || values == 18){
	  
	  if(values==4 && Units==1){
	    MeanWeight = 4.0/(3.*H_MASSFRAC+1.+4.*H_MASSFRAC*Ne[n]) * PROTONMASS;
	    convert = MeanWeight / BOLTZMANN * (GAMMA - 1.) * UnitEnergy_in_cgs/UnitMass_in_g;
	  }
	  MDATA(array,pc) = simdata[n] * convert;
	}
	
	else if((values==11 || values==14) && (header.flag_metals == 1))
	  MDATA(array,pc) = simdata[n];
	else if(values==11 && header.flag_metals > 1){
	  Z_tmp = 0.0;
	  for(int k=0; k<header.flag_metals; k++)
	    Z_tmp += metalarray[n][k];
	  Z_tmp *= METALFACTOR;
	  MDATA(array,pc) = Z_tmp;
	}
	else if(values==14 && header.flag_metals > 1){
	  for(int k=0; k<header.flag_metals; k++)
	    DATA(array,n,k) = metalarray[n][k];
	}
	
	pc++;
      }
    }

  //printf("CONVERT=%g, units=%d\n",convert,Units);

  if(pc!=nread_total)
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");

  return;
}


#endif
