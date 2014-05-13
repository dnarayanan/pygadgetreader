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

  //pid array
  unsigned int *piddata;
  
  //pos/vel array
  float posvel[header.npart[type]][3];

  //standard float value
  float *simdata;
  

  // hdf5 stuff
  hid_t hdf5_file, hdf5_grp, hdf5_dataset;  
  hdf5_file    = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  hdf5_grp     = H5Gopen1(hdf5_file, TYPE);


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
    piddata = (unsigned int*)malloc(nread*sizeof(int));
    
    hdf5_dataset = H5Dopen1(hdf5_grp, "ParticleIDs");
    H5Dread(hdf5_dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, piddata);
  }
  //standard float values
  else if((values > 2 && values < 10) || (values > 11 && values < 17)){
    ndim = 1;
    npy_intp dims[1]={nread_total};
    array   = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    simdata = (float*)malloc(nread*sizeof(float));
    
    if(values == 3) hdf5_dataset = H5Dopen1(hdf5_grp, "Masses");
    if(values == 4) hdf5_dataset = H5Dopen1(hdf5_grp, "InternalEnergy");
    if(values == 5) hdf5_dataset = H5Dopen1(hdf5_grp, "Density");
    if(values == 6) hdf5_dataset = H5Dopen1(hdf5_grp, "ElectronAbundance");
    if(values == 7) hdf5_dataset = H5Dopen1(hdf5_grp, "NeutralHydrogenAbundance");
    if(values == 8) hdf5_dataset = H5Dopen1(hdf5_grp, "SmoothingLength");
    if(values == 9) hdf5_dataset = H5Dopen1(hdf5_grp, "StarFormationRate");

    if(values == 12) hdf5_dataset = H5Dopen1(hdf5_grp, "FractionH2");
    if(values == 13) hdf5_dataset = H5Dopen1(hdf5_grp, "Sigma");
    if(values == 15) hdf5_dataset = H5Dopen1(hdf5_grp, "TemperatureMax");
    if(values == 16) hdf5_dataset = H5Dopen1(hdf5_grp, "DelayTime");

    H5Dread(hdf5_dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, simdata);    
  }

  H5Dclose(hdf5_dataset);
  H5Gclose(hdf5_grp);
  H5Fclose(hdf5_file);

  for(n=0;n<nread;n++)
    {

      if(values == 0 || values == 1){
	DATA(array,pc,0) = posvel[n][0];
	DATA(array,pc,1) = posvel[n][1];
	DATA(array,pc,2) = posvel[n][2];
      }
      else if(values == 2)
	PIDDATA(array,pc) = piddata[n];
      else if((values > 2 && values < 10) || (values > 11 && values < 17))
	MDATA(array,pc) = simdata[n];

      pc++;
    }

  if(pc!=nread_total)
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");

  return;
}


#endif
