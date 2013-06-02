#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

void gas_props()
{  
  float *simdata;
  float *Ne;
  int ndim  = 1;

  int i;
  unsigned int n;
  unsigned int pc = 0;
  unsigned int skip1,skip2;

  char* blocklabel = "READ BLOCK";

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
  }

  for(j=0;j<NumFiles;j++){
    skip_blocks(values);
    if(Ngas==0)
      PyErr_Format(PyExc_IndexError,"Ngas=0 - No gas to read!");

    if(j==0){
      if(type!=0)
	PyErr_Format(PyExc_IndexError,"%s can only be read for gas!!",Values);
      
      npy_intp dims[1]={header.npartTotal[type]};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    }
    
    simdata=(float*)malloc(header.npart[type]*sizeof(float)); 
    if(values==4) 
      Ne=(float*)malloc(header.npart[type]*sizeof(float)); 

    // read block
    fread(&skip1,sizeof(int),1,infp);
    fread(simdata,header.npart[type]*sizeof(float),1,infp);
    fread(&skip2,sizeof(int),1,infp);
    errorcheck(skip1,skip2,blocklabel);
    
    // read in Ne for U
    if(values==4 && Units==1){
      skiprho();
      fread(&skip1,sizeof(int),1,infp);
      fread(Ne, header.npart[type]*sizeof(float),1,infp);
      fread(&skip2,sizeof(int),1,infp);
      blocklabel="NEforU";
      errorcheck(skip1,skip2,blocklabel);
    }
    fclose(infp);

    for(n=0;n<header.npart[type];n++)
      {
	if(values==4 && Units==1){
	  MeanWeight = 4.0/(3.*H_MASSFRAC+1.+4.*H_MASSFRAC*Ne[n]) * PROTONMASS;
	  convert = MeanWeight / BOLTZMANN * (GAMMA - 1.) * UnitEnergy_in_cgs/UnitMass_in_g;
	}

	MDATA(array,pc) = simdata[n] * convert;
	pc++;
      }
  }
  if(pc!=header.npartTotal[type]){
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");
  }
  return;
}
