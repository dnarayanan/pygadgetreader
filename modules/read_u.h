#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/*######################### INTERNAL ENERGY ########################################*/
readu()
{  
  float *simdata;
  int ndim  = 1;
  
  double convert;
  if(Units==1){
    double boltzmann   = 1.380688e-16;   //erg/kelvin
    double proton_mass = 1.67262158e-24;  //grams
    double kmtocm      = 1.e5;
    double gammaminus1 = (5./3.)-1.;
    convert     = gammaminus1*(proton_mass/boltzmann)*pow(kmtocm,2);
  }
  int i;
  int n;
  int pc = 0; 

  /* GADGET */
  if(Tipsy==0){
    for(j=0;j<NumFiles;j++){
      read_header();
      if(Ngas==0){
	PyErr_Format(PyExc_IndexError,"Ngas=0 - No internal energy to read!");
	//return NULL;
      }
      if(j==0){
	if(type!=0){
	  PyErr_Format(PyExc_IndexError,"Internal Energy can only be read for gas!!");
	  //return NULL;
	}
	npy_intp dims[1]={header.npartTotal[type]};
	array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
      }
      
      simdata=(float*)malloc(header.npart[type]*sizeof(float));
      
      skippos();
      skipvel();
      skippid();
      skipmass();
      
      Skip; //skip before U
      fread(simdata,header.npart[type]*sizeof(float),1,infp);
      Skip; //skip after U
      fclose(infp);
      
      //count = count + header.npart[type];
      for(n=0;n<header.npart[type];n++)
	{
	  if(Units==0) MDATA(array,pc) = simdata[n];
	  if(Units==1) MDATA(array,pc) = simdata[n]*convert;
	  pc++;
	}
    }
    if(pc!=header.npartTotal[type]){
      PyErr_Format(PyExc_IndexError,"particle count mismatch!");
      //return NULL;
    }
  }

  /* TIPSY */
  else{
    long nselect = 0;
    read_header();
    if(type==0) nselect = t_header.ngas;
    else PyErr_Format(PyExc_IndexError,"Must select gas!\n");

    npy_intp dims[1]={nselect};
    array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);

    float tmp=0.;
    for(n=0;n<nselect;n++){
      fseek(infp,sizeof(float)*8,SEEK_CUR);
      fread(&tmp,sizeof(float),1,infp);
      fseek(infp,sizeof(float)*3,SEEK_CUR);

      MDATA(array,n) = tmp;
    }
    fclose(infp);
  }
}
