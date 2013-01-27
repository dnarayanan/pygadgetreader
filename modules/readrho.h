#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/*######################### RHO ########################################*/
readrho()
{  
  float *simdata;
  int ndim  = 1;

  double convert;
  if(Units==1){
    //values used for converting from code units to cgs
    double solarmass = 1.98892e33;
    double kpctocm   = 3.08568025e21;
    convert   = (pow(10.,10)*solarmass)/pow(kpctocm,3);
  }

  int i;
  int n;
  int pc = 0;

  /* GADGET */
  if(Tipsy==0){
    for(j=0;j<NumFiles;j++){
      read_header();
      if(Ngas==0){
	PyErr_Format(PyExc_IndexError,"Ngas=0 - No density to read!");
	//return NULL;
      }
      if(j==0){
	if(type!=0){
	  PyErr_Format(PyExc_IndexError,"Density can only be read for gas!!");
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
      skipu();
      
      Skip; //skip before RHO
      fread(simdata,header.npart[type]*sizeof(float),1,infp);
      Skip; //skip after RHO
      fclose(infp);
      
      //count = count + header.npart[type];
      for(n=0;n<header.npart[type];n++)
	{
	  if(header.OmegaLambda == 0.){
	    if(Units==0) MDATA(array,pc) = simdata[n];
	    if(Units==1) MDATA(array,pc) = simdata[n]*convert;
	    pc++;
	  }
	  else{
	    if(Units==0) MDATA(array,pc) = simdata[n]*pow(1.+header.redshift,3);
	    if(Units==1) MDATA(array,pc) = simdata[n]*pow(1.+header.redshift,3)*convert;
	    pc++;
	  }
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
      fseek(infp,sizeof(float)*7,SEEK_CUR);
      fread(&tmp,sizeof(float),1,infp);
      fseek(infp,sizeof(float)*4,SEEK_CUR);

      MDATA(array,n) = tmp;
    }
    fclose(infp);
  }
}
