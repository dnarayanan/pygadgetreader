#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/*##################### sigma ############################*/
readsigma()
{
  float *simdata;
  int ndim = 1;

  double convert;
  if(Units==1){
    //values used for converting from code units to cgs
    double solarmass = 1.98892e33;
    double kpctocm   = 3.08568025e21;
    convert   = (pow(10.,10)*solarmass)/pow(kpctocm,2);
  }

  int i;
  int n;
  int pc = 0;
  for(j=0;j<NumFiles;j++){
    read_header();
    if(Ngas==0){
      PyErr_Format(PyExc_IndexError,"Ngas=0 - No Sigma to read!");
      //return NULL;
    }
    if(j==0){
      if(type!=0){
	PyErr_Format(PyExc_IndexError,"Sigma can only be read for gas!!");
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
    skiprho();
    skipne();
    skipnh();
    skiphsml();
    skipsfr();
    skipage();
    skipz();
    skipfh2();

    Skip; //skip before sigma
    fread(simdata,header.npart[type]*sizeof(float),1,infp);
    Skip; //skip after sigma
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
