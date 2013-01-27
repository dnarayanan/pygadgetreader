#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/*##################### fH2 ############################*/
readfh2()
{
  float *simdata;
  int ndim = 1;

  int i;
  int n;
  int pc = 0;
  for(j=0;j<NumFiles;j++){
    read_header();
    if(Ngas==0){
      PyErr_Format(PyExc_IndexError,"Ngas=0 - No fH2 to read!");
      //return NULL;
    }
    if(j==0){
      if(type!=0){
	PyErr_Format(PyExc_IndexError,"fH2 can only be read for gas!!");
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

    Skip; //skip before fH2
    fread(simdata,header.npart[type]*sizeof(float),1,infp);
    Skip; //skip after fH2
    fclose(infp);
    
    //count = count + header.npart[type];
    for(n=0;n<header.npart[type];n++)
      {
	MDATA(array,pc) = simdata[n];
	pc++;
      }
  }
  if(pc!=header.npartTotal[type]){
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");
    //return NULL;
  }
}
