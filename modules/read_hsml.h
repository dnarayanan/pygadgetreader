#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/*######################### HSML ########################################*/
readHSML()
{  
  float *simdata;
  int ndim = 1;

  int i;
  int n;
  int pc = 0;
  
  /* GADGET */
  if(Tipsy==0){
    for(j=0;j<NumFiles;j++){
      read_header();
      if(Ngas==0){
	PyErr_Format(PyExc_IndexError,"Ngas=0 - No HSML to read!");
	//return NULL;
      }
      if(j==0){
	if(type!=0){
	  PyErr_Format(PyExc_IndexError,"HSML can only be read for gas!!");
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
      
      Skip; //skip before HSML
      fread(simdata,header.npart[type]*sizeof(float),1,infp);
      Skip; //skip after HSML
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
      fseek(infp,sizeof(float)*9,SEEK_CUR);
      fread(&tmp,sizeof(float),1,infp);
      fseek(infp,sizeof(float)*2,SEEK_CUR);
      
      MDATA(array,n) = tmp;
    }
    fclose(infp);
  }
}
