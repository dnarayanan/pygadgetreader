#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/*##################### T-Max ############################*/
void gadget_readtmax()
{
  float *simdata;
  int ndim = 1;

  int i;
  unsigned int n;
  unsigned int pc = 0;

  unsigned int skip1,skip2;
  char* blocklabel = "TMAX";

  for(j=0;j<NumFiles;j++){
    skip_blocks(values);
    if(header.flag_tmax==0)
      PyErr_Format(PyExc_IndexError,"flag_tmax=%d --> tmax not output",header.flag_tmax);
    
    if(Ngas==0 && Nstar==0)
      PyErr_Format(PyExc_IndexError,"Ngas=0 and Nstar=0 - No Tmax to read!");
    
    if(j==0){
      if(type!=0 && type!=4)
	PyErr_Format(PyExc_IndexError,"Tmax can only be read for gas/stars!!");
      
      npy_intp dims[1]={header.npartTotal[type]};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    }
    
    simdata=(float*)malloc(header.npart[type]*sizeof(float));
    
    fread(&skip1,sizeof(int),1,infp);
    if(type==4)
      fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
    fread(simdata,header.npart[type]*sizeof(float),1,infp);
    if(type==0)
      fseek(infp,header.npart[4]*sizeof(float),SEEK_CUR);
    fread(&skip2,sizeof(int),1,infp);
    errorcheck(skip1,skip2,blocklabel);
    fclose(infp);
    
    for(n=0;n<header.npart[type];n++)
      {
	if(Units==1)
	  MDATA(array,pc) = simdata[n] * tconvert;
	else
	  MDATA(array,pc) = simdata[n];
	pc++;
      }
  }
  if(pc!=header.npartTotal[type])
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");
  return;
}
