#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/*######################### PID ########################################*/
readpid()
{
  unsigned int *simdata;
  int ndim = 1;

  int i;
  int n;
  unsigned int pc = 0;
  for(j=0;j<NumFiles;j++){
    read_header();
    if(j==0){
      npy_intp dims[1]={header.npartTotal[type]};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_UINT32);
    }
    simdata=(unsigned int*)malloc(header.npart[type]*sizeof(int));
    
    skippos();
    skipvel();

    Skip; //skip before PID
    //seek past particle groups not interested in
    for(i=1;i<=type;i++){
      fseek(infp,header.npart[i-1]*sizeof(int),SEEK_CUR);
    }
    fread(simdata,header.npart[type]*sizeof(int),1,infp);
    Skip; //skip after PID
    fclose(infp);

    //count = count + header.npart[type];
    for(n=0;n<header.npart[type];n++)
      {
	PIDDATA(array,pc) = simdata[n];
	pc++;
      }
  }
  if(pc!=header.npartTotal[type]){
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");
    //return NULL;
  }
}
