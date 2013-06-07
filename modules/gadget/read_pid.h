#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/*######################### PID ########################################*/
void readpid()
{
  unsigned int *simdata;
  int ndim = 1;

  int k;
  int i;
  unsigned int n;
  unsigned int pc = 0;

  unsigned int skip1,skip2;
  char* blocklabel = "PID";

  for(j=0;j<NumFiles;j++){
    skip_blocks(values);
    if(j==0){
      npy_intp dims[1]={header.npartTotal[type]};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_UINT32);
    }
    simdata=(unsigned int*)malloc(header.npart[type]*sizeof(int));
    
    fread(&skip1,sizeof(int),1,infp);
    //seek past particle groups not interested in
    if(type>0){
      for(i=1;i<=type;i++){
	fseek(infp,header.npart[i-1]*sizeof(int),SEEK_CUR);
      }
    }

    fread(simdata,header.npart[type]*sizeof(int),1,infp);

    if(type<5){
      for(i=type+1; i<6; i++)
	fseek(infp, header.npart[i]*sizeof(int),SEEK_CUR);
    }    

    fread(&skip2,sizeof(int),1,infp);
    errorcheck(skip1,skip2,blocklabel);
    fclose(infp);

    for(n=0;n<header.npart[type];n++)
      {
	PIDDATA(array,pc) = simdata[n];
	pc++;
      }
  }
  if(pc!=header.npartTotal[type])
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");

  return;
}
