#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

void gadget_posvel(){
  skip_blocks(values);
  
  float *simdata;
  int ndim = 2;
  
  int i;
  int n;
  unsigned int pc = 0;
  double factor = 1.;
  unsigned int skip1, skip2;

  if(values==0)
    char* blocklabel = "POS";
  if(values==1)
    char* blocklabel = "VEL";

  /* GADGET */
  for(j=0;j<NumFiles;j++){
    skip_blocks(values);
    if(j==0){
      npy_intp dims[2]={header.npartTotal[type],3};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    }
    simdata=(float*)malloc(header.npart[type]*sizeof(float)*3);
    
    fread(&skip1,sizeof(int),1,infp);
    //seek past particle groups not interested in
    if(type>0){
      for(i=1;i<=type;i++)
	fseek(infp,header.npart[i-1]*3*sizeof(float),SEEK_CUR);
    }
    fread(simdata,header.npart[type]*sizeof(float)*3,1,infp);
    //skip past additional particle types
    if(type<5){
      for(i=type+1; i<6; i++)
	fseek(infp, header.npart[i]*sizeof(float)*3,SEEK_CUR);
    }
    fread(&skip2,sizeof(int),1,infp);
    errorcheck(skip1,skip2,blocklabel);
    fclose(infp);
    
    //count = count + header.npart[type];
    for(n=0;n<header.npart[type];n++)
      {
	//correct velocities
	if(Units==1 && values==1) factor = sqrt(header.time);
	
	DATA(array,pc,0) = simdata[3*n]   * factor;
	DATA(array,pc,1) = simdata[3*n+1] * factor;
	DATA(array,pc,2) = simdata[3*n+2] * factor;
	pc++;
      }
  }
  if(pc!=header.npartTotal[type])
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");
  return;
}
