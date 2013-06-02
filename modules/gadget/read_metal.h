#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/*######################### Z ########################################*/
void gadget_readZ()
{  
  float *simdata;
  int ndim = 1;

  int k;
  int i;
  unsigned int n;
  unsigned int pc = 0;
  
  unsigned int skip1,skip2;
  char* blocklabel = "Z";
  
  /* GADGET */
  for(j=0;j<NumFiles;j++){
    skip_blocks(values);
    if(header.flag_metals==0)
      PyErr_Format(PyExc_IndexError,"flag_metals=%d --> METALS NOT TRACKED",header.flag_metals);
    
    if(Ngas==0 && Nstar==0 || header.flag_metals==0){
      if(Ngas==0 && Nstar==0) PyErr_Format(PyExc_IndexError,"Nstar=0 and Ngas=0 - No metallicity to read!");
      if(header.flag_metals==0) PyErr_Format(PyExc_IndexError,"flag_metals=%d - No metallicity to read!",header.flag_metals);
    }
    if(j==0){
      if(type!=0 && type!=4)
	PyErr_Format(PyExc_IndexError,"Z can only be read for gas or stars!!");
      
      npy_intp dims[1]={header.npartTotal[type]};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    }
    
    //simdata=(float*)malloc(header.npart[type]*sizeof(float));
    
    /* OLD METHOD FOR 1D Z ARRAY
       Skip; //skip before Z
       if(type==0){
       fread(simdata,header.npart[type]*sizeof(float),1,infp);
       }
       else{
       fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
       fread(simdata,header.npart[type]*sizeof(float),1,infp);
       }
       Skip; //skip after Z
       fclose(infp);
       
       //count = count + header.npart[type];
       for(n=0;n<header.npart[type];n++)
       {
       MDATA(array,pc) = simdata[n];
       pc++;
       }
    */      
    
    fread(&skip1,sizeof(int),1,infp);
    if(type==4)
      fseek(infp,header.npart[0]*sizeof(float)*header.flag_metals,SEEK_CUR);
    
    float tmp;
    float Z_tmp;
    for(n=0;n<header.npart[type];n++)
      {
	Z_tmp = 0.;
	for(k=0; k < header.flag_metals; k++)
	  {
	    fread(&tmp,sizeof(float),1,infp);
	    Z_tmp += tmp;
	  }
	Z_tmp *= METALFACTOR;
	
	MDATA(array,pc) = Z_tmp;
	pc++;
      }
    if(type==0)
      fseek(infp,header.npart[4]*sizeof(float)*header.flag_metals,SEEK_CUR);
    fread(&skip2,sizeof(int),1,infp);
    errorcheck(skip1,skip2,blocklabel);
    fclose(infp);
  }
  if(pc!=header.npartTotal[type])
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");
  
  fclose(infp);
  return;
}

/*######################### METAL ARRAY ########################################*/
void gadget_readmetals()
{
  int ndim = 2;
  unsigned int n;
  int i;
  unsigned int pc = 0;

  unsigned int skip1,skip2;
  char* blocklabel = "METALARRAY";

  //#################
  for(j=0;j<NumFiles;j++){
    skip_blocks(values);
    if(header.flag_metals==0)
      PyErr_Format(PyExc_IndexError,"flag_metals=%d --> METALS NOT TRACKED",header.flag_metals);

    if(Ngas==0 && Nstar==0 || header.flag_metals==0){
      if(Ngas==0 && Nstar==0) PyErr_Format(PyExc_IndexError,"Nstar=0 and Ngas=0 - No metallicity to read!");
      if(header.flag_metals==0) PyErr_Format(PyExc_IndexError,"flag_metals=%d - No metallicity to read!",header.flag_metals);
    }
    if(j==0){
      if(type!=0 && type!=4)
	PyErr_Format(PyExc_IndexError,"Z can only be read for gas or stars!!");

      npy_intp dims[2]={header.npartTotal[type],header.flag_metals};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    }

    fread(&skip1,sizeof(int),1,infp);
    if(type==4)
      fseek(infp,header.npart[0]*sizeof(float)*header.flag_metals,SEEK_CUR);
    
    float tmp = 0.;
    for(n=0;n<header.npart[type];n++){
      for(i=0;i<header.flag_metals;i++){
	fread(&tmp,sizeof(float),1,infp);
	DATA(array,n,i) = tmp;
      }
      pc++;
    }
    if(type==0)
      fseek(infp,header.npart[4]*sizeof(float)*header.flag_metals,SEEK_CUR);
    fread(&skip2,sizeof(int),1,infp);
    errorcheck(skip1,skip2,blocklabel);

    /*
      float tmp0,tmp1,tmp2,tmp3=0.;
      for(n=0;n<header.npart[type];n++)
      {
      fread(&tmp0,sizeof(float),1,infp);
      fread(&tmp1,sizeof(float),1,infp);
      fread(&tmp2,sizeof(float),1,infp);
      fread(&tmp3,sizeof(float),1,infp);
      
      DATA(array,n,0) = tmp0; //C
      DATA(array,n,1) = tmp1; //O
      DATA(array,n,2) = tmp2; //Si
      DATA(array,n,3) = tmp3; //Fe	  
      pc++;
      }
    */
    
  } //end j loop
  if(pc!=header.npartTotal[type])
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");

  fclose(infp);
  return;
}
