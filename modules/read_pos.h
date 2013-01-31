#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/*######################### POS ########################################*/
readpos()
{  
  float *simdata;
  int ndim = 2;
  
  int i;
  int n;
  int pc = 0;

  /* GADGET */
  if(Tipsy==0){
    for(j=0;j<NumFiles;j++){
      read_header();
      if(j==0){
	npy_intp dims[2]={header.npartTotal[type],3};
	array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
      }
      simdata=(float*)malloc(header.npart[type]*sizeof(float)*3);
      
      Skip; //skip before POS
      //seek past particle groups not interested in
      for(i=1;i<=type;i++){
	fseek(infp,header.npart[i-1]*3*sizeof(float),SEEK_CUR);
      }
      fread(simdata,header.npart[type]*sizeof(float)*3,1,infp);
      Skip; //skip after POS
      fclose(infp);
      
      //count = count + header.npart[type];
      for(n=0;n<header.npart[type];n++)
	{
	  DATA(array,pc,0) = simdata[3*n];
	  DATA(array,pc,1) = simdata[3*n+1];
	  DATA(array,pc,2) = simdata[3*n+2];
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
    else if(type==1) nselect = t_header.ndark;
    else if(type==4) nselect = t_header.nstar;
    else PyErr_Format(PyExc_IndexError,"Must select gas, dm, or stars!\n");
    
    npy_intp dims[2]={nselect,3};
    array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    simdata=(float*)malloc(nselect*sizeof(float)*3);
    
    float tmp0,tmp1,tmp2=0.;
    if(type==0){  //GAS
      for(n=0;n<nselect;n++){
	fseek(infp,  sizeof(float),SEEK_CUR);
	
	fread(&tmp0, sizeof(float),1,infp);
	fread(&tmp1, sizeof(float),1,infp);
	fread(&tmp2, sizeof(float),1,infp);
	
	fseek(infp, sizeof(float)*8,SEEK_CUR);
	
	DATA(array,n,0) = tmp0;
	DATA(array,n,1) = tmp1;
	DATA(array,n,2) = tmp2;
      }
    }
    else if(type==1){  //DM
      fseek(infp,t_header.ngas * sizeof(struct tipsy_gas),SEEK_CUR); //skip gas
      for(n=0;n<nselect;n++){
	fseek(infp, sizeof(float),SEEK_CUR);
	
	fread(&tmp0, sizeof(float),1,infp);
	fread(&tmp1, sizeof(float),1,infp);
	fread(&tmp2, sizeof(float),1,infp);
	
	fseek(infp, sizeof(float)*3,SEEK_CUR); //vels
	fseek(infp, sizeof(float)*2,SEEK_CUR); //eps & phi

	DATA(array,n,0) = tmp0;
	DATA(array,n,1) = tmp1;
	DATA(array,n,2) = tmp2;
      }
    }
    else if(type==4){  //STAR
      fseek(infp,t_header.ngas * sizeof(struct tipsy_gas),SEEK_CUR); //skip gas
      fseek(infp,t_header.ndark * sizeof(struct tipsy_dm),SEEK_CUR); //skip DM
      for(n=0;n<nselect;n++){
	fseek(infp, sizeof(float),SEEK_CUR);

	fread(&tmp0, sizeof(float),1,infp);
	fread(&tmp1, sizeof(float),1,infp);
	fread(&tmp2, sizeof(float),1,infp);

	fseek(infp, sizeof(float)*3,SEEK_CUR); //vels
	fseek(infp, sizeof(float)*4,SEEK_CUR); //metals, tform, eps, phi

	DATA(array,n,0) = tmp0;
	DATA(array,n,1) = tmp1;
	DATA(array,n,2) = tmp2;
      }
    }
    fclose(infp); 
  }
}
