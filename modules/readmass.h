#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/*######################### MASSES ########################################*/
readmass()
{  
  float *simdata;
  int ndim = 1;

  int i;
  int n;
  int k;
  int pc = 0;

  /* GADGET */
  if(Tipsy==0){
    for(j=0;j<NumFiles;j++){
      read_header();
      if(j==0){
	npy_intp dims[1]={header.npartTotal[type]};
	array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
      }
      /*
	if(Units==1){
	printf("### RETURNING MASS IN GADGET UNITS, MULTIPLY BY 1.98892e43 TO CONVER TO GRAMS ###\n");
	}
      */
      
      if(header.mass[type]>0 && header.npart[type]>0){
	printf("non-zero header mass detected - using header mass for %s\n",Type);
	for(n=0;n<header.npart[type];n++)
	  {
	    if(Units==0) MDATA(array,pc)=header.mass[type];
	    if(Units==1) MDATA(array,pc)=header.mass[type]*1e10;//1.98892e43;
	    pc++;
	  }
      }
      else{
	printf("reading mass block for %s\n",Type);
	simdata=(float*)malloc(header.npart[type]*sizeof(float));
	
	skippos();
	skipvel();
	skippid();
	
	Skip; //skip before MASS
	//seek past particle groups not interested in
	if(type==1){
	  if(header.mass[0]==0)
	    fseek(infp, header.npart[0]*sizeof(float),SEEK_CUR);
	}
	if(type==2){
	  for(k=0;k<2;k++){
	    if(header.mass[k]==0)
	      fseek(infp,header.npart[k]*sizeof(float),SEEK_CUR);
	  }
	}      
	if(type==3){
	  for(k=0;k<3;k++){
	    if(header.mass[k]==0)
	      fseek(infp,header.npart[k]*sizeof(float),SEEK_CUR);
	  }
	} 
	if(type==4){
	  for(k=0;k<4;k++){
	    if(header.mass[k]==0)
	      fseek(infp,header.npart[k]*sizeof(float),SEEK_CUR);
	}
	}
	if(type==5){
	  for(k=0;k<5;k++){
	    if(header.mass[k]==0)
	      fseek(infp,header.npart[k]*sizeof(float),SEEK_CUR);
	  }
	}
	//if(type==4) fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
	fread(simdata,header.npart[type]*sizeof(float),1,infp);
	Skip; //skip after MASS
	fclose(infp);
	
	//count = count + header.npart[type];
	for(n=0;n<header.npart[type];n++)
	  {
	    if(Units==0) MDATA(array,pc) = simdata[n];
	    if(Units==1) MDATA(array,pc) = simdata[n]*1e10;//1.98892e43;
	    pc++;
	  }
      }
  }
    if(pc!=header.npartTotal[type]){
      PyErr_Format(PyExc_IndexError,"particle count mismatch! pc=%d  npartTotal[%d]=%d",pc,type,header.npartTotal[type]);
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
    else PyErr_Format(PyExc_IndexError,"Must select gas!\n");
    
    npy_intp dims[1]={nselect};
    array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);

    float tmp=0.;

    if(type==0){ //GAS
      for(n=0;n<nselect;n++){
	fread(&tmp, sizeof(float),1, infp);
	
	fseek(infp,sizeof(float)*MAXDIM*2,SEEK_CUR); //pos & vel
	fseek(infp,sizeof(float),SEEK_CUR);          //rho
	fseek(infp,sizeof(float),SEEK_CUR);          //temp
	fseek(infp,sizeof(float),SEEK_CUR);          //hsmooth
	fseek(infp,sizeof(float),SEEK_CUR);          //metals
	fseek(infp,sizeof(float),SEEK_CUR);          //phi
	
	MDATA(array,n) = tmp;
      }
    }
    else if(type==1){ //DM
      fseek(infp,t_header.ngas * sizeof(struct tipsy_gas),SEEK_CUR); //skip gas
      for(n=0;n<nselect;n++){
	fread(&tmp, sizeof(float),1,infp);
	fseek(infp,sizeof(float)*8,SEEK_CUR);

	MDATA(array,n) = tmp;
      }
    }
    else if(type==4){
      fseek(infp,t_header.ngas * sizeof(struct tipsy_gas),SEEK_CUR); //skip gas
      fseek(infp,t_header.ndark * sizeof(struct tipsy_dm),SEEK_CUR); //skip DM
      for(n=0;n<nselect;n++){
	fread(&tmp, sizeof(float),1,infp);
	fseek(infp, sizeof(float)*10,SEEK_CUR);

	MDATA(array,n) = tmp;
      }
    }
    fclose(infp);
  }
}
