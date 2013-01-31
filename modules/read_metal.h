#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/*######################### Z ########################################*/
readZ()
{  
  float *simdata;
  int ndim = 1;

  int k;
  int i;
  int n;
  int pc = 0;

  /* GADGET */
  if(Tipsy==0){
    for(j=0;j<NumFiles;j++){
      read_header();
      if(header.flag_metals==0){
	PyErr_Format(PyExc_IndexError,"flag_metals=%d --> METALS NOT TRACKED",header.flag_metals);
	//return NULL;
      }
      if(Ngas==0 && Nstar==0 || header.flag_metals==0){
	if(Ngas==0 && Nstar==0) PyErr_Format(PyExc_IndexError,"Nstar=0 and Ngas=0 - No metallicity to read!");
	if(header.flag_metals==0) PyErr_Format(PyExc_IndexError,"flag_metals=%d - No metallicity to read!",header.flag_metals);
	//return NULL;
      }
      if(j==0){
	if(type!=0 && type!=4){
	  PyErr_Format(PyExc_IndexError,"Z can only be read for gas or stars!!");
	  //return NULL;
	}
	npy_intp dims[1]={header.npartTotal[type]};
	array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
      }
      
      //simdata=(float*)malloc(header.npart[type]*sizeof(float));
      
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
      skipdelaytime();
      skipfh2();
      skipsigma();
      skipage();
      
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

      Skip; //skip before Z
      if(type==4)
	fseek(infp,header.npart[0]*sizeof(float)*header.flag_metals,SEEK_CUR);

      float tmp;
      float Z_tmp;
      for(n=0;n<header.npart[type];n++)
	{
	  tmp = 0.;
	  Z_tmp = 0.;
	  for(k=0; k < header.flag_metals; k++)
	    {
	      fread(&tmp,sizeof(float),1,infp);
	      Z_tmp += tmp;
	    }
	  Z_tmp *= (0.0189/0.0147);

	  MDATA(array,pc) = Z_tmp;
	  pc++;
	}
      Skip;
      fclose(infp);
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
      fseek(infp,sizeof(float)*10,SEEK_CUR);
      fread(&tmp,sizeof(float),1,infp);
      fseek(infp,sizeof(float),SEEK_CUR);

      MDATA(array,n) = tmp;
    }
    fclose(infp);
  }
}


/* TIPSY */
/*######################### METAL ARRAY ########################################*/
readmetals()
{
  if(Tipsy==0){
    int ndim = 2;
    int n;
    int i;
    int pc = 0;
    //#################
    for(j=0;j<NumFiles;j++){
      read_header();
      if(header.flag_metals==0){
	PyErr_Format(PyExc_IndexError,"flag_metals=%d --> METALS NOT TRACKED",header.flag_metals);
	//return NULL;
      }
      if(Ngas==0 && Nstar==0 || header.flag_metals==0){
	if(Ngas==0 && Nstar==0) PyErr_Format(PyExc_IndexError,"Nstar=0 and Ngas=0 - No metallicity to read!");
	if(header.flag_metals==0) PyErr_Format(PyExc_IndexError,"flag_metals=%d - No metallicity to read!",header.flag_metals);
	//return NULL;
      }
      if(j==0){
	if(type!=0 && type!=4){
	  PyErr_Format(PyExc_IndexError,"Z can only be read for gas or stars!!");
	  //return NULL;
	}
	npy_intp dims[2]={header.npartTotal[type],header.flag_metals};
	array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
      }
      
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
      skipdelaytime();
      skipfh2();
      skipsigma();
      skipage();
        
      Skip; //skip before Z
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
    if(pc!=header.npartTotal[type]){
      PyErr_Format(PyExc_IndexError,"particle count mismatch!");
      //return NULL;
    }
  }
  //##############

    //PyErr_Format(PyExc_IndexError,"Not a Tipsy file!\n");
  //}
  else{
    //float *simdata;
    int ndim = 2;
   
    int n;

    char auxfile[500];
    long nselect = 0;
    read_header();
    if(type==0) nselect = t_header.ngas;
    else PyErr_Format(PyExc_IndexError,"Must select gas!\n");
    fclose(infp);
    
    strcpy(auxfile, filename);
    
    if(auxfile[strlen(auxfile)-3] == 'b'){
      //printf("MATCHED A  b\n");
      auxfile[strlen(auxfile)-3] = 'a';
    }
    if(auxfile[strlen(auxfile)-2] == 'i'){
      //printf("MATCHED AN i\n");
      auxfile[strlen(auxfile)-2] = 'u';
    }
    if(auxfile[strlen(auxfile)-1] == 'n'){
      //printf("MATCHED A  n\n");
      auxfile[strlen(auxfile)-1] = 'x';
    }
    
    printf("reading METAL ARRAY from %s\n", auxfile);
    sprintf(infile,"%s",auxfile);
    if(!(auxfp=fopen(infile,"r"))) {
      ERR = 1;
      PyErr_Format(PyExc_TypeError,"can't open file : '%s'",infile);
      //return NULL;
    }
  
    npy_intp dims[2]={nselect,4};
    array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);

    float tmp0,tmp1,tmp2,tmp3=0.;
    for(n=0;n<nselect;n++){
      fread(&tmp0, sizeof(float),1,infp);        //metals
      fread(&tmp1, sizeof(float),1,infp);
      fread(&tmp2, sizeof(float),1,infp);
      fread(&tmp3, sizeof(float),1,infp);

      fseek(infp,sizeof(float),SEEK_CUR);         //sfr
      fseek(infp,sizeof(float),SEEK_CUR);         //tmax
      fseek(infp,sizeof(float),SEEK_CUR);         //delaytime
      fseek(infp,sizeof(float),SEEK_CUR);         //ne
      fseek(infp,sizeof(float),SEEK_CUR);         //nh
      fseek(infp,sizeof(int),SEEK_CUR);           //nspawn

      //printf("%e   %e   %e   %e\n",tmp0,tmp1,tmp2,tmp3);

      DATA(array,n,0) = tmp0; //C
      DATA(array,n,1) = tmp1; //O
      DATA(array,n,2) = tmp2; //Si
      DATA(array,n,3) = tmp3; //Fe
    }
    fclose(auxfp);
  }
}
