#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

//used for reading in phi(18), tmax(15), & nspawn(17)
read_tipsy(int value)
{
  if(Tipsy==0)
    PyErr_Format(PyExc_IndexError,"Not a Tipsy file!\n");
  else{
    
    if(type==0 && values==19)
      PyErr_Format(PyExc_IndexError,"MUST SELECT STARS FOR S_AGE\n");
    
    int ndim = 1;
    int n;
    
    char auxfile[500];
    long nselect = 0;
    read_header();
    if(type==0) nselect = t_header.ngas;
    else if(type==1) nselect = t_header.ndark;
    else if(type==4) nselect = t_header.nstar;
    else PyErr_Format(PyExc_IndexError,"Must select gas,dm,star!\n");
    
    npy_intp dims[1]={nselect};
    
    //phi from bin
    if(values==18){
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
      
      float tmp=0.;
      
      if(type==0){ //GAS
	for(n=0;n<nselect;n++){
	  fseek(infp,sizeof(float)*11,SEEK_CUR);
	  fread(&tmp,sizeof(float),1,infp);
	  
	  MDATA(array,n) = tmp;
	}
      }
      else if(type==1){ //DM
	fseek(infp,t_header.ngas * sizeof(struct tipsy_gas),SEEK_CUR); //skip gas
	for(n=0;n<nselect;n++){
	  fseek(infp,sizeof(float)*8,SEEK_CUR);
	  fread(&tmp,sizeof(float),1,infp);
	  
	  MDATA(array,n) = tmp;
	}
      }
      else if(type==4){ //STAR
	fseek(infp,t_header.ngas * sizeof(struct tipsy_gas),SEEK_CUR); //skip gas
	fseek(infp,t_header.ndark * sizeof(struct tipsy_dm),SEEK_CUR); //skip DM
	for(n=0;n<nselect;n++){
	  fseek(infp,sizeof(float)*10,SEEK_CUR);
	  fread(&tmp,sizeof(float),1,infp);

	  MDATA(array,n) = tmp;
	}
      }
      fclose(infp);
    }
   
    //tmax, delaytime, nspawn from aux
    else{
      fclose(infp);
      strcpy(auxfile, filename);
      
      if(auxfile[strlen(auxfile)-3] == 'b')
	auxfile[strlen(auxfile)-3] = 'a';
      if(auxfile[strlen(auxfile)-2] == 'i')
	auxfile[strlen(auxfile)-2] = 'u';
      if(auxfile[strlen(auxfile)-1] == 'n')
	auxfile[strlen(auxfile)-1] = 'x';
      
      if(values==15)
	printf("reading tmax from %s\n", auxfile);
      if(values==17)
	printf("reading nspawn from %s\n", auxfile);
      
      
      sprintf(infile,"%s",auxfile);
      if(!(auxfp=fopen(infile,"r"))) {
	ERR = 1;
	PyErr_Format(PyExc_TypeError,"can't open file : '%s'",infile);
	//return NULL;
      }
      
      if(values==17)
	array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_UINT32);
      else
	array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
      
      int tmp1=0;
      float tmp=0.;
      if(type==0){
	for(n=0;n<nselect;n++){
	  fseek(infp,sizeof(float)*NMETALS,SEEK_CUR); //metals
	  fseek(infp,sizeof(float),SEEK_CUR);         //sfr
	  
	  if(values==15)
	    fread(&tmp,sizeof(float),1,infp);
	  else
	    fseek(infp,sizeof(float),SEEK_CUR);         //tmax
	  
	  fseek(infp,sizeof(float),SEEK_CUR);         //delaytime
	  
	  fseek(infp,sizeof(float),SEEK_CUR);         //ne
	  fseek(infp,sizeof(float),SEEK_CUR);         //nh
	  
	  if(values==17)
	    fread(&tmp1,sizeof(int),1,infp);
	  else
	    fseek(infp,sizeof(int),SEEK_CUR);           //nspawn
	  
	  if(values==17)
	    NSPAWNDATA(array,n) = tmp1;
	  else
	    MDATA(array,n) = tmp;
	}
      }
      else if(type==4){
	fseek(infp,t_header.ngas * sizeof(struct tipsy_gas),SEEK_CUR);
	fseek(infp,t_header.ndark * sizeof(struct tipsy_dm),SEEK_CUR);
	for(n=0;n<nselect;n++){
	  fseek(infp,sizeof(float)*4,SEEK_CUR);
	  
	  if(values==19)
	    fread(&tmp,sizeof(float),1,infp);
	  else
	    fseek(infp,sizeof(float),SEEK_CUR);

	  if(values==15)
	    fread(&tmp,sizeof(float),1,infp);
	  else
	    fseek(infp,sizeof(float),SEEK_CUR);

	  if(values==17)
	    fread(&tmp1,sizeof(int),1,infp);
	  else
	    fseek(infp,sizeof(int),SEEK_CUR);

	  if(values==17)
	    NSPAWNDATA(array,n) = tmp1;
	  else
	    MDATA(array,n) = tmp;
	}
      }
      fclose(auxfp);
    }
  }
}


//read from the envira file.  20=Mhalo, 21=wind form age
read_tipsy_envira(int value)
{
  printf("IN READ TIPSY ENVIRA\n");

  if(Tipsy==0)
    PyErr_Format(PyExc_IndexError,"Not a Tipsy file!\n");
  else{
    int ndim = 1;
    int n;
    
    char auxfile[500];

    long nselect = 0;
    read_header();
    if(type==0) nselect = t_header.ngas;
    //else if(type==1) nselect = t_header.ndark;
    //else if(type==4) nselect = t_header.nstar;
    else PyErr_Format(PyExc_IndexError,"Must select gas,dm,star!\n");
    
    npy_intp dims[1]={nselect};

    fclose(infp);
    strcpy(auxfile, filename);
    
    if(auxfile[strlen(auxfile)-3] == 'b')
      auxfile[strlen(auxfile)-3] = 'e';
    if(auxfile[strlen(auxfile)-2] == 'i')
      auxfile[strlen(auxfile)-2] = 'n';
    if(auxfile[strlen(auxfile)-1] == 'n')
      auxfile[strlen(auxfile)-1] = 'v';

    strcat(auxfile,"ira");

    if(values==20)
      printf("reading Mhalo from %s\n", auxfile);
    if(values==21)
      printf("reading wind_age from %s\n", auxfile);

    /*
    if(values==17)
      printf("reading nspawn from %s\n", auxfile);
    */

    sprintf(infile,"%s",auxfile);
    if(!(auxfp=fopen(infile,"r"))) {
      ERR = 1;
      PyErr_Format(PyExc_TypeError,"can't open file : '%s'",infile);
      //return NULL;
    }

    array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);

    printf("file read, array allocated, starting loop\n");

    float tmp=0.;
    if(type==0){
      for(n=0;n<nselect;n++){
	fseek(infp,sizeof(int),SEEK_CUR); //ID
	fseek(infp,sizeof(int),SEEK_CUR); //phase
	
	if(value==20)
	  fread(&tmp,sizeof(float),1,infp); //mhalo
	else
	  fseek(infp,sizeof(float),SEEK_CUR);

	fseek(infp,sizeof(float),SEEK_CUR); //Rvir
	fseek(infp,sizeof(float),SEEK_CUR); //Vvir
	fseek(infp,sizeof(float),SEEK_CUR); //Vang
	fseek(infp,sizeof(float),SEEK_CUR); //Vmag
	fseek(infp,sizeof(int),SEEK_CUR); //satswitch
	fseek(infp,sizeof(int),SEEK_CUR); //i

	if(value==21)
	  fread(&tmp,sizeof(float),1,infp); //wind form age
	else
	  fseek(infp,sizeof(float),SEEK_CUR); 

	fseek(infp,sizeof(int),SEEK_CUR); //skid_ID
	fseek(infp,sizeof(float),SEEK_CUR); //skid_rad
	
	MDATA(array,n) = tmp;
      }
      printf("finished gas loop\n");
    }
    fclose(auxfp);
  }
}
