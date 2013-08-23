#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

void tipsy_auxmetals(){
  int ndim = 2;
  unsigned int n;
  int i;
  
  char auxfile[500];
  unsigned int nselect = 0;
  read_header();
  if(type==0) nselect = t_header.ngas;
  if(type==4) nselect = t_header.nstar;
  else PyErr_Format(PyExc_IndexError,"Must select gas/star!\n");
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
  }
  
  npy_intp dims[2]={nselect,4};
  array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
  
  if(type==4)
    fseek(infp,t_header.ngas * sizeof(struct tipsy_gas_aux),SEEK_CUR);

  float tmp;
  for(n=0;n<nselect;n++){
    if(type==0){
      for(i=0;i<NMETALS;i++){
	fread(&tmp, sizeof(float),1,infp);
	DATA(array,n,i) = tmp;
      }

      fseek(infp,sizeof(float),SEEK_CUR);         //sfr
      fseek(infp,sizeof(float),SEEK_CUR);         //tmax
      fseek(infp,sizeof(float),SEEK_CUR);         //delaytime
      fseek(infp,sizeof(float),SEEK_CUR);         //ne
      fseek(infp,sizeof(float),SEEK_CUR);         //nh
      fseek(infp,sizeof(int),SEEK_CUR);           //nspawn      
    }
    else{
      for(i=0;i<NMETALS;i++){
	fread(&tmp, sizeof(float),1,infp);
	DATA(array,n,i) = tmp;
      }
      fseek(infp,sizeof(float)*2,SEEK_CUR);
      fseek(infp,sizeof(int),SEEK_CUR);
    }
  }
  fclose(auxfp);
  return;
}

// READ VALUES FROM AUX
void tipsy_aux(){
  int ndim = 1;
  unsigned int n;

  char auxfile[500];
  unsigned int nselect = 0;
  read_header();
  if(type==0) nselect = t_header.ngas;
  else if(type==4) nselect = t_header.nstar;
  else PyErr_Format(PyExc_IndexError,"Must select gas/star!  you selected %d\n",type);
  fclose(infp);
  
  strcpy(auxfile, filename);
  
  /*
    printf("last 3 chars=%c  %c   %c \n",
    auxfile[(strlen(auxfile)-3)],
    auxfile[(strlen(auxfile)-2)],
    auxfile[(strlen(auxfile)-1)]);
  */
  
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
  
  //printf("FILENAME=%s, auxfile=%s,  nselect=%d \n",filename,auxfile,nselect);
  
  printf("reading from %s\n",auxfile);
  sprintf(infile,"%s",auxfile);
  if(!(auxfp=fopen(infile,"r"))) {
    ERR = 1;
    PyErr_Format(PyExc_TypeError,"can't open file : '%s'",infile);
  }
  
  npy_intp dims[1]={nselect};
  if(values==17)
    array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_INT32);
  else
    array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
  
  float tmp=0.;
  int itmp = 0;


  if(type==0){
    for(n=0;n<nselect;n++){
      fseek(infp,sizeof(float)*NMETALS,SEEK_CUR); //metals
      
      //SFR
      if(values==9)
	fread(&tmp,sizeof(float),1,infp);
      else
	fseek(infp,sizeof(float),SEEK_CUR);
      
      //TMAX
      if(values==15)
	fread(&tmp,sizeof(float),1,infp);
      else
	fseek(infp,sizeof(float),SEEK_CUR);
      
      //DELAYTIME
      if(values==16)
	fread(&tmp,sizeof(float),1,infp);
      else
	fseek(infp,sizeof(float),SEEK_CUR);
      
      //NE
      if(values==6)
	fread(&tmp,sizeof(float),1,infp);
      else
	fseek(infp,sizeof(float),SEEK_CUR);
      
      //NH
      if(values==7)
	fread(&tmp,sizeof(float),1,infp);
      else
	fseek(infp,sizeof(float),SEEK_CUR);
      
      //NSPAWN
      if(values==17)
	fread(&itmp,sizeof(int),1,infp);
      else
	fseek(infp,sizeof(int),SEEK_CUR);
      
      if(values==17)
	NSPAWNDATA(array,n) = itmp;
      else
	MDATA(array,n) = tmp;
    }
  }

  if(type==4){
    fseek(infp,t_header.ngas * sizeof(struct tipsy_gas_aux),SEEK_CUR);    
    for(n=0;n<nselect;n++){
      fseek(infp,sizeof(float)*NMETALS,SEEK_CUR);

      //AGE
      if(values==10)
	fread(&tmp,sizeof(float),1,infp);
      else
	fseek(infp,sizeof(float),SEEK_CUR);

      //TMAX
      if(values==15)
	fread(&tmp,sizeof(float),1,infp);
      else
	fseek(infp,sizeof(float),SEEK_CUR);
      
      //NSPAWN
      if(values==17)
	fread(&itmp,sizeof(int),1,infp);
      else
	fseek(infp,sizeof(int),SEEK_CUR);

      if(values==17)
	NSPAWNDATA(array,n) = itmp;
      else
	MDATA(array,n) = tmp;
    }
  }

  fclose(auxfp);
  return;
}
