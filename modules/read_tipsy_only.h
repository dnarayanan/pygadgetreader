#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

//used for reading in s_age(19)
read_tipsy()
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
    if(type==4) nselect = t_header.nstar;
    else PyErr_Format(PyExc_IndexError,"Must select star!\n");
    
    npy_intp dims[1]={nselect};
    
    fclose(infp);
    strcpy(auxfile, filename);
    
    if(auxfile[strlen(auxfile)-3] == 'b')
      auxfile[strlen(auxfile)-3] = 'a';
    if(auxfile[strlen(auxfile)-2] == 'i')
      auxfile[strlen(auxfile)-2] = 'u';
    if(auxfile[strlen(auxfile)-1] == 'n')
      auxfile[strlen(auxfile)-1] = 'x';
    

    printf("reading s_age from %s\n", auxfile);


    sprintf(infile,"%s",auxfile);
    if(!(auxfp=fopen(infile,"r"))) {
      ERR = 1;
      PyErr_Format(PyExc_TypeError,"can't open file : '%s'",infile);
      //return NULL;
    }
    
    array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    
    float tmp=0.;
    if(type==4){
      fseek(infp,t_header.ngas * sizeof(struct tipsy_gas_aux),SEEK_CUR);
      for(n=0;n<nselect;n++){
	fseek(infp,sizeof(float)*NMETALS,SEEK_CUR);
	fread(&tmp,sizeof(float),1,infp);
	fseek(infp,sizeof(float),SEEK_CUR);
	fseek(infp,sizeof(int),SEEK_CUR);
	MDATA(array,n) = tmp;
      }
    }
    fclose(auxfp);
  }
}


//read from the envira file.  20=Mhalo, 21=wind form age
read_tipsy_envira()
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
    if(values==22)
      printf("reading rvir from %s\n", auxfile);
    if(values==23)
      printf("reading vvir from %s\n", auxfile);
    

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
	
	if(values==20)
	  fread(&tmp,sizeof(float),1,infp); //mhalo
	else
	  fseek(infp,sizeof(float),SEEK_CUR);

	if(values==22)
	  fread(&tmp,sizeof(float),1,infp);
	else
	  fseek(infp,sizeof(float),SEEK_CUR); //Rvir

	if(values==23)
	  fread(&tmp,sizeof(float),1,infp);
	else
	  fseek(infp,sizeof(float),SEEK_CUR); //Vvir

	fseek(infp,sizeof(float),SEEK_CUR); //Vang
	fseek(infp,sizeof(float),SEEK_CUR); //Vmag
	fseek(infp,sizeof(int),SEEK_CUR); //satswitch
	fseek(infp,sizeof(int),SEEK_CUR); //i

	if(values==21)
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


//read from the future file.
read_tipsy_future()
{
  printf("IN READ TIPSY FUT%03d\n",Future);

  if(Tipsy==0)
    PyErr_Format(PyExc_IndexError,"Not a Tipsy file!\n");
  else{
    int ndim = 1;
    int n;
    
    char futurefile[500];

    long nselect = 0;
    read_header();
    fclose(infp);

    if(type==0) nselect = t_header.ngas;
    else PyErr_Format(PyExc_IndexError,"Must select gas!\n");
    
    npy_intp dims[1]={nselect};

    strcpy(futurefile, filename);
    
    if(futurefile[strlen(futurefile)-3] == 'b')
      futurefile[strlen(futurefile)-3] = 'f';
    if(futurefile[strlen(futurefile)-2] == 'i')
      futurefile[strlen(futurefile)-2] = 'u';
    if(futurefile[strlen(futurefile)-1] == 'n')
      futurefile[strlen(futurefile)-1] = 't';

    printf("reading %s%03d \n",futurefile,Future);
    sprintf(infile,"%s%03d",futurefile,Future);
    if(!(futfp=fopen(infile,"r"))) {
      ERR = 1;
      PyErr_Format(PyExc_TypeError,"can't open file : '%s'",infile);
    }

    array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);

    printf("file read, array allocated, starting loop\n");

    float tmp=0.;
    float itmp=0;
    if(type==0){
      for(n=0;n<nselect;n++){
	fseek(infp,sizeof(int),SEEK_CUR); //ID

	if(values==5) //RHO
	  fread(&tmp,sizeof(float),1,infp);
	else
	  fseek(infp,sizeof(float),SEEK_CUR);

	if(values==4) //temp
	  fread(&tmp,sizeof(float),1,infp);
	else
	  fseek(infp,sizeof(float),SEEK_CUR);

	if(values==11) //metals
	  fread(&tmp,sizeof(float),1,infp);
	else
	  fseek(infp,sizeof(float),SEEK_CUR);

	if(values==24)
	  fread(&tmp,sizeof(float),1,infp);
	else
	  fseek(infp,sizeof(float),SEEK_CUR); //starfrac
	
	if(values==25)  //AGESTARFORM
	  fread(&tmp,sizeof(float),1,infp);
	else
	  fseek(infp,sizeof(float),SEEK_CUR);
	
	/*
	if(values==10) //star age
	  fread(&tmp,sizeof(float),1,infp);
	else
	  fseek(infp,sizeof(float),SEEK_CUR);
	*/

	if(values==26)
	  fread(&itmp,sizeof(int),1,infp);
	else
	  fseek(infp,sizeof(int),SEEK_CUR); //relaunch
	
	fseek(infp,sizeof(int),SEEK_CUR); //i

	//assign to array
	if(values==26)
	  NSPAWNDATA(array,n) = itmp;
	else
	  MDATA(array,n) = tmp;
      }
      printf("finished gas loop\n");
    }
    fclose(futfp);
  }
}
