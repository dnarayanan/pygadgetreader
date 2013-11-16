#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

//used for reading in s_age(19)
void read_tipsy()
{
  if(Tipsy==0)
    PyErr_Format(PyExc_IndexError,"Not a Tipsy file!\n");
  else{
    
    if(type==0 && values==19)
      PyErr_Format(PyExc_IndexError,"MUST SELECT STARS FOR S_AGE\n");

    int ndim = 1;
    int n;
    unsigned int nread = 0;
    
    char auxfile[500];
    unsigned int nselect = 0;
    read_header();
    if(type==4) nselect = t_header.nstar;
    else PyErr_Format(PyExc_IndexError,"Must select star!\n");

    if(nth_Particle)
      nread = (float)nselect / (float)nth_Particle;
    else{
      nread = nselect;
      nth_Particle = 1;
    }
  
    if(Debug && nth_Particle && Supress==0)
      printf("particles being read in %d/%d\n",nselect, nread);
    
    npy_intp dims[1]={nread};
    
    fclose(infp);
    strcpy(auxfile, filename);
    
    if(auxfile[strlen(auxfile)-3] == 'b')
      auxfile[strlen(auxfile)-3] = 'a';
    if(auxfile[strlen(auxfile)-2] == 'i')
      auxfile[strlen(auxfile)-2] = 'u';
    if(auxfile[strlen(auxfile)-1] == 'n')
      auxfile[strlen(auxfile)-1] = 'x';
    
    if(Supress == 0)
      printf("reading s_age from %s\n", auxfile);

    sprintf(infile,"%s",auxfile);
    if(!(auxfp=fopen(infile,"r"))) {
      ERR = 1;
      PyErr_Format(PyExc_TypeError,"can't open file : '%s'",infile);
      //return NULL;
    }
    
    array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    
    float tmp=0.;
    unsigned int pc = 0;
    if(type==4){
      fseek(infp,t_header.ngas * sizeof(struct tipsy_gas_aux),SEEK_CUR);
      for(n=0;n<nselect;n++){
	fseek(infp,sizeof(float)*NMETALS,SEEK_CUR);
	fread(&tmp,sizeof(float),1,infp);
	fseek(infp,sizeof(float),SEEK_CUR);
	fseek(infp,sizeof(int),SEEK_CUR);

	if(n % nth_Particle == 0){
	  MDATA(array,pc) = tmp;
	  pc++;
	}
      }
    }
    fclose(auxfp);
  }
}


//read from the envira file.  20=Mhalo, 21=wind form age
void read_tipsy_envira()
{
  if(Debug)
    printf("IN READ TIPSY ENVIRA\n");

  if(Tipsy==0)
    PyErr_Format(PyExc_IndexError,"Not a Tipsy file!\n");
  else{
    int ndim = 1;
    int n;
    
    char auxfile[500];

    unsigned int nread = 0;
    unsigned int nselect = 0;
    read_header();
    if(type==0) nselect = t_header.ngas;
    //else if(type==1) nselect = t_header.ndark;
    //else if(type==4) nselect = t_header.nstar;
    else PyErr_Format(PyExc_IndexError,"Must select gas,dm,star!\n");
    
    if(nth_Particle)
      nread = (float)nselect / (float)nth_Particle;
    else{
      nread = nselect;
      nth_Particle = 1;
    }
    
    if(Debug && nth_Particle && Supress==0)
      printf("particles being read in %d/%d\n",nselect, nread);

    npy_intp dims[1]={nread};

    fclose(infp);
    strcpy(auxfile, filename);
    
    if(auxfile[strlen(auxfile)-3] == 'b')
      auxfile[strlen(auxfile)-3] = 'e';
    if(auxfile[strlen(auxfile)-2] == 'i')
      auxfile[strlen(auxfile)-2] = 'n';
    if(auxfile[strlen(auxfile)-1] == 'n')
      auxfile[strlen(auxfile)-1] = 'v';

    strcat(auxfile,"ira");

    if(Supress == 0){
      if(values==20)
	printf("reading Mhalo from %s\n", auxfile);
      if(values==21)
	printf("reading wind_age from %s\n", auxfile);
      if(values==22)
	printf("reading rvir from %s\n", auxfile);
      if(values==23)
	printf("reading vvir from %s\n", auxfile);
      if(values==27)
	printf("reading satswitch from %s\n", auxfile);
      if(values==28)
	printf("reading skid_id from %s\n", auxfile);
    }
    

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

    if(values==28)
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_INT32);
    else
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);

    if(Debug)
      printf("file read, array allocated, starting loop\n");

    float tmp=0.;
    int itmp=0;
    unsigned int pc = 0;
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
	if(values==27)
	  fread(&itmp,sizeof(int),1,infp);
	else
	  fseek(infp,sizeof(int),SEEK_CUR); //satswitch
	fseek(infp,sizeof(int),SEEK_CUR); //i

	if(values==21)
	  fread(&tmp,sizeof(float),1,infp); //wind form age
	else
	  fseek(infp,sizeof(float),SEEK_CUR); 

	if(values==28)
	  fread(&itmp,sizeof(int),1,infp);
	else
	  fseek(infp,sizeof(int),SEEK_CUR); //skid_ID

	fseek(infp,sizeof(float),SEEK_CUR); //skid_rad
	
	if(n % nth_Particle == 0){
	  if(values==27 || values==28)
	    NSPAWNDATA(array,pc) = itmp;
	  else
	    MDATA(array,pc) = tmp;
	  pc++;
	}
      }
      if(Debug)
	printf("finished gas loop\n");
    }
    fclose(auxfp);
  }
  return;
}


//read from the future file.
void read_tipsy_future(int Future,int values)
{
  if(Debug)
    printf("IN READ TIPSY FUT%03d\n",Future);

  if(Tipsy==0)
    PyErr_Format(PyExc_IndexError,"Not a Tipsy file!\n");
  else{
    int ndim = 1;
    int n;
    
    char futurefile[500];

    unsigned int nread;
    unsigned int nselect = 0;
    read_header();
    fclose(infp);

    if(type==0) nselect = t_header.ngas;
    else PyErr_Format(PyExc_IndexError,"Must select gas!\n");
    
    if(nth_Particle)
      nread = (float)nselect / (float)nth_Particle;
    else{
      nread = nselect;
      nth_Particle = 1;
    }
    
    if(Debug && nth_Particle && Supress==0)
      printf("particles being read in %d/%d\n", nread, nselect);

    npy_intp dims[1]={nread};

    strcpy(futurefile, filename);
    
    if(futurefile[strlen(futurefile)-3] == 'b')
      futurefile[strlen(futurefile)-3] = 'f';
    if(futurefile[strlen(futurefile)-2] == 'i')
      futurefile[strlen(futurefile)-2] = 'u';
    if(futurefile[strlen(futurefile)-1] == 'n')
      futurefile[strlen(futurefile)-1] = 't';
    
    if(Supress == 0)
      printf("reading %s%03d \n",futurefile,Future);
    sprintf(infile,"%s%03d",futurefile,Future);
    if(!(futfp=fopen(infile,"r"))) {
      ERR = 1;
      PyErr_Format(PyExc_TypeError,"can't open file : '%s'",infile);
    }

    if(values==26) 
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_INT32);
    else
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);

    //printf("file read, array allocated, starting loop\n");

    unsigned int pc = 0;
    float tmp=0.;
    int itmp=0;
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

	if(n % nth_Particle == 0){
	  //assign to array
	  if(values==26)
	    NSPAWNDATA(array,pc) = itmp;
	  else
	    MDATA(array,pc) = tmp;
	  pc++;
	}
      }
      //printf("finished gas loop\n");
    }
    fclose(futfp);
  }
  return;
}


//read from the future file.
void read_tipsy_ID()
{
  if(Debug)
    printf("IN READ TIPSY ID\n");

  if(Tipsy==0)
    PyErr_Format(PyExc_IndexError,"Not a Tipsy file!\n");
  else{
    int ndim = 1;
    int n;
    
    char IDfile[500];

    unsigned int nread;
    unsigned int nselect = 0;
    read_header();
    fclose(infp);

    if(type==0) nselect = t_header.ngas;
    else if(type==1) nselect = t_header.ndark;
    else if(type==4) nselect = t_header.nstar;
    else PyErr_Format(PyExc_IndexError,"Must select gas/dm/star!\n");
    
    if(nth_Particle)
      nread = (float)nselect / (float)nth_Particle;
    else{
      nread = nselect;
      nth_Particle = 1;
    }
    
    if(Debug && nth_Particle && Supress==0)
      printf("particles being read in %d/%d\n",nselect, nread);

    npy_intp dims[1]={nread};

    strcpy(IDfile, filename);
    
    if(IDfile[strlen(IDfile)-3] == 'b')
      IDfile[strlen(IDfile)-3] = 'i';
    if(IDfile[strlen(IDfile)-2] == 'i')
      IDfile[strlen(IDfile)-2] = 'd';
    if(IDfile[strlen(IDfile)-1] == 'n')
      IDfile[strlen(IDfile)-1] = 'n';
    

    if(Supress == 0)
      printf("reading %sum \n",IDfile);
    sprintf(infile,"%sum",IDfile);
    if(!(futfp=fopen(infile,"r"))) {
      ERR = 1;
      PyErr_Format(PyExc_TypeError,"can't open file : '%s'",infile);
    }

    array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_UINT32);

    unsigned int pc = 0;
    unsigned int itmp=0;

    if(type==1)
      fseek(infp,sizeof(int)*t_header.ngas,SEEK_CUR);
    if(type==4)
      fseek(infp,sizeof(int)*(t_header.ngas+t_header.ndark),SEEK_CUR);

    for(n=0;n<nselect;n++){
      fread(&itmp,sizeof(int),1,infp);
      
      if(n % nth_Particle == 0){
	PIDDATA(array,pc) = itmp;
	pc++;
      }
      //printf("finished gas loop\n");
    }
    fclose(futfp);
  }
  return;
}
