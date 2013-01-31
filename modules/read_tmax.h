#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/*##################### T-Max ############################*/
readtmax()
{
  float *simdata;
  int ndim = 1;

  int i;
  int n;
  int pc = 0;

  /* GADGET */
  if(Tipsy==0){ 
    for(j=0;j<NumFiles;j++){
      read_header();
      if(header.flag_tmax==0){
	PyErr_Format(PyExc_IndexError,"flag_tmax=%d --> tmax not output",header.flag_tmax);
      }
      if(Ngas==0 && Nstar==0){
	PyErr_Format(PyExc_IndexError,"Ngas=0 and Nstar=0 - No Tmax to read!");
	//return NULL;
      }
      if(j==0){
	if(type!=0 && type!=4){
	  PyErr_Format(PyExc_IndexError,"Tmax can only be read for gas/stars!!");
	  //return NULL;
	}
	npy_intp dims[1]={header.npartTotal[type]};
	array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
      }
      
      simdata=(float*)malloc(header.npart[type]*sizeof(float));
      
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
      skipz();

      Skip; //skip before Tmax
      if(type==4)
	fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
      fread(simdata,header.npart[type]*sizeof(float),1,infp);
      Skip; //skip after Tmax
      fclose(infp);
      
      //count = count + header.npart[type];
      for(n=0;n<header.npart[type];n++)
	{
	  if(Units==1)
	    MDATA(array,pc) = simdata[n] * tconvert;
	  else
	    MDATA(array,pc) = simdata[n];
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
    char auxfile[500];
    long nselect = 0;
    read_header();
    if(type==0) nselect = t_header.ngas;
    else if(type==4) nselect = t_header.nstar;
    else PyErr_Format(PyExc_IndexError,"Must select gas,star!\n");
    fclose(infp);
    strcpy(auxfile, filename);
    
    if(auxfile[strlen(auxfile)-3] == 'b')
      auxfile[strlen(auxfile)-3] = 'a';
    if(auxfile[strlen(auxfile)-2] == 'i')
      auxfile[strlen(auxfile)-2] = 'u';
    if(auxfile[strlen(auxfile)-1] == 'n')
      auxfile[strlen(auxfile)-1] = 'x';
    
    printf("reading tmax from %s\n", auxfile);
      
    sprintf(infile,"%s",auxfile);
    if(!(auxfp=fopen(infile,"r"))) {
      ERR = 1;
      PyErr_Format(PyExc_TypeError,"can't open file : '%s'",infile);
      //return NULL;
    }
      
    npy_intp dims[1]={nselect};
    array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    
    float tmp=0.;
    if(type==4)
      fseek(infp,t_header.ngas * sizeof(struct tipsy_gas_aux),SEEK_CUR);

    for(n=0;n<nselect;n++){
      if(type==0){	
	fseek(infp,sizeof(float)*NMETALS,SEEK_CUR); //metals
	fseek(infp,sizeof(float),SEEK_CUR);         //sfr
	
	fread(&tmp,sizeof(float),1,infp);           //tmax
	
	fseek(infp,sizeof(float),SEEK_CUR);         //delaytime      
	fseek(infp,sizeof(float),SEEK_CUR);         //ne
	fseek(infp,sizeof(float),SEEK_CUR);         //nh
	fseek(infp,sizeof(int),SEEK_CUR);           //nspawn
      }
      else{
	fseek(infp,sizeof(float)*NMETALS,SEEK_CUR); //metals
	fseek(infp,sizeof(float),SEEK_CUR);         //age
	fread(&tmp,sizeof(float),1,infp);           //tmax
	fseek(infp,sizeof(int),SEEK_CUR);               //nspawn
      }
      MDATA(array,n) = tmp;
    }
  }
}
