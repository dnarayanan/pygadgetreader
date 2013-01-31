#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/*##################### Potential ############################*/
readpotential()
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
      if(header.flag_potential==0){
	PyErr_Format(PyExc_IndexError,"flag_potential=%d --> potentials not output",header.flag_potential);
      }
      if(j==0){
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
      skiptmax();
      skipnspawn();

      Skip; //skip before Potentials
      for(i=1;i<type;i++)
	fseek(infp,header.npart[i-1]*sizeof(float),SEEK_CUR);
      fread(simdata,header.npart[type]*sizeof(float),1,infp);
      Skip; //skip after Potentials
      fclose(infp);
      
      //count = count + header.npart[type];
      for(n=0;n<header.npart[type];n++)
	{
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
    else if(type==1) nselect = t_header.ndark;
    else if(type==4) nselect = t_header.nstar;
    else PyErr_Format(PyExc_IndexError,"Must select gas,dm,star!\n");

    printf("reading potentials from %s\n", auxfile);
      
    npy_intp dims[1]={nselect};
    array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    
    float tmp=0.;
    if(type==1 || type==4)
      fseek(infp,t_header.ngas*sizeof(struct tipsy_gas),SEEK_CUR);
    if(type==4)
      fseek(infp,t_header.ndark*sizeof(struct tipsy_dm),SEEK_CUR);

    for(n=0;n<nselect;n++){
      if(type==0){ //GAS
	fseek(infp,sizeof(float)*11,SEEK_CUR);
	fread(&tmp,sizeof(float),1,infp);
      }
      else if(type==1){ //DM
	fseek(infp,sizeof(float)*8,SEEK_CUR);
	fread(&tmp,sizeof(float),1,infp);
      }
      else if(type==4){ //STAR
	fseek(infp,sizeof(float)*10,SEEK_CUR);
	fread(&tmp,sizeof(float),1,infp);
      }
      MDATA(array,n) = tmp;
    }
  }
}
