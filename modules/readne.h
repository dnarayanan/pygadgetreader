#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/*######################### NE ########################################*/
readNE()
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
      if(Ngas==0){
	PyErr_Format(PyExc_IndexError,"Ngas=0 - No NE to read!");
	//return NULL;
      }
      if(j==0){
	if(type!=0){
	  PyErr_Format(PyExc_IndexError,"NE can only be read for gas!!");
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
      
      Skip; //skip before NE
      fread(simdata,header.npart[type]*sizeof(float),1,infp);
      Skip; //skip after NE
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
    else PyErr_Format(PyExc_IndexError,"Must select gas!\n");
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

    printf("reading ne from %s\n", auxfile);
    sprintf(infile,"%s",auxfile);
    if(!(auxfp=fopen(infile,"r"))) {
      ERR = 1;
      PyErr_Format(PyExc_TypeError,"can't open file : '%s'",infile);
      //return NULL;
    }
  
    npy_intp dims[1]={nselect};
    array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);

    float tmp=0.;
    for(n=0;n<nselect;n++){
      fseek(infp,sizeof(float)*NMETALS,SEEK_CUR); //metals
      fseek(infp,sizeof(float),SEEK_CUR);         //sfr
      fseek(infp,sizeof(float),SEEK_CUR);         //tmax
      fseek(infp,sizeof(float),SEEK_CUR);         //delaytime

      fread(&tmp,sizeof(float),1,infp);           //ne

      fseek(infp,sizeof(float),SEEK_CUR);         //nh
      fseek(infp,sizeof(int),SEEK_CUR);           //nspawn

      MDATA(array,n) = tmp;
    }
    fclose(auxfp);
  }
}
