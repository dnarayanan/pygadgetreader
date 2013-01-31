#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/*######################### AGE ########################################*/
readage()
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
      if(header.flag_stellarage==0){
	PyErr_Format(PyExc_IndexError,"flag_stellarage=%d --> AGE NOT TRACKED",header.flag_stellarage);
	//return NULL;
      }
      if(Nstar==0){
	PyErr_Format(PyExc_IndexError,"Nstar=0 - No stars to read!");
	//return NULL;
      }
      if(j==0){
	
	if(type!=4){
	  PyErr_Format(PyExc_IndexError,"Age can only be read for stars!!");
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
      
      Skip; //skip before AGE
      fread(simdata,header.npart[type]*sizeof(float),1,infp);
      Skip; //skip after AGE
      fclose(infp);
      
      //NEW
      /*
      double MPC = 3.0857e+24;
      double KM  = 1.e+5;
      H0 = sqrt(8*PI/3);
      t0 = 2./(3.*H0);
      unit_Time = H0*MPC / (100.*header.hubbleparam*KM);
      */

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
    long nselect = 0;
    read_header();
    if(type==4) nselect = t_header.nstar;
    else PyErr_Format(PyExc_IndexError,"Must select stars!\n");

    npy_intp dims[1]={nselect};
    array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);

    fseek(infp,t_header.ngas * sizeof(struct tipsy_gas),SEEK_CUR);
    float tmp=0.;
    for(n=0;n<nselect;n++){
      fseek(infp,sizeof(float)*8,SEEK_CUR);
      fread(&tmp,sizeof(float),1,infp);
      fseek(infp,sizeof(float)*2,SEEK_CUR);

      MDATA(array,n) = tmp;
    }
    fclose(infp);
  }
}
