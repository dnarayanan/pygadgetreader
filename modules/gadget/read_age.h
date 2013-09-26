#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/*######################### AGE ########################################*/
void gadget_readage()
{  
  float *simdata;
  int ndim = 1;
  int i;
  unsigned int n;
  unsigned int pc = 0;
  unsigned int nread;
  unsigned int cnt = 0;

  unsigned int skip1,skip2;
  char* blocklabel = "AGE";

  for(j=0;j<NumFiles;j++){
    skip_blocks(values);
    if(header.flag_stellarage==0)
      PyErr_Format(PyExc_IndexError,"flag_stellarage=%d --> AGE NOT TRACKED",header.flag_stellarage);
    
    if(Nstar==0)
      PyErr_Format(PyExc_IndexError,"Nstar=0 - No stars to read!");
    
    if(j==0){
      if(type!=4)
	PyErr_Format(PyExc_IndexError,"Age can only be read for stars!!");

      npy_intp dims[1]={nread_total};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    }

    if(nth_Particle)
      nread = ceil((float)header.npart[type]/(float)nth_Particle);
    else
      nread = header.npart[type];
    
    if(Debug && nth_Particle && Supress==0)
      printf("particles being read in %d/%d\n",nread,header.npart[type]);

    simdata=(float*)malloc(nread*sizeof(float));
    
    fread(&skip1,sizeof(int),1,infp);
    if(nth_Particle){
      cnt = 0;
      for(i=0;i<header.npart[type];i++){
	if(i % nth_Particle == 0){
	  fread(&simdata[cnt],sizeof(float),1,infp);
	  cnt++;
	}
	else
	  fseek(infp,sizeof(float),SEEK_CUR);
      }
    }
    else
      fread(simdata,header.npart[type]*sizeof(float),1,infp);
    fread(&skip2,sizeof(int),1,infp);
    errorcheck(skip1,skip2,blocklabel);
    fclose(infp);
    
    //NEW
    /*
      double MPC = 3.0857e+24;
      double KM  = 1.e+5;
      H0 = sqrt(8*PI/3);
      t0 = 2./(3.*H0);
      unit_Time = H0*MPC / (100.*header.hubbleparam*KM);
    */
    
    for(n=0;n<nread;n++)
      {
	MDATA(array,pc) = simdata[n];
	pc++;
      }
  }
  if(pc!=nread_total)
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");

  return;
}
