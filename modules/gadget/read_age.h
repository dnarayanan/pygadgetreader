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
  //int i;
  unsigned int n;
  unsigned int pc = 0;

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

      npy_intp dims[1]={header.npartTotal[type]};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    }
    
    simdata=(float*)malloc(header.npart[type]*sizeof(float));
    
    fread(&skip1,sizeof(int),1,infp);
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
    
    for(n=0;n<header.npart[type];n++)
      {
	MDATA(array,pc) = simdata[n];
	pc++;
      }
  }
  if(pc!=header.npartTotal[type])
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");

  return;
}
