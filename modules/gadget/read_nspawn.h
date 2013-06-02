#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/*##################### N spawn ############################*/
void gadget_readnspawn()
{
  int *simdata;
  int ndim = 1;

  int i;
  unsigned int n;
  unsigned int pc = 0;

  unsigned int skip1,skip2;
  char* blocklabel = "NSPAWN";

  for(j=0;j<NumFiles;j++){
    skip_blocks(values);
    if(header.flag_stellarage==0)
      PyErr_Format(PyExc_IndexError,"flag_stellarage=%d --> Nspawn not output",header.flag_stellarage);
    
    if(Ngas==0 && Nstar==0)
      PyErr_Format(PyExc_IndexError,"Ngas=0 and Nstar=0 - No Nspawn to read!");
    
    if(j==0){
      if(type!=0 && type!=4)
	PyErr_Format(PyExc_IndexError,"Nspawn can only be read for gas/stars!!");
      
      npy_intp dims[1]={header.npartTotal[type]};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_INT32);
    }
    
    simdata=(int*)malloc(header.npart[type]*sizeof(int));
    
    fread(&skip1,sizeof(int),1,infp);
    if(type==4)
      fseek(infp,header.npart[0]*sizeof(int),SEEK_CUR);
    fread(simdata,header.npart[type]*sizeof(int),1,infp);
    if(type==0)
      fseek(infp,header.npart[4]*sizeof(int),SEEK_CUR);
    fread(&skip2,sizeof(int),1,infp);
    errorcheck(skip1,skip2,blocklabel);
    fclose(infp);
    
    for(n=0;n<header.npart[type];n++)
      {
	NSPAWNDATA(array,pc) = simdata[n];
	pc++;
      }
  }
  if(pc!=header.npartTotal[type])
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");
  return;
}
