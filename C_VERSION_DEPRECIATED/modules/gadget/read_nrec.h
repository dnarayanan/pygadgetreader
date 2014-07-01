#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/*##################### Nrec ############################*/
void gadget_readnrec()
{
  float *simdata;
  int ndim = 1;

  int i;
  unsigned int n;
  unsigned int pc = 0;
  unsigned int nread;
  unsigned int cnt = 0;

  unsigned int skip1,skip2;
  char* blocklabel = "NREC";

  //no NREC in Ken's code
#ifdef ALTBLOCK
  PyErr_Format(PyExc_IndexError,"ALTBLOCK DEFINED, no NREC");
  return;
#endif

  for(j=0;j<NumFiles;j++){
    skip_blocks(values);

    if(Ngas==0)
      PyErr_Format(PyExc_IndexError,"Ngas=0 - No Nrec to read!");
    
    if(j==0){
      if(type!=0 && type!=4)
	PyErr_Format(PyExc_IndexError,"Nrec can only be read for gas!!");
      
      npy_intp dims[1]={nread_total};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_SHORT);
    }
    
    /*
    if(nth_Particle)
      nread = ceil((float)header.npart[type]/(float)nth_Particle);
    else
      nread = header.npart[type];
    
    if(Debug && nth_Particle && Supress==0)
      printf("particles being read in %d/%d\n",nread,header.npart[type]);
      */

    nread = Nth(nth_Particle,header.npart[type]);

    simdata=(short int*)malloc(nread*sizeof(short int));
    
    fread(&skip1,sizeof(int),1,infp);
    if(nth_Particle){
      cnt = 0;
      for(i=0;i<header.npart[type];i++){
	if(i % nth_Particle == 0){
	  fread(&simdata[cnt],sizeof(short int),1,infp);
	  cnt++;
	}
	else
	  fseek(infp,sizeof(short int),SEEK_CUR);
      }
    }
    else
      fread(simdata,header.npart[type]*sizeof(short int),1,infp);

    //skip the empty stuff
    fseek(infp,sizeof(short int)*header.npart[type],SEEK_CUR);
    fread(&skip2,sizeof(int),1,infp);    
    errorcheck(skip1,skip2,blocklabel);
    fclose(infp);
    
    for(n=0;n<nread;n++)
      {
	NRECDATA(array,pc) = simdata[n];
	pc++;
      }
  }
  if(pc!=nread_total)
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");
  return;
}
