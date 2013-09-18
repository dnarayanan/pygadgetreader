#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/*##################### Potential ############################*/
void gadget_readpotential()
{
  float *simdata;
  int ndim = 1;

  int i;
  unsigned int n;
  unsigned int pc = 0;
  unsigned int nread;

  unsigned int skip1,skip2;
  char* blocklabel = "POTENTIAL";

  if(nth_Particle)
    nread = ceil((float)header.npart[type]/(float)nth_Particle);
  else
    nread = header.npart[type];

  if(Debug && nth_Particle && Supress==0)
    printf("particles being read in %d/%d\n",nread,header.npart[type]);

  for(j=0;j<NumFiles;j++){
    skip_blocks(values);
    if(header.flag_potential==0)
      PyErr_Format(PyExc_IndexError,"flag_potential=%d --> potentials not output",header.flag_potential);
    
    if(j==0){
      npy_intp dims[1]={nread};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    }
    
    simdata=(float*)malloc(nread*sizeof(float));
    
    fread(&skip1,sizeof(int),1,infp);
    if(type>0){
      for(i=1;i<type;i++)
	fseek(infp,header.npart[i-1]*sizeof(float),SEEK_CUR);
    }

    if(nth_Particle){
      if(Supress==0)
	printf("READING EVERY %dnth PARTICLE\n",nth_Particle);
      unsigned int cnt = 0;
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

    if(type<5){
      for(i=type+1; i<6; i++)
	fseek(infp,header.npart[i]*sizeof(float),SEEK_CUR);
    }
    fread(&skip2,sizeof(int),1,infp);
    errorcheck(skip1,skip2,blocklabel);
    fclose(infp);
    
    for(n=0;n<nread;n++)
      {
	MDATA(array,pc) = simdata[n];
	pc++;
      }
  }
  if(pc!=nread)
      PyErr_Format(PyExc_IndexError,"particle count mismatch!");
  return;
}


