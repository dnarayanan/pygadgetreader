#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/*######################### MASSES ########################################*/
void gadget_mass()
{  
  float *simdata;
  int ndim = 1;

  int i;
  unsigned int n;
  int k;
  unsigned int pc = 0;
  unsigned int nread;
  unsigned int cnt = 0;

  unsigned int skip1, skip2;
  char* blocklabel = "MASS";

  double convert = UnitMass_in_g / SOLARMASS;

  for(j=0;j<NumFiles;j++){
    skip_blocks(values);
    if(j==0){
      npy_intp dims[1]={nread_total};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    }
    /*
      if(Units==1){
      printf("### RETURNING MASS IN GADGET UNITS, MULTIPLY BY 1.98892e43 TO CONVER TO GRAMS ###\n");
      }
    */

    if(nth_Particle)
      nread = ceil((float)header.npart[type]/(float)nth_Particle);
    else
      nread = header.npart[type];
    
    if(Debug && nth_Particle && Supress==0)
      printf("particles being read in %d/%d\n",nread,header.npart[type]);
    
    if(header.mass[type]>0 && header.npart[type]>0){
      if(Supress==0) printf("non-zero header mass detected - using header mass for %s\n",Type);
      for(n=0;n<nread;n++)
	{
	  if(Units==0) MDATA(array,pc)=header.mass[type];
	  if(Units==1) MDATA(array,pc)=header.mass[type] * convert;
	  pc++;
	}
    }
    else{
      if(Supress==0) printf("reading mass block for %s\n",Type);
      simdata=(float*)malloc(nread*sizeof(float));
      
      fread(&skip1,sizeof(int),1,infp);

      if(nth_Particle){
	cnt = 0;

	if(type==0){
	  for(i=0;i<header.npart[type];i++){
	    if(i % nth_Particle == 0){
	      fread(&simdata[cnt],sizeof(float),1,infp);
	      cnt++;
	    }
	    else
	      fseek(infp,sizeof(float),SEEK_CUR);
	  }
	  for(k=type+1;k<6;k++)
	    if(header.mass[k]==0 && header.npart[k]>0){
	      fseek(infp, header.npart[k]*sizeof(float),SEEK_CUR);
	    }
	}
	else{
	  for(k=0;k<type;k++)
	    if(header.mass[k]==0 && header.npart[k]>0)
	      fseek(infp, header.npart[k]*sizeof(float),SEEK_CUR);

	  for(i=0;i<header.npart[type];i++){
	    if(i % nth_Particle == 0){
	      fread(&simdata[cnt],sizeof(float),1,infp);
	      cnt++;
	    }
	    else
	      fseek(infp,sizeof(float),SEEK_CUR);
	  }

	  for(k=type+1;k<6;k++)
	    if(header.mass[k]==0 && header.npart[k]>0)
	      fseek(infp, header.npart[k]*sizeof(float),SEEK_CUR);
	}

      }
      else{
	if(type==0){
	  fread(simdata,header.npart[type]*sizeof(float),1,infp);
	  for(k=type+1;k<6;k++)
	    if(header.mass[k]==0 && header.npart[k]>0){
	      fseek(infp, header.npart[k]*sizeof(float),SEEK_CUR);
	    }
	}
	else{
	  for(k=0;k<type;k++)
	    if(header.mass[k]==0 && header.npart[k]>0)
	      fseek(infp, header.npart[k]*sizeof(float),SEEK_CUR);
	  fread(simdata,header.npart[type]*sizeof(float),1,infp);
	  for(k=type+1;k<6;k++)
	    if(header.mass[k]==0 && header.npart[k]>0)
	      fseek(infp, header.npart[k]*sizeof(float),SEEK_CUR);
	}
      }
      /*
      //seek past particle groups not interested in
      if(type>0){
	for(k=1;k<=type;k++){
	  if(header.mass[k-1]==0 && header.npart[k-1]>0){
	    if(Debug) printf("skipping past mass block of type type %d before\n",k-1);	    
	    fseek(infp, header.npart[k-1]*sizeof(float),SEEK_CUR);
	  }
	}
      }



      if(type<5){
	for(k=type+1; k<6; k++){
	  if(header.mass[k]==0 && header.npart[k]>0){
	    if(Debug) printf("skipping past mass block of type type %d after\n",k);	    
	    fseek(infp, header.npart[k]*sizeof(float),SEEK_CUR);
	  }
	}
      }
      */

      fread(&skip2,sizeof(int),1,infp);
      errorcheck(skip1,skip2,blocklabel);
      fclose(infp);
      
      //count = count + header.npart[type];
      for(n=0;n<nread;n++)
	{
	  if(Units==0) MDATA(array,pc) = simdata[n];
	  if(Units==1) MDATA(array,pc) = simdata[n] * convert;
	  pc++;
	}
    }
  }
  if(pc!=nread_total){
    PyErr_Format(PyExc_IndexError,"particle count mismatch! pc=%d  npartTotal[%d]=%d",pc,type,header.npartTotal[type]);
    //return NULL;
  }
  return;
}
