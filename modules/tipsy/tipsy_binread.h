#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

// POSITIONS
void tipsy_posvel(){
  int ndim = 2;
  unsigned int n;

  unsigned int nread = 0;
  unsigned int nselect = 0;
  read_header();
  if(type==0) nselect = t_header.ngas;
  else if(type==1) nselect = t_header.ndark;
  else if(type==4) nselect = t_header.nstar;
  else PyErr_Format(PyExc_IndexError,"Must select gas, dm, or stars!\n");

  /*
  if(nth_Particle)
    nread = (float)nselect / (float)nth_Particle;
  else{
    nread = nselect;
    nth_Particle = 1;
  }
  
  if(Debug && nth_Particle && Supress==0)
    printf("particles being read in %d/%d\n",nread, nselect);
    */

  nread = Nth(nth_Particle,nselect);

  npy_intp dims[2]={nread,3};
  array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
  
  unsigned int pc = 0;
  float tmp0,tmp1,tmp2=0.;
  if(type==0){  //GAS
    for(n=0;n<nselect;n++){
      fseek(infp,  sizeof(float),SEEK_CUR);
      
      //POS
      if(values==0){
	fread(&tmp0, sizeof(float),1,infp);
	fread(&tmp1, sizeof(float),1,infp);
	fread(&tmp2, sizeof(float),1,infp);
      }
      else
	fseek(infp,sizeof(float)*3,SEEK_CUR);

      //VEL
      if(values==1){
	fread(&tmp0, sizeof(float),1,infp);
	fread(&tmp1, sizeof(float),1,infp);
	fread(&tmp2, sizeof(float),1,infp);
      }
      else
	fseek(infp,sizeof(float)*3,SEEK_CUR);

      fseek(infp,sizeof(float)*5,SEEK_CUR);

      if(n % nth_Particle == 0){
	DATA(array,pc,0) = tmp0;
	DATA(array,pc,1) = tmp1;
	DATA(array,pc,2) = tmp2;
	pc++;
      }
    }
  }
  else if(type==1){  //DM
    fseek(infp,t_header.ngas * sizeof(struct tipsy_gas),SEEK_CUR); //skip gas
    for(n=0;n<nselect;n++){
      fseek(infp, sizeof(float),SEEK_CUR);
      
      //POS
      if(values==0){
	fread(&tmp0, sizeof(float),1,infp);
	fread(&tmp1, sizeof(float),1,infp);
	fread(&tmp2, sizeof(float),1,infp);
      }
      else
	fseek(infp,sizeof(float)*3,SEEK_CUR);
      
      //VEL
      if(values==1){
	fread(&tmp0, sizeof(float),1,infp);
	fread(&tmp1, sizeof(float),1,infp);
	fread(&tmp2, sizeof(float),1,infp);
      }
      else
	fseek(infp,sizeof(float)*3,SEEK_CUR);

      fseek(infp, sizeof(float)*2,SEEK_CUR); //eps & phi
      
      if(n % nth_Particle == 0){
	DATA(array,pc,0) = tmp0;
	DATA(array,pc,1) = tmp1;
	DATA(array,pc,2) = tmp2;
	pc++;
      }
    }
  }
  else if(type==4){  //STAR
    fseek(infp,t_header.ngas * sizeof(struct tipsy_gas),SEEK_CUR); //skip gas
    fseek(infp,t_header.ndark * sizeof(struct tipsy_dm),SEEK_CUR); //skip DM
    for(n=0;n<nselect;n++){
      fseek(infp, sizeof(float),SEEK_CUR);
      
      //POS
      if(values==0){
	fread(&tmp0, sizeof(float),1,infp);
	fread(&tmp1, sizeof(float),1,infp);
	fread(&tmp2, sizeof(float),1,infp);
      }
      else
	fseek(infp,sizeof(float)*3,SEEK_CUR);

      //VEL
      if(values==1){
	fread(&tmp0, sizeof(float),1,infp);
	fread(&tmp1, sizeof(float),1,infp);
	fread(&tmp2, sizeof(float),1,infp);
      }
      else
	fseek(infp,sizeof(float)*3,SEEK_CUR);

      fseek(infp, sizeof(float)*4,SEEK_CUR); //metals, tform, eps, phi
      
      if(n % nth_Particle == 0){
	DATA(array,pc,0) = tmp0;
	DATA(array,pc,1) = tmp1;
	DATA(array,pc,2) = tmp2;
	pc++;
      }
    }
  }
  fclose(infp); 
  return;
}


/* 
   BIN READING! 
   gas:  mass, rho, temp, hsml, metals, phi
   dm:   mass, eps, phi
   star: mass, metals, tform, eps, phi
*/ 

void tipsy_bin(){
  int ndim = 1;
  unsigned int n;

  unsigned int nread = 0;
  unsigned int nselect = 0;
  read_header();
  if(type==0) nselect = t_header.ngas;
  else if(type==1) nselect = t_header.ndark;
  else if(type==4) nselect = t_header.nstar;
  else PyErr_Format(PyExc_IndexError,"Must select gas/dm/star!\n");
  
  /*
  if(nth_Particle)
    nread = (float)nselect / (float)nth_Particle;
  else{
    nread = nselect;
    nth_Particle = 1;
  }
  
  if(Debug && nth_Particle && Supress==0)
    printf("particles being read in %d/%d\n",nselect, nread);
    */

      nread = Nth(nth_Particle,nselect);

  npy_intp dims[1]={nread};
  array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
  
  float tmp=0.;
  unsigned int pc = 0;
  if(type==0){ //GAS
    for(n=0;n<nselect;n++){
      //MASS
      if(values==3)
	fread(&tmp, sizeof(float),1, infp);
      else
	fseek(infp,sizeof(float),SEEK_CUR);
      
      fseek(infp,sizeof(float)*MAXDIM*2,SEEK_CUR); //pos & vel

      //RHO
      if(values==5)
	fread(&tmp, sizeof(float),1, infp);
      else
	fseek(infp,sizeof(float),SEEK_CUR);

      //TEMP
      if(values==4)
	fread(&tmp, sizeof(float),1, infp);
      else
	fseek(infp,sizeof(float),SEEK_CUR);

      //HSML
      if(values==8)
	fread(&tmp, sizeof(float),1, infp);
      else
	fseek(infp,sizeof(float),SEEK_CUR);

      //TOTAL METALS
      if(values==11)
	fread(&tmp, sizeof(float),1, infp);
      else
	fseek(infp,sizeof(float),SEEK_CUR);

      //PHI
      if(values==18)
	fread(&tmp, sizeof(float),1, infp);
      else
	fseek(infp,sizeof(float),SEEK_CUR);

      if(n % nth_Particle == 0){
	MDATA(array,pc) = tmp;
	pc++;
      }
    }
  }
  else if(type==1){ //DM
    fseek(infp,t_header.ngas * sizeof(struct tipsy_gas),SEEK_CUR); //skip gas
    for(n=0;n<nselect;n++){
      //MASS
      if(values==3)
	fread(&tmp, sizeof(float),1, infp);
      else
	fseek(infp,sizeof(float),SEEK_CUR);

      fseek(infp,sizeof(float)*7,SEEK_CUR);

      //PHI
      if(values==18)
	fread(&tmp, sizeof(float),1, infp);
      else
	fseek(infp,sizeof(float),SEEK_CUR);
            
      if(n % nth_Particle == 0){
	MDATA(array,pc) = tmp;
	pc++;
      }
    }
  }
  else if(type==4){
    fseek(infp,t_header.ngas * sizeof(struct tipsy_gas),SEEK_CUR); //skip gas
    fseek(infp,t_header.ndark * sizeof(struct tipsy_dm),SEEK_CUR); //skip DM
    for(n=0;n<nselect;n++){
      //MASS
      if(values==3)
	fread(&tmp, sizeof(float),1, infp);
      else
	fseek(infp,sizeof(float),SEEK_CUR);

      fseek(infp,sizeof(float)*6,SEEK_CUR);

      //TOTAL METALS
      if(values==11)
	fread(&tmp, sizeof(float),1, infp);
      else
	fseek(infp,sizeof(float),SEEK_CUR);

      //TFORM
      if(values==10)
	fread(&tmp, sizeof(float),1, infp);
      else
	fseek(infp,sizeof(float),SEEK_CUR);

      fseek(infp,sizeof(float),SEEK_CUR);  //eps

      //PHI
      if(values==18)
	fread(&tmp, sizeof(float),1, infp);
      else
	fseek(infp,sizeof(float),SEEK_CUR);
      
      if(n % nth_Particle == 0){
	MDATA(array,pc) = tmp;
	pc++;
      }
    }
  }
  fclose(infp);
  return;
}
