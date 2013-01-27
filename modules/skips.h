#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

//DEFINE SKIPS
skippos(){ //skip positions
  Skip;
  fseek(infp,Ntotal_local*sizeof(float)*3,SEEK_CUR);
  Skip;
}
skipvel(){ //skip velocities
  Skip;
  fseek(infp,Ntotal_local*sizeof(float)*3,SEEK_CUR);
  Skip;
}
skippid(){ //skip PIDs
  Skip;
  fseek(infp,Ntotal_local*sizeof(int),SEEK_CUR);
  Skip;
}
skipmass(){ //skip MASS
  if(header.mass[0]==0 && header.npart[0]>0 || 
     header.mass[1]==0 && header.npart[1]>0 ||
     header.mass[2]==0 && header.npart[2]>0 ||
     header.mass[3]==0 && header.npart[3]>0 ||
     header.mass[4]==0 && header.npart[4]>0 ||
     header.mass[5]==0 && header.npart[5]>0){
    Skip;
    if(header.mass[0]==0 && header.npart[0]>0) fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
    if(header.mass[1]==0 && header.npart[1]>0) fseek(infp,header.npart[1]*sizeof(float),SEEK_CUR);
    if(header.mass[2]==0 && header.npart[2]>0) fseek(infp,header.npart[2]*sizeof(float),SEEK_CUR);
    if(header.mass[3]==0 && header.npart[3]>0) fseek(infp,header.npart[3]*sizeof(float),SEEK_CUR);
    if(header.mass[4]==0 && header.npart[4]>0) fseek(infp,header.npart[4]*sizeof(float),SEEK_CUR);
    if(header.mass[5]==0 && header.npart[5]>0) fseek(infp,header.npart[5]*sizeof(float),SEEK_CUR);
    Skip;
  }
}
skipu(){ //skip U
  Skip;
  fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
  Skip;
}
skiprho(){ //skip RHO
  Skip;
  fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
  Skip;
}
skipne(){ //skip NE
  Skip;
  fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
  Skip;
}
skipnh(){ //skip NH
  Skip;
  fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
  Skip;
}
skiphsml(){ //skip HSML
  Skip;
  fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
  Skip;
}
skipsfr(){ //skip SFR
  Skip;
  fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
  Skip;
}
skipage(){ //skip AGE if stars exist
  if(Nstar>0 && header.flag_stellarage==1){
    Skip;
    fseek(infp,header.npart[4]*sizeof(float),SEEK_CUR);
    Skip;
  }
}
skipz(){ //skip Metallicity
  if(header.flag_metals==1){
    Skip;    
    if(Nstar > 0){
      fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
      fseek(infp,header.npart[4]*sizeof(float),SEEK_CUR);
    }
    else{
      fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
    }
    Skip;		 
  }
}
skipfh2(){ //skip fH2
  Skip;
  fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
  Skip;
}
