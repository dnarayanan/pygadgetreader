#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

int i = 0;

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
  skipgas();
}
skiprho(){ //skip RHO
  skipgas();
}
skipne(){ //skip NE
  skipgas();
}
skipnh(){ //skip NH
  skipgas();
}
skiphsml(){ //skip HSML
  skipgas();
}
skipsfr(){ //skip SFR
  skipgas();
}
skipdelaytime(){ //skip delay time
  if(header.flag_delaytime)
    skipgas();
}
skipfh2(){ //skip fH2
  if(header.flag_fH2)
    skipgas();
}
skipsigma(){ //skip sigma
  if(header.flag_fH2)
    skipgas();
}
skipage(){ //skip AGE if stars exist
  if(Nstar>0 && header.flag_stellarage==1){
    Skip;
    fseek(infp,header.npart[4]*sizeof(float),SEEK_CUR);
    Skip;
  }
}
skipz(){ //skip Metallicity
  if(header.flag_metals>0){
    Skip;
    fseek(infp,header.npart[0]*sizeof(float)*header.flag_metals,SEEK_CUR);
    if(header.npart[4]>0)
      fseek(infp,header.npart[4]*sizeof(float)*header.flag_metals,SEEK_CUR);
    Skip;
  }   
}
skiptmax(){ //skip max gas temp
  if(header.flag_tmax){
    Skip;
    fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
    fseek(infp,header.npart[4]*sizeof(float),SEEK_CUR);
    Skip;
  }
}
skipnspawn(){
  if(header.flag_stellarage){
    Skip;
    fseek(infp,header.npart[0]*sizeof(int),SEEK_CUR);
    fseek(infp,header.npart[4]*sizeof(int),SEEK_CUR);
    Skip;
  }
}
skippot(){ //skip potentials
  if(header.flag_potential){
    Skip;
    for(i=0; i<6; i++)
      fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
    Skip;
  }
}

/*
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
*/

skipgas(){
  Skip;
  fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
  Skip;
}


/*
#define SKIPGAS(y, x) { Skip; fseek(y,x*sizeof(float),SEEK_CUR); Skip; }
then u call SKIPGAS(infp, header.npart[0]);
*/
