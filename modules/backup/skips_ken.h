#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

int i = 0;


errorcheck(unsigned int skip1, unsigned int skip2, char *blocklabel){
  if(Debug) printf("checking block %s -- %d  vs  %d\n",blocklabel,skip1,skip2);
  if(skip1 != skip2)
    PyErr_Format(PyExc_IndexError,"skips before and after %s don't match!  %d vs %d",blocklabel,skip1,skip2);  
  else
    return;
}


void skip_blocks(int blockval){

  unsigned int skip1, skip2;
  int k;
  char* blocklabel;

  read_header();

  // POS
  if(blockval==0) return;
  fread(&skip1,sizeof(int),1,infp);
  for(k=0;k<6;k++)
    fseek(infp,3*header.npart[k]*sizeof(float),SEEK_CUR);
  fread(&skip2,sizeof(int),1,infp);
  blocklabel="POS";
  errorcheck(skip1,skip2,blocklabel);

  // VEL
  if(blockval==1) return;
  fread(&skip1,sizeof(int),1,infp);
  for(k=0;k<6;k++)
    fseek(infp,3*header.npart[k]*sizeof(float),SEEK_CUR);
  fread(&skip2,sizeof(int),1,infp);
  blocklabel="VEL";
  errorcheck(skip1,skip2,blocklabel);


  // PID
  if(blockval==2) return;
  fread(&skip1,sizeof(int),1,infp);
  for(k=0;k<6;k++)
    fseek(infp,header.npart[k]*sizeof(int),SEEK_CUR);
  fread(&skip2,sizeof(int),1,infp);
  blocklabel="PID";
  errorcheck(skip1,skip2,blocklabel);

  // MASS
  if(blockval==3) return;
  if(header.mass[0]==0 && header.npart[0]>0 || 
     header.mass[1]==0 && header.npart[1]>0 ||
     header.mass[2]==0 && header.npart[2]>0 ||
     header.mass[3]==0 && header.npart[3]>0 ||
     header.mass[4]==0 && header.npart[4]>0 ||
     header.mass[5]==0 && header.npart[5]>0){
    fread(&skip1,sizeof(int),1,infp);
    if(header.mass[0]==0 && header.npart[0]>0) fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
    if(header.mass[1]==0 && header.npart[1]>0) fseek(infp,header.npart[1]*sizeof(float),SEEK_CUR);
    if(header.mass[2]==0 && header.npart[2]>0) fseek(infp,header.npart[2]*sizeof(float),SEEK_CUR);
    if(header.mass[3]==0 && header.npart[3]>0) fseek(infp,header.npart[3]*sizeof(float),SEEK_CUR);
    if(header.mass[4]==0 && header.npart[4]>0) fseek(infp,header.npart[4]*sizeof(float),SEEK_CUR);
    if(header.mass[5]==0 && header.npart[5]>0) fseek(infp,header.npart[5]*sizeof(float),SEEK_CUR);
    fread(&skip2,sizeof(int),1,infp);
    blocklabel="MASS";
    errorcheck(skip1,skip2,blocklabel);    
  }

  // U
  if(blockval==4) return;
  fread(&skip1,sizeof(int),1,infp);
  fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
  fread(&skip2,sizeof(int),1,infp);
  blocklabel="U";
  errorcheck(skip1,skip2,blocklabel);

  // RHO
  if(blockval==5) return;
  fread(&skip1,sizeof(int),1,infp);
  fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
  fread(&skip2,sizeof(int),1,infp);
  blocklabel="RHO";
  errorcheck(skip1,skip2,blocklabel);

  // NE
  if(blockval==6) return;
  fread(&skip1,sizeof(int),1,infp);
  fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
  fread(&skip2,sizeof(int),1,infp);
  blocklabel="NE";
  errorcheck(skip1,skip2,blocklabel);

  // NH
  if(blockval==7) return;
  fread(&skip1,sizeof(int),1,infp);
  fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
  fread(&skip2,sizeof(int),1,infp);
  blocklabel="NH";
  errorcheck(skip1,skip2,blocklabel);

  // HSML
  if(blockval==8) return;
  fread(&skip1,sizeof(int),1,infp);
  fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
  fread(&skip2,sizeof(int),1,infp);
  blocklabel="HSML";
  errorcheck(skip1,skip2,blocklabel);

  // SFR
  if(blockval==9) return;
  if(header.flag_sfr){
    fread(&skip1,sizeof(int),1,infp);
    fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
    fread(&skip2,sizeof(int),1,infp);
    blocklabel="SFR";
    errorcheck(skip1,skip2,blocklabel);
  }

  /*
  // DELAYTIME
  if(blockval==16) return;
  if(header.flag_delaytime){
    fread(&skip1,sizeof(int),1,infp);
    fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
    fread(&skip2,sizeof(int),1,infp);
    blocklabel="DELAYTIME";
    errorcheck(skip1,skip2,blocklabel);
  }
  */

  // AGE
  if(blockval==10) return;
  if(header.flag_stellarage){
    fread(&skip1,sizeof(int),1,infp);
    fseek(infp,header.npart[4]*sizeof(float),SEEK_CUR);
    fread(&skip2,sizeof(int),1,infp);
    blocklabel="AGE";
    errorcheck(skip1,skip2,blocklabel);
  }

  // Z
  if(blockval==11 || blockval==14) return;
  if(header.flag_metals){
    fread(&skip1,sizeof(int),1,infp);
    fseek(infp,header.npart[0]*sizeof(float)*header.flag_metals,SEEK_CUR);
    if(header.npart[4]>0)
      fseek(infp,header.npart[4]*sizeof(float)*header.flag_metals,SEEK_CUR);
    fread(&skip2,sizeof(int),1,infp);
    blocklabel="Z";
    errorcheck(skip1,skip2,blocklabel);
  }

  // FH2
  if(blockval==12) return;
  if(header.flag_fH2){
    fread(&skip1,sizeof(int),1,infp);
    fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
    fread(&skip2,sizeof(int),1,infp);
    blocklabel="FH2";
    errorcheck(skip1,skip2,blocklabel);
  }

  // SIGMA
  if(blockval==13) return;
  if(header.flag_fH2){
    fread(&skip1,sizeof(int),1,infp);
    fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
    fread(&skip2,sizeof(int),1,infp);
    blocklabel="SIGMA";
    errorcheck(skip1,skip2,blocklabel);
  }


  /*
  // TMAX
  if(blockval==15) return;
  if(header.flag_tmax){
    fread(&skip1,sizeof(int),1,infp);
    fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
    if(header.npart[4]>0)
      fseek(infp,header.npart[4]*sizeof(float),SEEK_CUR);
    fread(&skip2,sizeof(int),1,infp);
    blocklabel="TMAX";
    errorcheck(skip1,skip2,blocklabel);
  }

  // NSPAWN
  if(blockval==17) return;
  fread(&skip1,sizeof(int),1,infp);
  fseek(infp,header.npart[0]*sizeof(int),SEEK_CUR);
  if(header.npart[4]>0)
    fseek(infp,header.npart[4]*sizeof(int),SEEK_CUR);
  fread(&skip2,sizeof(int),1,infp);
  blocklabel="NSPAWN";
  errorcheck(skip1,skip2,blocklabel);
  */

  // POTENTIAL
  if(blockval==18) return;
  if(header.flag_potential){
    fread(&skip1,sizeof(int),1,infp);
    for(k=0;k<6;k++)
      fseek(infp,header.npart[k]*sizeof(float),SEEK_CUR);
    fread(&skip2,sizeof(int),1,infp);
    blocklabel="POTENTIAL";
    errorcheck(skip1,skip2,blocklabel);
  }

  return;
}



/* SEMI OBSOLETE */

//DEFINE SKIPS
/*
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
*/

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

void skiprho(){ //skip RHO
  skipgas();
  return;
}
void skipgas(){
  Skip;
  fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
  Skip;
  return;
}


/*
#define SKIPGAS(y, x) { Skip; fseek(y,x*sizeof(float),SEEK_CUR); Skip; }
then u call SKIPGAS(infp, header.npart[0]);
*/
