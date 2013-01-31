#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>
#define Skip fread(&dummy,sizeof(dummy),1,infp)
#define DATA(a,i,j)*((double *) PyArray_GETPTR2(a,i,j))
#define PIDDATA(a,i)*((unsigned int *) PyArray_GETPTR1(a,i))
#define MDATA(a,i)*((double *) PyArray_GETPTR1(a,i))
#define GALIDDATA(a,i,j)*((int *) PyArray_GETPTR2(a,i,j))

#define NSPAWNDATA(a,i)*((int *) PyArray_GETPTR1(a,i))

const char *filename;  
FILE *infp;
FILE *auxfp;
char infile[500];
int NumFiles, Units, j, dummy;
int ERR;
int Tipsy = 0;
int MAXDIM = 3;
int NMETALS = 4;

PyArrayObject *array;

struct io_header
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  int      flag_stellarage;
  int      flag_metals;
  int      npartTotalHighWord[6];
  int      flag_entropy_instead_u;
  int      flag_doubleprecision;

  int flag_lpt_ics;
  int flag_lpt_scalingcator;

  int flag_tmax;
  int flag_potential;
  int flag_delaytime;
  int flag_fH2;

  char     fill[32];  /* fills to 256 Bytes */
} header;

struct tipsy_io_header
{
  double time;
  int ntotal;
  int ndim;
  int ngas;
  int ndark;
  int nstar;
  int pad;
} t_header;

struct tipsy_gas
{
  float mass;
  float pos[3];
  float vel[3];
  float rho;
  float temp;
  float hsmooth;
  float metals;
  float phi;
} t_gas;

struct tipsy_gas_aux
{
  float metals[4];
  float sfr;
  float tmax;
  float delaytime;
  float ne;
  float nh;
  int nspawn;
} t_gas_aux;

struct tipsy_dm
{
  float mass;
  float pos[3];
  float vel[3];
  float eps;
  float phi;
} t_dm;



int Ngas,Ndm,Ndisk,Nbulge,Nstar,Nbdry,Ntotal;
int Ngas_local,Ndm_local,Ndisk_local,Nbulge_local,Nstar_local,Nbdry_local,Ntotal_local;
read_header()
{
  
  /*
  if(NumFiles>1) 
    sprintf(infile,"%s.%d",filename,j);
  else 
    sprintf(infile,"%s",filename);
  */

  
  sprintf(infile,"%s",filename);
  if(!(infp=fopen(infile,"r"))){
    sprintf(infile,"%s.%d",filename,j);
    if(!(infp=fopen(infile,"r"))) {
      ERR = 1;
      PyErr_Format(PyExc_TypeError,"can't open file : '%s'",infile);
      //return NULL;
    }
  }
  

  if(!(infp=fopen(infile,"r"))){
    ERR = 1;
    PyErr_Format(PyExc_TypeError,"can't open file: '%s'",infile);
    //return NULL;
  }

  if(Tipsy==1){
    fread(&t_header, sizeof(t_header),1,infp);
    MAXDIM = t_header.ndim;
    Ngas   = t_header.ngas;
    Ndm    = t_header.ndark;
    Nstar  = t_header.nstar;
    Ntotal = t_header.ntotal;
  }
  else{
    
    // READ HEADER
    Skip;
    fread(&header,sizeof(header),1,infp);
    Skip;
    
    //printf("numfiles=%d \t\t header.num_files=%d \n",NumFiles,header.num_files);
    
    NumFiles = header.num_files;
    
    if(NumFiles != header.num_files){
      PyErr_Format(PyExc_IndexError,"NumFiles(%d) != header.num_files(%d)!",NumFiles,header.num_files);
      //return NULL;
    }
    
    if(j==0){
      // Assign TOTAL particle counts
      Ngas   = header.npartTotal[0];
      Ndm    = header.npartTotal[1];
      Ndisk  = header.npartTotal[2];
      Nbulge = header.npartTotal[3];
      Nstar  = header.npartTotal[4];
      Nbdry  = header.npartTotal[5];
      Ntotal = Ngas+Ndm+Ndisk+Nbulge+Nstar+Nbdry;
    }
    // Assign LOCAL particle counts (if multiple files)
    Ngas_local   = header.npart[0];
    Ndm_local    = header.npart[1];
    Ndisk_local  = header.npart[2];
    Nbulge_local = header.npart[3];
    Nstar_local  = header.npart[4];
    Nbdry_local  = header.npart[5];
    Ntotal_local = Ngas_local+Ndm_local+Ndisk_local+Nbulge_local+Nstar_local+Nbdry_local;
  }
}


char* gas   = "gas";
char* dm    = "dm";
char* disk  = "disk";
char* bulge = "bulge";
char* star  = "star";
char* stars = "stars";
char* bndry = "bndry";
char* Type;
int type;

char* POS   = "pos";
char* VEL   = "vel";
char* PID   = "pid";
char* MASS  = "mass";
char* U     = "u";
char* RHO   = "rho";
char* NE    = "ne";
char* NH    = "nh";
char* HSML  = "hsml";
char* SFR   = "sfr";
char* AGE   = "age";
char* Z     = "z";
char* fH2   = "fH2";
char* fh2   = "fh2";
char* SIGMA = "Sigma";
char* sigma = "sigma";
char* DELAYT= "delaytime";
char* Values;
int values;

/* TIPSY VARS */
char* METALS = "metals";
char* TMAX   = "tmax";
char* NSPAWN = "nspawn";
char* PHI    = "potential";
char* PHI2   = "pot";
char* PHI3   = "phi";
char* S_AGE  = "s_age";
char* MHALO  = "mhalo";
char* W_AGE  = "windage";

assign_type()
{
  // match particle type and assign corresponding integer value
  if(strcmp(Type,gas)==0)        type = 0;
  else if(strcmp(Type,dm)==0)    type = 1;
  else if(strcmp(Type,disk)==0)  type = 2;
  else if(strcmp(Type,bulge)==0) type = 3;
  else if(strcmp(Type,star)==0)  type = 4; //allow for 'star' or 'stars'
  else if(strcmp(Type,stars)==0) type = 4;
  else if(strcmp(Type,bndry)==0) type = 5;
  else{
    PyErr_Format(PyExc_IndexError,"wrong particle type selected");
    //return NULL;
  }
  // match requested data type and assign corresponding integer
  if(strcmp(Values,POS)==0)           values = 0;
  else if(strcmp(Values,VEL)==0)      values = 1;
  else if(strcmp(Values,PID)==0)      values = 2;
  else if(strcmp(Values,MASS)==0)     values = 3;
  else if(strcmp(Values,U)==0)        values = 4;
  else if(strcmp(Values,RHO)==0)      values = 5;
  else if(strcmp(Values,NE)==0)       values = 6;
  else if(strcmp(Values,NH)==0)       values = 7;
  else if(strcmp(Values,HSML)==0)     values = 8;
  else if(strcmp(Values,SFR)==0)      values = 9;
  else if(strcmp(Values,AGE)==0)      values = 10;
  else if(strcmp(Values,Z)==0)        values = 11;
  else if(strcmp(Values,fH2)==0)      values = 12;
  else if(strcmp(Values,fh2)==0)      values = 12;
  else if(strcmp(Values,SIGMA)==0)    values = 13;
  else if(strcmp(Values,sigma)==0)    values = 13;
  else if(strcmp(Values,METALS)==0)   values = 14;
  else if(strcmp(Values,TMAX)==0)     values = 15;
  else if(strcmp(Values,DELAYT)==0)   values = 16;
  else if(strcmp(Values,NSPAWN)==0)   values = 17;
  else if(strcmp(Values,PHI)==0)      values = 18;
  else if(strcmp(Values,PHI2)==0)     values = 18;
  else if(strcmp(Values,PHI3)==0)     values = 18;
  else if(strcmp(Values,S_AGE)==0)    values = 19;
  else if(strcmp(Values,MHALO)==0)    values = 20;
  else if(strcmp(Values,W_AGE)==0)    values = 21;
  else{
    PyErr_Format(PyExc_IndexError,"wrong values type selected");
    //return NULL;
  }

  //  printf("type=%d \t values=%d\n",type,values);

  if(type==0 && Ngas==0){
    PyErr_Format(PyExc_IndexError,"No %s particles!",gas);
    //return NULL;
  }
  if(type==1 && Ndm==0){
    PyErr_Format(PyExc_IndexError,"No %s particles!",dm);
    //return NULL;
  }
  if(type==2 && Ndisk==0){
    PyErr_Format(PyExc_IndexError,"No %s particles!",disk);
    //return NULL;
  }
  if(type==3 && Nbulge==0){
    PyErr_Format(PyExc_IndexError,"No %s particles!",bulge);
    //return NULL;
  }
  if(type==4 && Nstar==0){
    PyErr_Format(PyExc_IndexError,"No %s!",stars);
    //return NULL;
  }
  if(type==5 && Nbdry==0){
    PyErr_Format(PyExc_IndexError,"No %s particles!",bndry);
    //return NULL;
  }

  /*
  if(type==0) printf("%d gas selected, extracting %s data\n",Ngas,Values);  
  if(type==1) printf("%d dm selected, extracting %s data\n",Ndm,Values);
  if(type==2) printf("%d disk selected, extracting %s data\n",Ndisk,Values);
  if(type==3) printf("%d bulge selected, extracting %s data\n",Nbulge,Values);
  if(type==4) printf("%d stars selected, extracting %s data\n",Nstar,Values);
  if(type==5) printf("%d bndry selected, extracting %s data\n",Nbdry,Values);
  */
}
