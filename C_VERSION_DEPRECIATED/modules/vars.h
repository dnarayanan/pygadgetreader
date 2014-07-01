
//DEPRECIATED!!!!  NO TOUCHING!!!
/* -- enable this if you are in Ken's group --*/
//#define ALTBLOCK
/*--------------------------------------------*/

/* -- enable this to read tipsy by default  --*/
//#define TIPSY_DEFAULT
/*--------------------------------------------*/

//#define ENABLE_HDF5
//#define HDF5_DEFAULT

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>
#define Skip fread(&dummy,sizeof(dummy),1,infp)
#define DATA(a,i,j)*((double *) PyArray_GETPTR2(a,i,j))
#define DMDATA(a,i,j)*((unsigned int *) PyArray_GETPTR2(a,i,j))
#define PIDDATA(a,i)*((unsigned int *) PyArray_GETPTR1(a,i))
#define MDATA(a,i)*((double *) PyArray_GETPTR1(a,i))
#define GALIDDATA(a,i,j)*((int *) PyArray_GETPTR2(a,i,j))

#define NSPAWNDATA(a,i)*((int *) PyArray_GETPTR1(a,i))

#define NRECDATA(a,i)*((short int *) PyArray_GETPTR1(a,i))

#define H_MASSFRAC 0.76
#define BOLTZMANN 1.3806e-16
#define PROTONMASS 1.6726e-24
#define GAMMA 5./3.
#define CM_PER_MPC 3.085678e24
#define UnitLength_in_cm 3.085678e21
#define UnitMass_in_g 1.989e43
#define UnitVelocity_in_cm_per_s 1.e5
#define SOLARMASS 1.989e33

#ifndef ALTBLOCK  //romeel's group
//#define METALFACTOR 0.0189/0.0147
float METALFACTOR = 0.0189/0.0147;
int NMETALS = 4;
#else           //ken's group
//#define METALFACTOR 1.
float METALFACTOR = 1.;
int NMETALS = 1;
#endif

#ifdef TIPSY_DEFAULT
int Tipsy = 1;
#else
int Tipsy = 0;
#endif
int Future = 0;
int MAXDIM = 3;
int i = 0;

const char *filename;  
FILE *infp;
FILE *auxfp;
FILE *futfp;
char infile[500];
int NumFiles, Units, j, dummy;
int HDF5_FILE;
int Debug;
int Supress;
int ERR;
int nth_Particle;
float nth_P;
int nMetals;
unsigned int nread_total;

PyArrayObject *array;

#ifndef ALTBLOCK
struct io_header
{
  unsigned int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  unsigned int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  int      flag_stellarage;
  int      flag_metals;
  unsigned int      npartTotalHighWord[6];
  int      flag_entropy_instead_u;
  int      flag_doubleprecision;

  int flag_potential;
  int flag_fH2;

  int flag_tmax;
  int flag_delaytime;

  int flag_lpt_ics;
  float flag_lpt_scalingcator;

  char     fill[32];  /* fills to 256 Bytes */
} header;
#else //Ken's header
struct io_header
{
  unsigned int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  unsigned int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  int      flag_stellarage;
  int      flag_metals;
  unsigned int      npartTotalHighWord[6];
  int      flag_entropy_instead_u;
  int      flag_doubleprecision;

  int flag_potential;
  int flag_fH2;

  int flag_delaytime;  //DUMMY
  int flag_tmax;       //DUMMY

  char     fill[40];  /* fills to 256 Bytes */
} header;
#endif


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


// ROCKSTAR STUFF
#define BINARY_HEADER_SIZE 256
#define VERSION_MAX_SIZE 12
struct rs_binary_output_header {
  uint64_t magic;
  int64_t snap, chunk;
  float scale, Om, Ol, h0;
  float bounds[6];
  int64_t num_halos, num_particles;
  float box_size, particle_mass;
  int64_t particle_type;
  int32_t format_revision;
  char rockstar_version[VERSION_MAX_SIZE];
  char unused[BINARY_HEADER_SIZE - (sizeof(char)*VERSION_MAX_SIZE) - (sizeof(float)*12) - sizeof(int32_t) - (sizeof(int64_t)*6)];
} rs_header;
/*
struct rs_binary_output_header {
  uint64_t magic;
  int64_t snap, chunk;
  float scale, Om, Ol, h0;
  float bounds[6];
  int64_t num_halos, num_particles;
  float box_size, particle_mass;
  int64_t particle_type;
  char unused[BINARY_HEADER_SIZE - (sizeof(float)*12) - (sizeof(int64_t)*6)];
} rs_header;
*/
struct rshalo {
  int64_t id;
  float pos[6], corevel[3], bulkvel[3];
  float m, r, child_r, vmax_r, mgrav, vmax, rvmax, rs, klypin_rs, vrms,
    J[3], energy, spin, alt_m[4], Xoff, Voff, b_to_a, c_to_a, A[3],
    b_to_a2, c_to_a2, A2[3],
    bullock_spin, kin_to_pot, m_pe_b, m_pe_d;
  int64_t num_p, num_child_particles, p_start, desc, flags, n_core;
  float min_pos_err, min_vel_err, min_bulkvel_err;
} *rs_halos;
struct rs_extra_halo_info {
  int64_t child, next_cochild, prev_cochild;
  int64_t sub_of, ph;
  float max_metric;
} rs_extra_halo;
struct rspart {
  int64_t id;
  float pos[6];
} *rs_particle;


int errorcheck(unsigned int skip1, unsigned int skip2, char *blocklabel){
  if(Debug) printf("checking block %s -- %u  vs  %u\n",blocklabel,skip1,skip2);
  if(skip1 != skip2){
    PyErr_Format(PyExc_IndexError,"skips before and after %s don't match!  %u vs %u",blocklabel,skip1,skip2); 
    return 1;
    //printf("EXITING\n");
    //exit(0);
  }
  else
    return 0;
}


int Nth(double n_to_read, unsigned int total){
  unsigned int nread;

  if(n_to_read)
    nread = ceil((float)total/(float)n_to_read);
  else{
    nread = total;
    nth_Particle = 1;
  }

  if(Debug && Supress==0 && n_to_read)
    printf("particles being read in %d/%d (%0.2f%%)\n",nread,total, \
	   ((float)nread/(float)total) * 100.);
  return nread;
}


int Ngas,Ndm,Ndisk,Nbulge,Nstar,Nbdry,Ntotal;
int Ngas_local,Ndm_local,Ndisk_local,Nbulge_local,Nstar_local,Nbdry_local,Ntotal_local;
//read_header went here!

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
char* NREC  = "nrec";
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
char* RVIR   = "rvir";
char* VVIR   = "vvir";
char* SATSW  = "satswitch";
char* SKIDID = "skidid";

/* TIPSY FUTURE */
char* STARFRAC    = "starfrac";
char* AGESTARFORM = "agestarform";
char* RELAUNCH    = "relaunch";

void assign_type()
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
  else if(strcmp(Values,RVIR)==0)     values = 22;
  else if(strcmp(Values,VVIR)==0)     values = 23;
  else if(strcmp(Values,STARFRAC)==0) values = 24;
  else if(strcmp(Values,AGESTARFORM)==0) values = 25;
  else if(strcmp(Values,RELAUNCH)==0) values = 26;
  else if(strcmp(Values,SATSW)==0)    values = 27;
  else if(strcmp(Values,SKIDID)==0)   values = 28;
  else if(strcmp(Values,NREC)==0)     values = 29;
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
  return;
}


double tconvert;
void init_tconvert(){
  double boltzmann   = 1.380688e-16;   //erg/kelvin
  double proton_mass = 1.67262158e-24;  //grams
  double kmtocm      = 1.e5;
  double gammaminus1 = (5./3.)-1.;
  tconvert     = gammaminus1*(proton_mass/boltzmann)*pow(kmtocm,2);
  return;
 }