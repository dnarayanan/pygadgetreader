#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

#include "modules/funcs.h"
#include "modules/vars.h"
#include "modules/read_header.h"
#ifndef KENCODE
#include "modules/skips.h"
#else
#include "modules/skips_ken.h"
#endif
#include "modules/director.h"

//#ifdef ENABLE_HDF5
//#include <hdf5.h>
//#endif

#define NPY_NO_DEPRECIATED_API NPY_1_7_API_VERSION

// GADGET
#include "modules/gadget/galprop.h"
#include "modules/gadget/read_posvel.h"
#include "modules/gadget/read_gasprops.h"
#include "modules/gadget/read_pid.h"
#include "modules/gadget/read_mass.h"
#include "modules/gadget/read_nspawn.h"
#include "modules/gadget/read_metal.h"
#include "modules/gadget/read_age.h"
#include "modules/gadget/read_tmax.h"
#include "modules/gadget/read_potential.h"
#include "modules/gadget/read_nrec.h"

#ifdef ENABLE_HDF5
// HDF5 GADGET
#include "modules/gadget/read_hdf5.h"
#endif

// TIPSY
#include "modules/tipsy/tipsy_binread.h"
#include "modules/tipsy/tipsy_auxread.h"
#include "modules/tipsy/read_tipsy_only.h"

// ROCKSTAR
#include "modules/misc/read_rockstar.h"

// HDF5 STUFF
//PyObject *hdf5;
//PyObject *pModule,*pDict, *pFunc, *pValue;
//PyObject *pArgs;

void printer(){
    if(Supress==0){
      printf("\ninput (%d files): %s \n",NumFiles,filename);
      printf("extracting %s data for %s\n",Values,Type);
      if(Units==0){
	if(values==5)  printf("returning PHYSICAL density in CODE units\n\n");
	if(values==13) printf("returning PHYSICAL surface density in CODE units\n\n");
	else           printf("returning code units\n\n");
      }
      if(Units==1){
	if(values==3)       printf("returning units of Msun \n\n");
	else if(values==4)  printf("returning units of Kelvin \n\n");
	else if(values==5)  printf("returning PHYSICAL density in units of g/cm^3 \n\n");
	else if(values==13) printf("returning PHYSICAL surface density in units of g/cm^2 \n\n");
	else if(values==15) printf("returning units of Kelvin \n\n");
	else                printf("returning code units\n\n");
      }
    }
}


/*######################### READ HEADER ########################################*/
static PyObject *
readhead(PyObject *self, PyObject *args, PyObject *keywds)
{
  char* simtime   = "time";
  char* redshift  = "redshift";
  char* boxsize   = "boxsize";
  char* O0        = "O0";
  char* Ol        = "Ol";
  char* h         = "h";
  char* gascount  = "gascount";
  char* dmcount   = "dmcount";
  char* diskcount = "diskcount";
  char* bulgecount= "bulgecount";
  char* starcount = "starcount";
  char* bndrycount= "bndrycount";
  char* nfiles    = "num_files";
  //FLAGS
  char* f_sfr     = "f_sfr";
  char* f_fb      = "f_fb";
  char* f_cooling = "f_cooling";
  char* f_age     = "f_age";
  char* f_metals  = "f_metals";
  char* f_fh2     = "f_fh2";
  char* f_dt      = "f_delaytime";
  char* f_tmax    = "f_tmax";
  char* f_pot     = "f_potential";
  char* Value;
  int value;

  j=0;
  NumFiles=1;
  Units=0;
  Debug=0;
  HDF5_FILE = 0;
  
  static char *kwlist[]={"file","value","numfiles","tipsy","debug","hdf5",NULL};
  if(!PyArg_ParseTupleAndKeywords(args,keywds,"ss|iiii",kwlist,&filename,&Value,&NumFiles,&Tipsy,&Debug,&HDF5_FILE)){
    PyErr_Format(PyExc_TypeError,"incorrect input!  must provide filename and value of interest - see readme.txt");
    //return NULL;
  }

#ifndef ENABLE_HDF5
  if(HDF5_FILE==1){
    printf("selected HDF5 when HDF5 is not enabled, please uncomment '#define ENABLE_HDF5' in modules/vars.h\n");
    PyErr_Format(PyExc_IndexError,"HDF5 NOT ENABLED");
    return NULL;
  }
#endif
#if defined(HDF5_DEFAULT) && defined(ENABLE_HDF5)
  HDF5_FILE = 1;
#endif

  if(Debug) printf("metalfactor = %f\n",METALFACTOR);

  int filepresent = 0;
  filepresent = read_header();
  if(Debug) printf("filepresent=%d\n",filepresent);
  
  if(filepresent == 0) {
    printf("file not found...\n");
    PyErr_Format(PyExc_IndexError,"cannot open file : '%s!'",filename);
    return NULL;
  }
  //if the before and after skips don't match return error
  if(filepresent == 2){
    printf("header mismatch...\n");
    PyErr_Format(PyExc_IndexError,"header mismatch in '%s!' (tipsy=%d)",filename,Tipsy);
    return NULL;    
  }
  fclose(infp);

  if(Tipsy==1){
    if(strcmp(Value,simtime)==0)        return Py_BuildValue("d",t_header.time);
    else if(strcmp(Value,redshift)==0)  return Py_BuildValue("d",(1./t_header.time)-1.);
    else if(strcmp(Value,gascount)==0)  return Py_BuildValue("i",t_header.ngas);
    else if(strcmp(Value,dmcount)==0)   return Py_BuildValue("i",t_header.ndark);
    else if(strcmp(Value,starcount)==0) return Py_BuildValue("i",t_header.nstar);
    else
      PyErr_Format(PyExc_TypeError,"incorrect header selection!");
  }
  else{
    if(strcmp(Value,simtime)==0)        return Py_BuildValue("d",header.time);
    else if(strcmp(Value,redshift)==0)  return Py_BuildValue("d",header.redshift);
    else if(strcmp(Value,boxsize)==0)   return Py_BuildValue("d",header.BoxSize);
    else if(strcmp(Value,O0)==0)        return Py_BuildValue("d",header.Omega0);
    else if(strcmp(Value,Ol)==0)        return Py_BuildValue("d",header.OmegaLambda);
    else if(strcmp(Value,h)==0)         return Py_BuildValue("d",header.HubbleParam);
    else if(strcmp(Value,f_sfr)==0)     return Py_BuildValue("i",header.flag_sfr);
    else if(strcmp(Value,f_fb)==0)      return Py_BuildValue("i",header.flag_feedback);
    else if(strcmp(Value,f_cooling)==0) return Py_BuildValue("i",header.flag_cooling);
    else if(strcmp(Value,f_age)==0)     return Py_BuildValue("i",header.flag_stellarage);
    else if(strcmp(Value,f_metals)==0)  return Py_BuildValue("i",header.flag_metals);
    else if(strcmp(Value,gascount)==0)  return Py_BuildValue("i",header.npartTotal[0]);
    else if(strcmp(Value,dmcount)==0)   return Py_BuildValue("i",header.npartTotal[1]);
    else if(strcmp(Value,diskcount)==0) return Py_BuildValue("i",header.npartTotal[2]);
    else if(strcmp(Value,bulgecount)==0)return Py_BuildValue("i",header.npartTotal[3]);
    else if(strcmp(Value,starcount)==0) return Py_BuildValue("i",header.npartTotal[4]);
    else if(strcmp(Value,bndrycount)==0)return Py_BuildValue("i",header.npartTotal[5]);
    else if(strcmp(Value,f_fh2)==0)     return Py_BuildValue("i",header.flag_fH2);
    else if(strcmp(Value,f_tmax)==0)    return Py_BuildValue("i",header.flag_tmax);
    else if(strcmp(Value,f_pot)==0)     return Py_BuildValue("i",header.flag_potential);
    else if(strcmp(Value,f_dt)==0)      return Py_BuildValue("i",header.flag_delaytime);
    else if(strcmp(Value,nfiles)==0)    return Py_BuildValue("i",header.num_files);
    else
      PyErr_Format(PyExc_TypeError,"incorrect header selection!");
  }
}

/*############################# DO SOME WORK! ###########################################*/
//PyArrayObject *array;
static PyObject *
readsnap(PyObject *self, PyObject *args, PyObject *keywds)
{
  j=0;
  NumFiles=1;
  Units=0;
  ERR=0;
  Debug=0;
  Supress=0;
  nth_Particle = 0;
  nth_P = 0.;
  nMetals = 0;
  Future = 0;
  int filepresent = 0;
  HDF5_FILE = 0;

  static char *kwlist[]={"file","data","type","numfiles","units","tipsy","future","debug","supress_output","nth_particle","nmetals","hdf5",NULL};
  if(!PyArg_ParseTupleAndKeywords(args,keywds,"sss|iiiiiifii",kwlist,&filename,&Values,&Type,&NumFiles,&Units,&Tipsy,&Future,&Debug,&Supress,&nth_P,&nMetals,&HDF5_FILE)){
    PyErr_Format(PyExc_TypeError,"wrong input!  must provide filename, data block, and particle type of interest - see readme.txt");
    //return NULL;
  }

#ifndef ENABLE_HDF5
  if(HDF5_FILE==1){
    printf("selected HDF5 when HDF5 is not enabled, please uncomment '#define ENABLE_HDF5' in modules/vars.h\n");
    PyErr_Format(PyExc_IndexError,"HDF5 NOT ENABLED");
    return NULL;
  }
#endif
#if defined(HDF5_DEFAULT) && defined(ENABLE_HDF5)
  HDF5_FILE = 1;
#endif

  nth_Particle = (int)nth_P;

  //if user specifies nMetals then switch it!
  if(nMetals != 0)
    NMETALS = nMetals;

  if(nMetals == 1 && METALFACTOR != 1.0){
    if(Debug)
      printf("switching METALFACTOR from %f to 1.0\n",METALFACTOR);
    METALFACTOR = 1.0;
  }

  if(Debug)
    printf("nth_particle=%d\n",nth_Particle);

  if(Units>1 || Units<0){
    PyErr_Format(PyExc_IndexError,"Units flag must be 0 or 1!");
    //return NULL;
  }

  filepresent = read_header();
  if(Debug) printf("filepresent=%d\n",filepresent);
  
  //if we don't find the file return error
  if(filepresent == 0) {
    printf("file not found...\n");
    PyErr_Format(PyExc_IndexError,"cannot open file : '%s'",filename);
    return NULL;
  }
  //if the before and after skips don't match return error
  if(filepresent == 2){
    printf("header mismatch...\n");
    PyErr_Format(PyExc_IndexError,"header mismatch in '%s!' (tipsy=%d)",filename,Tipsy);
    return NULL;    
  }
  fclose(infp);
  
  assign_type();

  init_tconvert();

  if(nth_Particle){
    nread_total = ceil((float)header.npart[type]/(float)nth_Particle);
    for(j=1;j<NumFiles;j++){
      filepresent = read_header();
      if(filepresent == 0){
	printf("%d file not found...\n",j);
	PyErr_Format(PyExc_IndexError,"cannot open file : '%s!'",filename);
      }
      nread_total += ceil((float)header.npart[type]/(float)nth_Particle);
      fclose(infp);
    }
  }
  else
    nread_total = header.npartTotal[type];

  if(Debug && nth_Particle && Supress==0)
    printf("total particles being read in %d/%d\n",nread_total, header.npartTotal[type]);
  if(nth_Particle && Supress==0)
    printf("READING EVERY %dth PARTICLE\n",nth_Particle);

  if(Tipsy==1){
    Units=1;
    if(Supress==0)
      printf("READING TIPSY\n");
  }

  printer();

  if(HDF5_FILE){
#ifdef ENABLE_HDF5
    read_gadget_HDF5();
#endif
  }
  else if(Future>0){
    read_tipsy_future(Future,values);
  }
  else{

    //protect against "FUTURE" values
    if(values==24 || values==25 || values==26)
      PyErr_Format(PyExc_IndexError,"Only valid for Future files!!  turn future=1 flag on!");

    /*
    if(Supress==0){
      printf("\ninput (%d files): %s \n",NumFiles,filename);
      printf("extracting %s data for %s\n",Values,Type);
      if(Units==0){
	if(values==5)  printf("returning PHYSICAL density in CODE units\n\n");
	if(values==13) printf("returning PHYSICAL surface density in CODE units\n\n");
	else           printf("returning code units\n\n");
      }
      if(Units==1){
	if(values==3)       printf("returning units of Msun \n\n");
	else if(values==4)  printf("returning units of Kelvin \n\n");
	else if(values==5)  printf("returning PHYSICAL density in units of g/cm^3 \n\n");
	else if(values==13) printf("returning PHYSICAL surface density in units of g/cm^2 \n\n");
	else if(values==15) printf("returning units of Kelvin \n\n");
	else                printf("returning code units\n\n");
      }
    }
    */      
    //printf("j=%d\n",j);
    
    //printf("values=%d\n",values);
    if(values==0)       readpos();
    else if(values==1)  readvel();
    else if(values==2)  readpid();
    else if(values==3)  readmass();
    else if(values==4)  readu();
    else if(values==5)  readrho();
    else if(values==6)  readNE();
    else if(values==7)  readNH();
    else if(values==8)  readHSML();
    else if(values==9)  readSFR();
    else if(values==10) readage();
    else if(values==11) readZ();
    else if(values==12) readfh2();
    else if(values==13) readsigma();
    else if(values==14) readmetals();
    else if(values==15) readtmax();
    else if(values==16) readdelaytime();
    else if(values==17) readnspawn();
    else if(values==18) readpotential();
    else if(values==19) read_tipsy();
    else if(values==29) readnrec();
    else if(values==20 || values==21 || values==22 || values==23 || values==27 || values==28) 
      read_tipsy_envira();
    else printf("houston we have a problem...no values returned\n");
    j=0;
  }
  if(Debug) printf("returning array...\n");
  return PyArray_Return(array);
}


static PyObject *
readrockstar(PyObject *self, PyObject *args, PyObject *keywds)
{
  j=0;
  int filepresent = 0;
  Debug = 0;
  Supress = 0;

  static char *kwlist[]={"file","data","debug","supress_output",NULL};
  if(!PyArg_ParseTupleAndKeywords(args,keywds,"ss|ii",kwlist,&filename,&Values,&Debug,&Supress)){
    PyErr_Format(PyExc_TypeError,"wrong input!  must provide filename - see readme.txt");
  }

  char* PARTICLES = "particles";
  char* HALOS     = "halos";

  if(strcmp(Values,HALOS)==0)          values = 0;
  else if(strcmp(Values,PARTICLES)==0) values = 1;
  else PyErr_Format(PyExc_IndexError,"wrong values type selected");

  if(Debug) printf("reading in %s.0.bin\n",filename);

  if(Supress==0){
    printf("reading rockstar file %s.0.bin\n",filename);
    if(values == 0){
      printf("returning halo array:\n");
      printf("ID[0], num_particles[1], mass[2], radius[3], xpos[4], ypos[5], zpos[6], Jx[7], Jy[8], Jz[9]\n");
    }
    if(values == 1){
      printf("returning particle array:\n");
      printf("PID[0], Halo[1]\n");
    }
  }

  rockstar_halos(values);

  if(Debug) printf("returning array...\n");
  return PyArray_Return(array);
}

  
//Initialize Module
PyMethodDef methods[] = {
  {"test",test, METH_VARARGS ,"test function"},
  {"readsnap",readsnap,METH_VARARGS | METH_KEYWORDS, "readsnap info"},
  {"readhead",readhead,METH_VARARGS | METH_KEYWORDS, "read header data"},
  {"readrockstar",readrockstar,METH_VARARGS | METH_KEYWORDS, "read rockstar data"},
  {"readfof",readfof,METH_VARARGS | METH_KEYWORDS, "read fof data"},
  {"galprop",galprop,METH_VARARGS | METH_KEYWORDS, "read galaxy property data"},
  {"galdata",galdata,METH_VARARGS | METH_KEYWORDS, "return galaxy particle info"},
  {"galdataonly",galdataonly,METH_VARARGS | METH_KEYWORDS, "return galaxy particle info (no indexes)"},
  {NULL,NULL,0,NULL}
};


PyMODINIT_FUNC
initreadgadget()
{
  (void) Py_InitModule("readgadget",methods);
  import_array();
}
