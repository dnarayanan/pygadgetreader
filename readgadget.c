#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>
#define Skip fread(&dummy,sizeof(dummy),1,infp)
#define DATA(a,i,j)*((float *) PyArray_GETPTR2(a,i,j))
#define PIDDATA(a,i)*((int *) PyArray_GETPTR1(a,i))
#define MDATA(a,i)*((float *) PyArray_GETPTR1(a,i))

const char *filename;  
FILE *infp;
char infile[200];
int NumFiles;
int dummy;

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

  char     fill[56];  /* fills to 256 Bytes */
} header;

int j;
int Ngas,Ndm,Ndisk,Nbulge,Nstar,Nbdry,Ntotal;
int Ngas_local,Ndm_local,Ndisk_local,Nbulge_local,Nstar_local,Nbdry_local,Ntotal_local;
read_header()
{
  if(NumFiles>1) sprintf(infile,"%s.%d",filename,j);
  else sprintf(infile,"%s",filename);
  if(!(infp=fopen(infile,"r"))){
    PyErr_Format(PyExc_IOError,"cannot open '%s' !",infile);
    return NULL;
  }
  // READ HEADER
  Skip;
  fread(&header,sizeof(header),1,infp);
  Skip;
  
  if(NumFiles != header.num_files){
    PyErr_Format(PyExc_IndexError,"NumFiles != header.num_files!");
    return NULL;
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
  printf("%d\n",Ngas_local);
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
    PyErr_Format(PyExc_IndexError,"wrong partile type selected");
    return NULL;
  }
  if(type==0 && Ngas==0){
    PyErr_Format(PyExc_IndexError,"No %s particles!",gas);
    return NULL;
  }
  if(type==1 && Ndm==0){
    PyErr_Format(PyExc_IndexError,"No %s particles!",dm);
    return NULL;
  }
  if(type==2 && Ndisk==0){
    PyErr_Format(PyExc_IndexError,"No %s particles!",disk);
    return NULL;
  }
  if(type==3 && Nbulge==0){
    PyErr_Format(PyExc_IndexError,"No %s particles!",bulge);
    return NULL;
  }
  if(type==4 && Nstar==0){
    PyErr_Format(PyExc_IndexError,"No %s!",stars);
    return NULL;
  }
  if(type==5 && Nbdry==0){
    PyErr_Format(PyExc_IndexError,"No %s particles!",bndry);
    return NULL;
  }
  /*
  if(type==0) printf("%d gas selected\n",Ngas);  
  if(type==1) printf("%d dm selected\n",Ndm);
  if(type==2) printf("%d disk selected\n",Ndisk);
  if(type==3) printf("%d bulge selected\n",Nbulge);
  if(type==4) printf("%d stars selected\n",Nstar);
  if(type==5) printf("%d bndry selected\n",Nbdry);
  */
}

/*######################### POS ########################################*/
static PyObject *
readpos(PyObject *self, PyObject *args)
{  
  PyArrayObject *array;
  float *data;
  int ndim = 2;

  if(!PyArg_ParseTuple(args,"sis",&filename,&NumFiles,&Type)){
    PyErr_Format(PyExc_TypeError,"incorrect number of arguments - correct syntax is (filename,# of Files,'particle type' (gas,dm,disk,bulge,stars,bndry)");
    return NULL;
  }
  int i;
  int n;
  int pc = 0;
  //  int index = 0;
  //int count = 0;
  for(j=0;j<NumFiles;j++){
    read_header();
    if(j==0){
      assign_type();
      npy_intp dims[2]={header.npartTotal[type],3};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_FLOAT);
    }
    data=(float*)malloc(header.npart[type]*sizeof(float)*3);
    
    Skip; //skip before POS
    //seek past particle groups not interested in
    for(i=1;i<=type;i++){
      fseek(infp,header.npart[i-1]*3*sizeof(float),SEEK_CUR);
    }
    fread(data,header.npart[type]*sizeof(float)*3,1,infp);
    Skip; //skip after POS
    fclose(infp);

    //count = count + header.npart[type];
    for(n=0;n<header.npart[type];n++)
      {
	DATA(array,pc,0) = data[3*n];
	DATA(array,pc,1) = data[3*n+1];
	DATA(array,pc,2) = data[3*n+2];
	pc++;
      }
  }
  return PyArray_Return(array);
}

/*######################### VEL ########################################*/
static PyObject *
readvel(PyObject *self, PyObject *args)
{  
  PyArrayObject *array;
  float *data;
  int ndim = 2;

  if(!PyArg_ParseTuple(args,"sis",&filename,&NumFiles,&Type)){
    PyErr_Format(PyExc_TypeError,"incorrect number of arguments - correct syntax is (filename,# of Files,'particle type' (gas,dm,disk,bulge,stars,bndry)");
    return NULL;
  }
  int i;
  int n;
  int pc = 0;
  //  int index = 0;
  //int count = 0;
  for(j=0;j<NumFiles;j++){
    read_header();
    if(j==0){
      assign_type();
      npy_intp dims[2]={header.npartTotal[type],3};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_FLOAT);
    }
    data=(float*)malloc(header.npart[type]*sizeof(float)*3);
    
    //skip positions
    Skip;
    fseek(infp,Ntotal_local*sizeof(float)*3,SEEK_CUR);
    Skip;

    Skip; //skip before VEL
    //seek past particle groups not interested in
    for(i=1;i<=type;i++){
      fseek(infp,header.npart[i-1]*3*sizeof(float),SEEK_CUR);
    }
    fread(data,header.npart[type]*sizeof(float)*3,1,infp);
    Skip; //skip after VEL
    fclose(infp);

    //count = count + header.npart[type];
    for(n=0;n<header.npart[type];n++)
      {
	DATA(array,pc,0) = data[3*n];
	DATA(array,pc,1) = data[3*n+1];
	DATA(array,pc,2) = data[3*n+2];
	pc++;
      }
  }
  return PyArray_Return(array);
}

/*######################### PID ########################################*/
static PyObject *
readpid(PyObject *self, PyObject *args)
{  
  PyArrayObject *array;
  int *data;
  int ndim = 1;

  if(!PyArg_ParseTuple(args,"sis",&filename,&NumFiles,&Type)){
    PyErr_Format(PyExc_TypeError,"incorrect number of arguments - correct syntax is (filename,# of Files,'particle type' (gas,dm,disk,bulge,stars,bndry)");
    return NULL;
  }
  int i;
  int n;
  int pc = 0;
  //  int index = 0;
  //int count = 0;
  for(j=0;j<NumFiles;j++){
    read_header();
    if(j==0){
      assign_type();
      npy_intp dims[1]={header.npartTotal[type]};
      //npy_intp dims[1]={1};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_INT);
    }
    data=(int*)malloc(header.npart[type]*sizeof(int));
    
    //skip positions
    Skip;
    fseek(infp,Ntotal_local*sizeof(float)*3,SEEK_CUR);
    Skip;

    //skip velocities
    Skip;
    fseek(infp,Ntotal_local*sizeof(float)*3,SEEK_CUR);
    Skip;

    Skip; //skip before PID
    //seek past particle groups not interested in
    for(i=1;i<=type;i++){
      fseek(infp,header.npart[i-1]*sizeof(int),SEEK_CUR);
    }
    fread(data,header.npart[type]*sizeof(int),1,infp);
    Skip; //skip after PID
    fclose(infp);

    //count = count + header.npart[type];
    for(n=0;n<header.npart[type];n++)
      {
	PIDDATA(array,pc) = data[n];
	pc++;
      }
  }
  return PyArray_Return(array);
}


/*######################### MASSES ########################################*/
static PyObject *
readmass(PyObject *self, PyObject *args)
{  
  PyArrayObject *array;
  float *data;
  int ndim = 1;

  if(!PyArg_ParseTuple(args,"sis",&filename,&NumFiles,&Type)){
    PyErr_Format(PyExc_TypeError,"incorrect number of arguments - correct syntax is (filename,# of Files,'particle type' (gas,dm,disk,bulge,stars,bndry)");
    return NULL;
  }
  int i;
  int n;
  int pc = 0;
  //  int index = 0;
  //int count = 0;
  for(j=0;j<NumFiles;j++){
    read_header();
    if(j==0){
      assign_type();
      npy_intp dims[1]={header.npartTotal[type]};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_FLOAT);
    }

    if(header.mass[type]>0){
      for(n=0;n<header.npart[type];n++){
	  MDATA(array,pc)=header.mass[type];
	  pc++;
	}
      printf("non-zero header mass detected - using header mass\n");
      return PyArray_Return(array);
    }
    else{
    data=(float*)malloc(header.npart[type]*sizeof(float));

    //skip positions
    Skip;
    fseek(infp,Ntotal_local*sizeof(float)*3,SEEK_CUR);
    Skip;

    //skip velocities
    Skip;
    fseek(infp,Ntotal_local*sizeof(float)*3,SEEK_CUR);
    Skip;
    
    //skip PIDs
    Skip;
    fseek(infp,Ntotal_local*sizeof(int),SEEK_CUR);
    Skip;

    Skip; //skip before MASS
    //seek past particle groups not interested in
    for(i=1;i<=type;i++){
      if(header.mass[i]==0.0){
      fseek(infp,header.npart[i-1]*sizeof(float),SEEK_CUR);
      }
    }
    fread(data,header.npart[type]*sizeof(float),1,infp);
    Skip; //skip after MASS
    fclose(infp);

    //count = count + header.npart[type];
    for(n=0;n<header.npart[type];n++)
      {
	MDATA(array,pc) = data[n];
	pc++;
      }
  }
  return PyArray_Return(array);
  }
}


/*######################### INTERNAL ENERGY ########################################*/
static PyObject *
readu(PyObject *self, PyObject *args)
{  
  PyArrayObject *array;
  float *data;
  int ndim = 1;

  if(!PyArg_ParseTuple(args,"sis",&filename,&NumFiles,&Type)){
    PyErr_Format(PyExc_TypeError,"incorrect number of arguments - correct syntax is (filename,# of Files,'particle type' (gas,dm,disk,bulge,stars,bndry)");
    return NULL;
  }
  int i;
  int n;
  int pc = 0;
  //  int index = 0;
  //int count = 0;
  for(j=0;j<NumFiles;j++){
    read_header();
    if(Ngas==0){
      PyErr_Format(PyExc_IndexError,"Ngas=0 - No internal energy to read!");
      return NULL;
    }
    if(j==0){
      assign_type();
      if(type!=0){
	PyErr_Format(PyExc_IndexError,"Internal Energy can only be read for gas!!");
	return NULL;
      }
      npy_intp dims[1]={header.npartTotal[type]};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_FLOAT);
    }

    data=(float*)malloc(header.npart[type]*sizeof(float));

    //skip positions
    Skip;
    fseek(infp,Ntotal_local*sizeof(float)*3,SEEK_CUR);
    Skip;

    //skip velocities
    Skip;
    fseek(infp,Ntotal_local*sizeof(float)*3,SEEK_CUR);
    Skip;
    
    //skip PIDs
    Skip;
    fseek(infp,Ntotal_local*sizeof(int),SEEK_CUR);
    Skip;

    //skip MASS 
    Skip;
    for(i=0;i<6;i++){
      if(header.npart[i]>0 && header.mass[i]==0.0){
	fseek(infp,header.npart[i]*sizeof(float),SEEK_CUR);
      }
    }
    Skip;

    Skip; //skip before U
    fread(data,header.npart[type]*sizeof(float),1,infp);
    Skip; //skip after U
    fclose(infp);

    //count = count + header.npart[type];
    for(n=0;n<header.npart[type];n++)
      {
	MDATA(array,pc) = data[n];
	pc++;
      }
  }
  return PyArray_Return(array);
}


PyMethodDef methods[] = {
  //  {"readheader",readheader,METH_VARARGS, "reads gadget snap header"},
  {"readpos",readpos,METH_VARARGS, "reads particle position data from gadget snapshot"},
  {"readvel",readvel,METH_VARARGS, "reads particle velocity data from gadget snapshot"},
  {"readpid",readpid,METH_VARARGS, "reads particle ID data from gadget snapshot"},
  {"readmass",readmass,METH_VARARGS, "reads particle mass data from gadget snapshot"},
  {"readu",readu,METH_VARARGS, "reads particle internal energy data from gadget snapshot"}
};

PyMODINIT_FUNC
initreadgadget()
{
  (void) Py_InitModule("readgadget",methods);
  import_array();
}
