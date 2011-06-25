#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>
#define Skip fread(&dummy,sizeof(dummy),1,infp)
#define DATA(a,i,j)*((double *) PyArray_GETPTR2(a,i,j))
#define PIDDATA(a,i)*((int *) PyArray_GETPTR1(a,i))
#define MDATA(a,i)*((double *) PyArray_GETPTR1(a,i))

const char *filename;  
FILE *infp;
char infile[500];
int NumFiles, Units, j, dummy;
int ERR;

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


int Ngas,Ndm,Ndisk,Nbulge,Nstar,Nbdry,Ntotal;
int Ngas_local,Ndm_local,Ndisk_local,Nbulge_local,Nstar_local,Nbdry_local,Ntotal_local;
read_header()
{
  if(NumFiles>1) 
    sprintf(infile,"%s.%d",filename,j);
  else 
    sprintf(infile,"%s",filename);
  if(!(infp=fopen(infile,"r"))){
    ERR = 1;
    PyErr_Format(PyExc_TypeError,"can't open file: '%s'",infile);
    return NULL;
  }
  // READ HEADER
  Skip;
  fread(&header,sizeof(header),1,infp);
  Skip;
  
  //printf("numfiles=%d \t\t header.num_files=%d \n",NumFiles,header.num_files);

  if(NumFiles != header.num_files){
    PyErr_Format(PyExc_IndexError,"NumFiles(%d) != header.num_files(%d)!",NumFiles,header.num_files);
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
char* Values;
int values;

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
  // match requested data type and assign corresponding integer
  if(strcmp(Values,POS)==0)        values = 0;
  else if(strcmp(Values,VEL)==0)   values = 1;
  else if(strcmp(Values,PID)==0)   values = 2;
  else if(strcmp(Values,MASS)==0)  values = 3;
  else if(strcmp(Values,U)==0)     values = 4;
  else if(strcmp(Values,RHO)==0)   values = 5;
  else if(strcmp(Values,NE)==0)    values = 6;
  else if(strcmp(Values,NH)==0)    values = 7;
  else if(strcmp(Values,HSML)==0)  values = 8;
  else if(strcmp(Values,SFR)==0)   values = 9;
  else if(strcmp(Values,AGE)==0)   values = 10;
  else if(strcmp(Values,Z)==0)     values = 11;
  else{
    PyErr_Format(PyExc_IndexError,"wrong values type selected");
    return NULL;
  }

  //  printf("type=%d \t values=%d\n",type,values);

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
  if(type==0) printf("%d gas selected, extracting %s data\n",Ngas,Values);  
  if(type==1) printf("%d dm selected, extracting %s data\n",Ndm,Values);
  if(type==2) printf("%d disk selected, extracting %s data\n",Ndisk,Values);
  if(type==3) printf("%d bulge selected, extracting %s data\n",Nbulge,Values);
  if(type==4) printf("%d stars selected, extracting %s data\n",Nstar,Values);
  if(type==5) printf("%d bndry selected, extracting %s data\n",Nbdry,Values);
  */
}

/*###################### GALPROP ########################*/
static PyObject *
galprop(PyObject *self, PyObject *args, PyObject *keywds)
{
  PyObject *array;
  int i;
  int Ngroup;
  char *Directory;
  int Snap;
  char *Value;
  int value;
  int Units = 0;

  char* MSTAR   = "mstar";
  char* BMAG    = "bmag";
  char* IMAG    = "imag";
  char* VMAG    = "vmag";
  char* KMAG    = "kmag";
  char* CM      = "cm";
  char* SFR     = "sfr";
  char* MGAS    = "mgas";
  char* SMETAL  = "zstar";
  char* GMETAL  = "zgas";

  double convert = 1e10;

  static char *kwlist[]={"dir","snapnumber","value","units",NULL};
  if(!PyArg_ParseTupleAndKeywords(args,keywds,"sis|i",kwlist,&Directory,&Snap,&Value,&Units)){
    PyErr_Format(PyExc_TypeError,"incorrect input!  must provide properties file directory, snap number, and value you're interested in - see readme.txt");
    return NULL;
  }

  sprintf(infile,"%s/properties_%03d",Directory,Snap);
  if(!(infp=fopen(infile,"r"))){
    PyErr_Format(PyExc_IOError,"can't open file: '%s'",infile);
    return NULL;
  }

  if(strcmp(Value,MSTAR)==0)       value = 0;
  else if(strcmp(Value,BMAG)==0)   value = 1;
  else if(strcmp(Value,IMAG)==0)   value = 2;
  else if(strcmp(Value,VMAG)==0)   value = 3;
  else if(strcmp(Value,KMAG)==0)   value = 4;
  else if(strcmp(Value,CM)==0)     value = 5;
  else if(strcmp(Value,SFR)==0)    value = 6;
  else if(strcmp(Value,MGAS)==0)   value = 7;
  else if(strcmp(Value,SMETAL)==0) value = 8;
  else if(strcmp(Value,GMETAL)==0) value = 9;
  else{
    PyErr_Format(PyExc_IndexError,"wrong values type selected");
    return NULL;
  }

  printf("\n\nReading %s \n",infile);

  struct gal_data
  {
    float mstar;
    float Bmag;
    float Imag;
    float Vmag;
    float Kmag;
    float cm[3];
    float sfr;
    float mgas;
    float metalstar;
    float metalgas;
  } *gal;
  
  fread(&Ngroup,sizeof(int),1,infp);
  //printf("Ngroup=%d\n",Ngroup);
  if(value==0){
    if(Units==0) printf("Returning %d stellar group masses in code units\n",Ngroup);
    if(Units==1) printf("Returning %d stellar group masses in units of Msun\n",Ngroup);
}
  if(value==1) printf("Returning %d group B-magnitudes\n",Ngroup);
  if(value==2) printf("Returning %d group I-magnitudes\n",Ngroup);
  if(value==3) printf("Returning %d group V-magnitudes\n",Ngroup);
  if(value==4) printf("Returning %d group K-magnitudes\n",Ngroup);
  if(value==5) printf("Returning %d group center of mass positions\n",Ngroup);
  if(value==6) printf("Returning %d group SFR in Msun/year\n",Ngroup);
  if(value==7){
    if(Units==0) printf("Returning %d group gas masses in code units\n",Ngroup);
    if(Units==1) printf("Returning %d group gas masses units of Msun\n",Ngroup);
  }
  if(value==8) printf("Returning %d group stellar metallicities\n",Ngroup);
  if(value==9) printf("Returning %d group gas metallicities\n",Ngroup);

  gal=malloc(sizeof(struct gal_data)*Ngroup);

  for(i=0;i<Ngroup;i++){
    fread(&gal[i].mstar,    sizeof(float),1,infp);
    fread(&gal[i].Bmag,     sizeof(float),1,infp);
    fread(&gal[i].Imag,     sizeof(float),1,infp);
    fread(&gal[i].Vmag,     sizeof(float),1,infp);
    fread(&gal[i].Kmag,     sizeof(float),1,infp);
    fread(&gal[i].cm,       sizeof(float),3,infp);
    fread(&gal[i].sfr,      sizeof(float),1,infp);
    fread(&gal[i].mgas,     sizeof(float),1,infp);
    fread(&gal[i].metalstar,sizeof(float),1,infp);
    fread(&gal[i].metalgas, sizeof(float),1,infp);
  }
  fclose(infp);

  if(value==5){
    int ndim = 2;
    npy_intp dims[2] = {Ngroup,3};
    array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    
    for(i=0;i<Ngroup;i++){
      DATA(array,i,0) = gal[i].cm[0];
      DATA(array,i,1) = gal[i].cm[1];
      DATA(array,i,2) = gal[i].cm[2];
    }
  }
  else{
    int ndim = 1;
    npy_intp dims[1]={Ngroup};
    array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);

    if(value==0){
      if(Units==0) for(i=0;i<Ngroup;i++) MDATA(array,i) = gal[i].mstar;
      if(Units==1) for(i=0;i<Ngroup;i++) MDATA(array,i) = gal[i].mstar*convert;
    }
    if(value==1) for(i=0;i<Ngroup;i++) MDATA(array,i) = gal[i].Bmag;
    if(value==2) for(i=0;i<Ngroup;i++) MDATA(array,i) = gal[i].Imag;
    if(value==3) for(i=0;i<Ngroup;i++) MDATA(array,i) = gal[i].Vmag;
    if(value==4) for(i=0;i<Ngroup;i++) MDATA(array,i) = gal[i].Kmag;
    if(value==6) for(i=0;i<Ngroup;i++) MDATA(array,i) = gal[i].sfr;
    if(value==7){
      if(Units==0) for(i=0;i<Ngroup;i++) MDATA(array,i) = gal[i].mgas;
      if(Units==1) for(i=0;i<Ngroup;i++) MDATA(array,i) = gal[i].mgas*convert;
    }
    if(value==8) for(i=0;i<Ngroup;i++) MDATA(array,i) = gal[i].metalstar;
    if(value==9) for(i=0;i<Ngroup;i++) MDATA(array,i) = gal[i].metalgas;
  }
  return PyArray_Return(array);
}


static PyObject *
readfof(PyObject *self, PyObject *args, PyObject *keywds)
{
  int i;
  int Ngroup;

  static char *kwlist[]={"file","numfiles",NULL};
  if(!PyArg_ParseTupleAndKeywords(args,keywds,"s|i",kwlist,&filename,&NumFiles)){
    PyErr_Format(PyExc_TypeError,"incorrect input!  must provide filename and value of interest - see readme.txt");
    return NULL;
  }
  
  if(NumFiles>1) 
    sprintf(infile,"%s.%d",filename,j);
  else 
    sprintf(infile,"%s",filename);
  if(!(infp=fopen(infile,"r"))){
    ERR = 1;
    PyErr_Format(PyExc_TypeError,"can't open file: '%s'",infile);
    return NULL;
  }

  struct count{
    int total, cum, gas, dm, star;
  } *cnt;

  struct groupmass{
    float total,gas,dm,star;
  } *gmass;

  struct center_of_mass{
    float xpos,ypos,zpos;
  } *cm;

  fread(&Ngroup,sizeof(int),1,infp);
  
  cnt   = malloc(sizeof(struct count) * Ngroup);
  gmass = malloc(sizeof(struct groupmass) * Ngroup);
  cm    = malloc(sizeof(struct center_of_mass) * Ngroup);

  for(i=0;i<Ngroup;i++){
    cnt[i].total=0, cnt[i].cum=0,cnt[i].gas=0,cnt[i].dm=0,cnt[i].star=0;
    gmass[i].total=0.,gmass[i].gas=0.,gmass[i].dm=0.,gmass[i].star=0.;
    cm[i].xpos=0.,cm[i].ypos=0.,cm[i].zpos=0.;
  }

  for(i=0;i<Ngroup;i++) fread(&cnt[i].total,sizeof(int),1,infp);
  for(i=0;i<Ngroup;i++) fread(&cnt[i].cum,sizeof(int),1,infp);
  for(i=0;i<Ngroup;i++) fread(&gmass[i].total,sizeof(float),1,infp);

  for(i=0;i<Ngroup;i++){
    fread(&cm[i].xpos,sizeof(float),1,infp);
    fread(&cm[i].ypos,sizeof(float),1,infp);
    fread(&cm[i].zpos,sizeof(float),1,infp);
  }

  for(i=0;i<Ngroup;i++){
    fread(&cnt[i].gas,sizeof(int),1,infp);
    fread(&cnt[i].dm,sizeof(int),1,infp);
    fread(&cnt[i].star,sizeof(int),1,infp);
  }

  for(i=0;i<Ngroup;i++){
    fread(&gmass[i].gas,sizeof(float),1,infp);
    fread(&gmass[i].dm,sizeof(float),1,infp);
    fread(&gmass[i].star,sizeof(float),1,infp);
  }
  fclose(infp);


 
//  for(i=0;i<Ngroup;i++){
//    printf("%d\n",cnt[i].dm);
    //printf("%f\n",gmass[i].dm);
//  }


  //  simdata=(float*)malloc(header.npart[type]*sizeof(float)*3);
  //  fread(simdata,header.npart[type]*sizeof(float)*3,1,infp);


  printf("ngroups=%d!\n",Ngroup);
  return Py_None;

}


static PyObject *
test(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args,"si",&filename,&NumFiles)){
    PyErr_Format(PyExc_TypeError,"wrong input");
    return NULL;
  }
  j=0;
  /* 
 if(NumFiles>1) 
    sprintf(infile,"%s.%d",filename,j);
  else 
    sprintf(infile,"%s",filename);
  if(!(infp=fopen(infile,"r"))){
    printf("cannot open file!!\n");
    
    PyErr_Format(PyExc_IOError,"shitbags");
    PyErr_Format(PyExc_IndexError,"cannot open '%s' !",infile);
    printf("after PyErr_Format\n");
    return NULL;
    
    PyErr_Format(PyExc_IOError,"gadgetPyIO module can't open file: '%s'",infile);
    return NULL;
    printf("after return null\n");
  }
  */
  read_header();
  return Py_None;
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
  //FLAGS
  char* f_sfr     = "f_sfr";
  char* f_fb      = "f_fb";
  char* f_cooling = "f_cooling";
  char* f_age     = "f_age";
  char* f_metals  = "f_metals";
  char* Value;
  int value;

  j=0;
  NumFiles=1;
  Units=0;

  static char *kwlist[]={"file","value","numfiles",NULL};
  if(!PyArg_ParseTupleAndKeywords(args,keywds,"ss|i",kwlist,&filename,&Value,&NumFiles)){
    PyErr_Format(PyExc_TypeError,"incorrect input!  must provide filename and value of interest - see readme.txt");
    return NULL;
  }
  read_header();
  fclose(infp);  

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
}

/*############################# DO SOME WORK! ###########################################*/
PyArrayObject *array;
static PyObject *
readsnap(PyObject *self, PyObject *args, PyObject *keywds)
{
  j=0;
  NumFiles=1;
  Units=0;
  ERR=0;

  static char *kwlist[]={"file","data","type","numfiles","units",NULL};
  if(!PyArg_ParseTupleAndKeywords(args,keywds,"sss|ii",kwlist,&filename,&Values,&Type,&NumFiles,&Units)){
    PyErr_Format(PyExc_TypeError,"wrong input!  must provide filename, data block, and particle type of interest - see readme.txt");
    return NULL;
  }

  if(Units>1 || Units<0){
    PyErr_Format(PyExc_IndexError,"Units flag must be 0 or 1!");
    return NULL;
  }

  read_header();
  if(ERR==1){
    PyErr_Format(PyExc_TypeError,"readsnap: can't open file: '%s' ERR=%d",infile,ERR);
    return NULL;
  }
  if(NumFiles!=header.num_files){
    PyErr_Format(PyExc_IndexError,"NumFiles(%d) != header.num_files(%d)!",NumFiles,header.num_files);
    return NULL;
  }

  fclose(infp);
  assign_type();

  printf("\ninput: %s \n",filename);
  printf("extracting %s data for %s\n",Values,Type);
  if(Units==0) printf("returning code units\n\n");
  if(Units==1){
    if(values==3) printf("returning units of Msun \n\n");
    else if(values==4) printf("returning units of Kelvin \n\n");
    else if(values==5) printf("returning units of g/cm^3 \n\n");
    else printf("returning code units\n\n");
  }
  
  //printf("j=%d\n",j);

  //  printf("values=%d\n",values);
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
  else printf("houston we have a problem...no values returned\n");
  j=0;
  return PyArray_Return(array);
}

/*######################### POS ########################################*/
readpos()
{  
  float *simdata;
  int ndim = 2;
  
  int i;
  int n;
  int pc = 0;
  for(j=0;j<NumFiles;j++){
    read_header();
    if(j==0){
      npy_intp dims[2]={header.npartTotal[type],3};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    }
    simdata=(float*)malloc(header.npart[type]*sizeof(float)*3);
    
    Skip; //skip before POS
    //seek past particle groups not interested in
    for(i=1;i<=type;i++){
      fseek(infp,header.npart[i-1]*3*sizeof(float),SEEK_CUR);
    }
    fread(simdata,header.npart[type]*sizeof(float)*3,1,infp);
    Skip; //skip after POS
    fclose(infp);

    //count = count + header.npart[type];
    for(n=0;n<header.npart[type];n++)
      {
	DATA(array,pc,0) = simdata[3*n];
	DATA(array,pc,1) = simdata[3*n+1];
	DATA(array,pc,2) = simdata[3*n+2];
	pc++;
      }
  }
  if(pc!=header.npartTotal[type]){
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");
    return NULL;
  }
}

/*######################### VEL ########################################*/

readvel()
{  
  float *simdata;
  int ndim = 2;

  int i;
  int n;
  int pc = 0;
  for(j=0;j<NumFiles;j++){
    read_header();
    if(j==0){
      npy_intp dims[2]={header.npartTotal[type],3};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    }
    simdata=(float*)malloc(header.npart[type]*sizeof(float)*3);
    
    skippos();

    Skip; //skip before VEL
    //seek past particle groups not interested in
    for(i=1;i<=type;i++){
      fseek(infp,header.npart[i-1]*3*sizeof(float),SEEK_CUR);
    }
    fread(simdata,header.npart[type]*sizeof(float)*3,1,infp);
    Skip; //skip after VEL
    fclose(infp);

    //count = count + header.npart[type];
    for(n=0;n<header.npart[type];n++)
      {
	if(Units==1){
	  DATA(array,pc,0) = simdata[3*n]   * sqrt(header.time);
	  DATA(array,pc,1) = simdata[3*n+1] * sqrt(header.time);
	  DATA(array,pc,2) = simdata[3*n+2] * sqrt(header.time);
	  pc++;
	}
	else{
	  DATA(array,pc,0) = simdata[3*n];
	  DATA(array,pc,1) = simdata[3*n+1];
	  DATA(array,pc,2) = simdata[3*n+2];
	  pc++;
	}
      }
  }
  if(pc!=header.npartTotal[type]){
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");
    return NULL;
  }
}

/*######################### PID ########################################*/
readpid()
{  
  int *simdata;
  int ndim = 1;

  int i;
  int n;
  int pc = 0;
  for(j=0;j<NumFiles;j++){
    read_header();
    if(j==0){
      npy_intp dims[1]={header.npartTotal[type]};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_INT);
    }
    simdata=(int*)malloc(header.npart[type]*sizeof(int));
    
    skippos();
    skipvel();

    Skip; //skip before PID
    //seek past particle groups not interested in
    for(i=1;i<=type;i++){
      fseek(infp,header.npart[i-1]*sizeof(int),SEEK_CUR);
    }
    fread(simdata,header.npart[type]*sizeof(int),1,infp);
    Skip; //skip after PID
    fclose(infp);

    //count = count + header.npart[type];
    for(n=0;n<header.npart[type];n++)
      {
	PIDDATA(array,pc) = simdata[n];
	pc++;
      }
  }
  if(pc!=header.npartTotal[type]){
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");
    return NULL;
  }
}


/*######################### MASSES ########################################*/
readmass()
{  
  float *simdata;
  int ndim = 1;

  int i;
  int n;
  int pc = 0;
  for(j=0;j<NumFiles;j++){
    read_header();
    if(j==0){
      npy_intp dims[1]={header.npartTotal[type]};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    }
    /*
    if(Units==1){
      printf("### RETURNING MASS IN GADGET UNITS, MULTIPLY BY 1.98892e43 TO CONVER TO GRAMS ###\n");
    }
    */

    if(header.mass[type]>0 && header.npart[type]>0){
      printf("non-zero header mass detected - using header mass for %s\n",Type);
      for(n=0;n<header.npart[type];n++)
	{
	  if(Units==0) MDATA(array,pc)=header.mass[type];
	  if(Units==1) MDATA(array,pc)=header.mass[type]*1e10;//1.98892e43;
	  pc++;
	}
    }
    else{
      printf("reading mass block for %s\n",Type);
      simdata=(float*)malloc(header.npart[type]*sizeof(float));
      
      skippos();
      skipvel();
      skippid();
      
      Skip; //skip before MASS
      //seek past particle groups not interested in
      if(type==4) fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
      fread(simdata,header.npart[type]*sizeof(float),1,infp);
      Skip; //skip after MASS
      fclose(infp);
      
      //count = count + header.npart[type];
      for(n=0;n<header.npart[type];n++)
	{
	  if(Units==0) MDATA(array,pc) = simdata[n];
	  if(Units==1) MDATA(array,pc) = simdata[n]*1e10;//1.98892e43;
	  pc++;
	}
    }
  }
  if(pc!=header.npartTotal[type]){
    PyErr_Format(PyExc_IndexError,"particle count mismatch! pc=%d  npartTotal[%d]=%d",pc,type,header.npartTotal[type]);
    return NULL;
  }
}

/*######################### INTERNAL ENERGY ########################################*/
readu()
{  
  float *simdata;
  int ndim  = 1;
  
  double convert;
  if(Units==1){
    double boltzmann   = 1.380688e-16;   //erg/kelvin
    double proton_mass = 1.67262158e-24;  //grams
    double kmtocm      = 1.e5;
    double gammaminus1 = (5./3.)-1.;
    convert     = gammaminus1*(proton_mass/boltzmann)*pow(kmtocm,2);
  }
  int i;
  int n;
  int pc = 0; 
  for(j=0;j<NumFiles;j++){
    read_header();
    if(Ngas==0){
      PyErr_Format(PyExc_IndexError,"Ngas=0 - No internal energy to read!");
      return NULL;
    }
    if(j==0){
      if(type!=0){
	PyErr_Format(PyExc_IndexError,"Internal Energy can only be read for gas!!");
	return NULL;
      }
      npy_intp dims[1]={header.npartTotal[type]};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    }

    simdata=(float*)malloc(header.npart[type]*sizeof(float));

    skippos();
    skipvel();
    skippid();
    skipmass();
 
    Skip; //skip before U
    fread(simdata,header.npart[type]*sizeof(float),1,infp);
    Skip; //skip after U
    fclose(infp);

    //count = count + header.npart[type];
    for(n=0;n<header.npart[type];n++)
      {
	if(Units==0) MDATA(array,pc) = simdata[n];
	if(Units==1) MDATA(array,pc) = simdata[n]*convert;
	pc++;
      }
  }
  if(pc!=header.npartTotal[type]){
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");
    return NULL;
  }
}


/*######################### RHO ########################################*/
readrho()
{  
  float *simdata;
  int ndim  = 1;

  double convert;
  if(Units==1){
    //values used for converting from code units to cgs
    double solarmass = 1.98892e33;
    double kpctocm   = 3.08568025e21;
    convert   = (pow(10.,10)*solarmass)/pow(kpctocm,3);
  }

  int i;
  int n;
  int pc = 0;
  for(j=0;j<NumFiles;j++){
    read_header();
    if(Ngas==0){
      PyErr_Format(PyExc_IndexError,"Ngas=0 - No density to read!");
      return NULL;
    }
    if(j==0){
      if(type!=0){
	PyErr_Format(PyExc_IndexError,"Density can only be read for gas!!");
	return NULL;
      }
      npy_intp dims[1]={header.npartTotal[type]};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    }

    simdata=(float*)malloc(header.npart[type]*sizeof(float));

    skippos();
    skipvel();
    skippid();
    skipmass();
    skipu();

    Skip; //skip before RHO
    fread(simdata,header.npart[type]*sizeof(float),1,infp);
    Skip; //skip after RHO
    fclose(infp);

    //count = count + header.npart[type];
    for(n=0;n<header.npart[type];n++)
      {
	if(Units==0) MDATA(array,pc) = simdata[n];
	if(Units==1) MDATA(array,pc) = simdata[n]*convert;
	pc++;
      }
  }
  if(pc!=header.npartTotal[type]){
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");
    return NULL;
  }
}



/*######################### NE ########################################*/
readNE()
{  
  float *simdata;
  int ndim = 1;

  int i;
  int n;
  int pc = 0;
  for(j=0;j<NumFiles;j++){
    read_header();
    if(Ngas==0){
      PyErr_Format(PyExc_IndexError,"Ngas=0 - No NE to read!");
      return NULL;
    }
    if(j==0){
      if(type!=0){
	PyErr_Format(PyExc_IndexError,"NE can only be read for gas!!");
	return NULL;
      }
      npy_intp dims[1]={header.npartTotal[type]};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    }

    simdata=(float*)malloc(header.npart[type]*sizeof(float));

    skippos();
    skipvel();
    skippid();
    skipmass();
    skipu();
    skiprho();

    Skip; //skip before NE
    fread(simdata,header.npart[type]*sizeof(float),1,infp);
    Skip; //skip after NE
    fclose(infp);

    //count = count + header.npart[type];
    for(n=0;n<header.npart[type];n++)
      {
	MDATA(array,pc) = simdata[n];
	pc++;
      }
  }
  if(pc!=header.npartTotal[type]){
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");
    return NULL;
  }
}


/*######################### NH ########################################*/
readNH()
{  
  float *simdata;
  int ndim = 1;

  int i;
  int n;
  int pc = 0;
  for(j=0;j<NumFiles;j++){
    read_header();
    if(Ngas==0){
      PyErr_Format(PyExc_IndexError,"Ngas=0 - No NH to read!");
      return NULL;
    }
    if(j==0){
      if(type!=0){
	PyErr_Format(PyExc_IndexError,"NH can only be read for gas!!");
	return NULL;
      }
      npy_intp dims[1]={header.npartTotal[type]};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    }

    simdata=(float*)malloc(header.npart[type]*sizeof(float));

    skippos();
    skipvel();
    skippid();
    skipmass();
    skipu();
    skiprho();
    skipne();

    Skip; //skip before NH
    fread(simdata,header.npart[type]*sizeof(float),1,infp);
    Skip; //skip after NH
    fclose(infp);

    //count = count + header.npart[type];
    for(n=0;n<header.npart[type];n++)
      {
	MDATA(array,pc) = simdata[n];
	pc++;
      }
  }
  if(pc!=header.npartTotal[type]){
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");
    return NULL;
  }
}


/*######################### HSML ########################################*/
readHSML()
{  
  float *simdata;
  int ndim = 1;

  int i;
  int n;
  int pc = 0;
  for(j=0;j<NumFiles;j++){
    read_header();
    if(Ngas==0){
      PyErr_Format(PyExc_IndexError,"Ngas=0 - No HSML to read!");
      return NULL;
    }
    if(j==0){
      if(type!=0){
	PyErr_Format(PyExc_IndexError,"HSML can only be read for gas!!");
	return NULL;
      }
      npy_intp dims[1]={header.npartTotal[type]};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    }

    simdata=(float*)malloc(header.npart[type]*sizeof(float));

    skippos();
    skipvel();
    skippid();
    skipmass();
    skipu();
    skiprho();
    skipne();
    skipnh();

    Skip; //skip before HSML
    fread(simdata,header.npart[type]*sizeof(float),1,infp);
    Skip; //skip after HSML
    fclose(infp);

    //count = count + header.npart[type];
    for(n=0;n<header.npart[type];n++)
      {
	MDATA(array,pc) = simdata[n];
	pc++;
      }
  }
  if(pc!=header.npartTotal[type]){
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");
    return NULL;
  }
}


/*######################### SFR ########################################*/
readSFR()
{  
  float *simdata;
  int ndim = 1;

  int i;
  int n;
  int pc = 0;
  for(j=0;j<NumFiles;j++){
    read_header();
    if(Ngas==0){
      PyErr_Format(PyExc_IndexError,"Ngas=0 - No SFR to read!");
      return NULL;
    }
    if(j==0){
      if(type!=0){
	PyErr_Format(PyExc_IndexError,"SFR can only be read for gas!!");
	return NULL;
      }
      npy_intp dims[1]={header.npartTotal[type]};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    }

    simdata=(float*)malloc(header.npart[type]*sizeof(float));

    skippos();
    skipvel();
    skippid();
    skipmass();
    skipu();
    skiprho();
    skipne();
    skipnh();
    skiphsml();

    Skip; //skip before SFR
    fread(simdata,header.npart[type]*sizeof(float),1,infp);
    Skip; //skip after SFR
    fclose(infp);

    //count = count + header.npart[type];
    for(n=0;n<header.npart[type];n++)
      {
	MDATA(array,pc) = simdata[n];
	pc++;
      }
  }
  if(pc!=header.npartTotal[type]){
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");
    return NULL;
  }
}


/*######################### AGE ########################################*/
readage()
{  
  float *simdata;
  int ndim = 1;
  int i;
  int n;
  int pc = 0;
  for(j=0;j<NumFiles;j++){
    read_header();
    if(header.flag_stellarage==0){
      PyErr_Format(PyExc_IndexError,"flag_stellarage=%d --> AGE NOT TRACKED",header.flag_stellarage);
      return NULL;
    }
    if(Nstar==0){
      PyErr_Format(PyExc_IndexError,"Nstar=0 - No stars to read!");
      return NULL;
    }
    if(j==0){
      
      if(type!=4){
	PyErr_Format(PyExc_IndexError,"Age can only be read for stars!!");
	return NULL;
      }
      npy_intp dims[1]={header.npartTotal[type]};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    }

    simdata=(float*)malloc(header.npart[type]*sizeof(float));

    skippos();
    skipvel();
    skippid();
    skipmass();
    skipu();
    skiprho();
    skipne();
    skipnh();
    skiphsml();
    skipsfr();

    Skip; //skip before AGE
    fread(simdata,header.npart[type]*sizeof(float),1,infp);
    Skip; //skip after AGE
    fclose(infp);

    //count = count + header.npart[type];
    for(n=0;n<header.npart[type];n++)
      {
	MDATA(array,pc) = simdata[n];
	pc++;
      }
  }
  if(pc!=header.npartTotal[type]){
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");
    return NULL;
  }
}


/*######################### Z ########################################*/
readZ()
{  
  float *simdata;
  int ndim = 1;

  int i;
  int n;
  int pc = 0;
  for(j=0;j<NumFiles;j++){
    read_header();
    if(header.flag_metals==0){
      PyErr_Format(PyExc_IndexError,"flag_metals=%d --> METALS NOT TRACKED",header.flag_metals);
      return NULL;
    }
    if(Ngas==0 && Nstar==0 || header.flag_metals==0){
      if(Ngas==0 && Nstar==0) PyErr_Format(PyExc_IndexError,"Nstar=0 and Ngas=0 - No metallicity to read!");
      if(header.flag_metals==0) PyErr_Format(PyExc_IndexError,"flag_metals=%d - No metallicity to read!",header.flag_metals);
      return NULL;
    }
    if(j==0){
      if(type!=0 && type!=4){
	PyErr_Format(PyExc_IndexError,"Z can only be read for gas or stars!!");
	return NULL;
      }
      npy_intp dims[1]={header.npartTotal[type]};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
    }

    simdata=(float*)malloc(header.npart[type]*sizeof(float));

    skippos();
    skipvel();
    skippid();
    skipmass();
    skipu();
    skiprho();
    skipne();
    skipnh();
    skiphsml();
    skipsfr();
    skipage();

    Skip; //skip before Z
    if(type==0){
      fread(simdata,header.npart[type]*sizeof(float),1,infp);
    }
    else{
      fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
      fread(simdata,header.npart[type]*sizeof(float),1,infp);
    }
    Skip; //skip after Z
    fclose(infp);

    //count = count + header.npart[type];
    for(n=0;n<header.npart[type];n++)
      {
	MDATA(array,pc) = simdata[n];
	pc++;
      }
  }
  if(pc!=header.npartTotal[type]){
    PyErr_Format(PyExc_IndexError,"particle count mismatch!");
    return NULL;
  }
}

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
  if(header.mass[0]==0 && header.npart[0]>0 || header.mass[4]==0 && header.npart[4]>0){
    Skip;
    if(header.mass[0]==0 && header.npart[0]>0) fseek(infp,header.npart[0]*sizeof(float),SEEK_CUR);
    if(header.mass[4]==0 && header.npart[0]>0) fseek(infp,header.npart[4]*sizeof(float),SEEK_CUR);
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

//Initialize Module
PyMethodDef methods[] = {
  {"test",test, METH_VARARGS ,"test function"},
  {"readsnap",readsnap,METH_VARARGS | METH_KEYWORDS, "readsnap info"},
  {"readhead",readhead,METH_VARARGS | METH_KEYWORDS, "read header data"},
  {"readfof",readfof,METH_VARARGS | METH_KEYWORDS, "read fof data"},
  {"galprop",galprop,METH_VARARGS | METH_KEYWORDS, "read galaxy property data"},
  {NULL,NULL,0,NULL}
};


PyMODINIT_FUNC
initreadgadget()
{
  (void) Py_InitModule("readgadget",methods);
  import_array();
}
