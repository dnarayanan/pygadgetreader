#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

const char *filename;
FILE *infp;
char infile[500];
int ERR;
int NumFiles, Units, j, dummy;

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
  int SFRZ  = 1;
  int HIH2  = 0;

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
  char* SFRMETAL= "sfrmetal";
  char* HIMASS  = "HImass";
  char* H2MASS  = "H2mass";

  double convert = 1e10;

  static char *kwlist[]={"dir","snapnumber","value","units","sfrZ","HIH2",NULL};
  if(!PyArg_ParseTupleAndKeywords(args,keywds,"sis|iii",kwlist,&Directory,&Snap,&Value,&Units,&SFRZ,&HIH2)){
    PyErr_Format(PyExc_TypeError,"incorrect input!  must provide properties file directory, snap number, and value you're interested in - see readme.txt");
    //return NULL;
  }

  sprintf(infile,"%s/properties_%03d",Directory,Snap);
  if(!(infp=fopen(infile,"r"))){
    PyErr_Format(PyExc_IOError,"can't open file: '%s'",infile);
    //return NULL;
  }

  if(HIH2==1 && SFRZ == 0)
    SFRZ = 1;

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
  else if(strcmp(Value,SFRMETAL)==0 && SFRZ==1) value = 10;
  else if(strcmp(Value,HIMASS)==0 && HIH2==1) value = 11;
  else if(strcmp(Value,H2MASS)==0 && HIH2==1) value = 12;
  else{
    PyErr_Format(PyExc_IndexError,"wrong values type selected");
    //return NULL;
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
    float sfrmetal;
    float HIgalmass;
    float H2galmass;
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
  if(value==9) printf("Returning %d group gas MASS-weighted metallicities\n",Ngroup);
  if(value==10) printf("Returning %d group gas SFR-weighted metallicities\n",Ngroup);
  if(value==11){
    if(Units==0) printf("Returning %d group HI mass in code units\n",Ngroup);
    if(Units==1) printf("Returning %d group HI mass in units of Msun\n",Ngroup);
  }
  if(value==12){
    if(Units==0) printf("Returning %d group H2 mass in code units\n",Ngroup);
    if(Units==1) printf("Returning %d group H2 mass in units of Msun\n",Ngroup);
  } 

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
    if(SFRZ==1)
      fread(&gal[i].sfrmetal, sizeof(float),1,infp);
    if(HIH2==1){
      fread(&gal[i].HIgalmass, sizeof(float),1,infp);
      fread(&gal[i].H2galmass, sizeof(float),1,infp);
    }
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
    if(value==10) for(i=0;i<Ngroup;i++) MDATA(array,i) = gal[i].sfrmetal;
    if(value==11){
      if(Units==0) for(i=0;i<Ngroup;i++) MDATA(array,i) = gal[i].HIgalmass;
      if(Units==1) for(i=0;i<Ngroup;i++) MDATA(array,i) = gal[i].HIgalmass*convert;
    }
    if(value==12){
      if(Units==0) for(i=0;i<Ngroup;i++) MDATA(array,i) = gal[i].H2galmass;
      if(Units==1) for(i=0;i<Ngroup;i++) MDATA(array,i) = gal[i].H2galmass*convert;
    } 
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
    //return NULL;
  }
  
  if(NumFiles>1) 
    sprintf(infile,"%s.%d",filename,j);
  else 
    sprintf(infile,"%s",filename);
  if(!(infp=fopen(infile,"r"))){
    ERR = 1;
    PyErr_Format(PyExc_TypeError,"can't open file: '%s'",infile);
    //return NULL;
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


/*###################### GALDATA ########################*/
static PyObject *
galdata(PyObject *self, PyObject *args, PyObject *keywds)
{
  PyObject *array;
  //PyObject *array2;
  char *Directory;
  int Snap;
  int ndim = 2;
  FILE *fp_pos, *fp_index, *fp_type, *fp_cat, *fp_prop;
  int tot_gal, len;
  int k, i;
  int galnum;
  float x,y,z;
  unsigned int id;
  int Numfiles;
  double tempIndex;  

  static char *kwlist[]={"file","dir","snapnumber","galnum",NULL};
  if(!PyArg_ParseTupleAndKeywords(args,keywds,"ssii",kwlist,&filename,&Directory,&Snap,&galnum)){
    PyErr_Format(PyExc_TypeError,"incorrect input!  snapshot, galprop data file directory, snap number, and galaxy number - see readme.txt");
    //return NULL;
  }

  //READ SNAPSHOT
  sprintf(infile,"%s",filename);
  if(!(infp=fopen(infile,"r"))){
    sprintf(infile,"%s.0",filename);
    if(!(infp=fopen(infile,"r"))) {
      ERR = 1;
      PyErr_Format(PyExc_TypeError,"can't open file : '%s'",infile);
      //return NULL;
    }
  }

  // READ HEADER
  Skip;
  fread(&header,sizeof(header),1,infp);
  Skip;
  fclose(infp);
  NumFiles = header.num_files;


  unsigned int *s_tallies, *g_tallies, *num_gas, *num_dm, *num_star;
  unsigned int *num_disk, *num_bulge, *num_bndry;
  s_tallies=(unsigned int*)malloc(NumFiles*sizeof(unsigned int));
  g_tallies=(unsigned int*)malloc(NumFiles*sizeof(unsigned int));
  num_gas  =(unsigned int*)malloc(NumFiles*sizeof(unsigned int));
  num_dm   =(unsigned int*)malloc(NumFiles*sizeof(unsigned int));
  num_disk =(unsigned int*)malloc(NumFiles*sizeof(unsigned int));
  num_bulge=(unsigned int*)malloc(NumFiles*sizeof(unsigned int));
  num_star =(unsigned int*)malloc(NumFiles*sizeof(unsigned int));
  num_bndry=(unsigned int*)malloc(NumFiles*sizeof(unsigned int));

  s_tallies[0] = header.npart[0]+header.npart[1]+header.npart[2]+header.npart[3];
  g_tallies[0] = 0;
  num_gas[0]   = header.npart[0];
  num_dm[0]    = header.npart[1];
  num_disk[0]  = header.npart[2];
  num_bulge[0] = header.npart[3];
  num_star[0]  = header.npart[4];
  num_bndry[0] = header.npart[5];

  if(NumFiles>1){
    for(i=1;i<NumFiles;i++){
      sprintf(infile,"%s.%d",filename,i);
      if(!(infp=fopen(infile,"r"))){
        ERR=1;
        PyErr_Format(PyExc_TypeError, "can't open file : '%s'",infile);
      }

      Skip;
      fread(&header,sizeof(header),1,infp);
      Skip;
      fclose(infp);

      num_gas[i]   = header.npart[0];
      num_dm[i]    = header.npart[1];
      num_disk[i]  = header.npart[2];
      num_bulge[i] = header.npart[3];
      num_star[i]  = header.npart[4];
      num_bndry[i] = header.npart[5];
    }
    for(i=1;i<NumFiles;i++){
      /*
      g_tallies[i] = g_tallies[i-1] + num_dm[i-1]  + num_star[i-1];
      s_tallies[i] = s_tallies[i-1] + num_gas[i] + num_dm[i];
      */
      g_tallies[i] = g_tallies[i-1] + num_dm[i-1]  + num_disk[i-1] + num_bulge[i-1] + num_star[i-1] + num_bndry[i-1];
      s_tallies[i] = s_tallies[i-1] + num_bndry[i-1] + num_gas[i] + num_dm[i] + num_disk[i] + num_bulge[i];
    }
    /*
    for(i=0;i<NumFiles;i++){
      printf("num_gas[%d]=%d\n",i,num_gas[i]);
      fflush(stdout);
    }
    for(i=0;i<NumFiles;i++){
      printf("g_tallies[%d]=%d\ts_tallies[%d]=%d\n",i,g_tallies[i],i,s_tallies[i]);
      fflush(stdout);
    }
    */
  }

  sprintf(infile,"%s/pos_%03d",Directory,Snap);
  if(!(fp_pos=fopen(infile,"rb"))){
    PyErr_Format(PyExc_IOError,"can't open file: '%s'",infile);
    //return NULL;
  }

  sprintf(infile,"%s/index_%03d",Directory,Snap);
  if(!(fp_index=fopen(infile,"rb"))){
    PyErr_Format(PyExc_IOError,"can't open file: '%s'",infile);
    //return NULL;
  }
  sprintf(infile,"%s/type_%03d",Directory,Snap);
  if(!(fp_type=fopen(infile,"rb"))){
    PyErr_Format(PyExc_IOError,"can't open file: '%s'",infile);
    //return NULL;
  }
  sprintf(infile,"%s/catalogue_%03d",Directory,Snap);
  if(!(fp_cat=fopen(infile,"rb"))){
    PyErr_Format(PyExc_IOError,"can't open file: '%s'",infile);
    //return NULL;
  }

  /*
  sprintf(infile,"%s/properties_%03d",Directory,Snap);
  if(!(fp_prop=fopen(infile,"rb"))){
    PyErr_Format(PyExc_IOError,"can't open file: '%s'",infile);
    return NULL;
  }
  */

  fseek( fp_pos,   sizeof(int), SEEK_CUR);
  fseek( fp_index, sizeof(int), SEEK_CUR);
  fseek( fp_type,  sizeof(int), SEEK_CUR);
  fread( &tot_gal, sizeof(int), 1, fp_cat);
  //fseek( fp_prop,  sizeof(int), SEEK_CUR);
  
  //printf("  total number of galaxies: %d\n", tot_gal);

  for(i=0; i<tot_gal; i++){
    fread( &len,         sizeof(int), 1, fp_cat);
    fseek( fp_cat,       sizeof(int), SEEK_CUR);
    //fseek( fp_prop, 12*sizeof(float), SEEK_CUR);

    if (i==galnum){
      //printf("particles in target galaxy %d: %d\n",galnum,len);
      npy_intp dims[2]={len,5};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);

      for(k=0; k<len; k++){
        fread( &x,    sizeof(float), 1, fp_pos);
        fread( &y,    sizeof(float), 1, fp_pos);
        fread( &z,    sizeof(float), 1, fp_pos);
        fread( &id,   sizeof(int),   1, fp_index);
        fread( &type, sizeof(int),   1, fp_type);

        tempIndex = 0.;
        tempIndex = (double)return_index(NumFiles,s_tallies,g_tallies,id,type);
        //printf("tempIndex=%f\n",tempIndex);

        DATA(array, k, 0) = x;
        DATA(array, k, 1) = y;
        DATA(array, k, 2) = z;
        DATA(array, k, 3) = tempIndex;
        DATA(array, k, 4) = (double)type;
      }
    }
    else{
      for(k=0; k<len; k++){
        fseek( fp_pos, 3*sizeof(float), SEEK_CUR);
        fseek( fp_index,   sizeof(int), SEEK_CUR);
        fseek( fp_type,    sizeof(int), SEEK_CUR);
      }
    }  
  }
  fclose(fp_pos);
  fclose(fp_index);
  fclose(fp_type);
  fclose(fp_cat);
  //fclose(fp_prop);

  return PyArray_Return(array);
}

return_index(int NumFiles, unsigned int *s_tallies, unsigned int *g_tallies, unsigned int id, int type)
{
  int i;
  id = id - 1;
  unsigned int new_id;
  int escape = 0;

  if(NumFiles > 1){
    for(i=0;i<NumFiles-1;i++){
      if(type==4){
        if(id < s_tallies[i+1] && escape==0){
          new_id = id - s_tallies[i];
          escape = 1;
        }
        if(id > s_tallies[NumFiles-1] && escape==0){
          new_id = id - s_tallies[NumFiles-1];
          escape = 1;
        }
      }
      else if(type==0){
        for(i=0;i<NumFiles;i++){
          if(id<s_tallies[i] && escape==0){
            new_id = id - g_tallies[i];
            escape = 1;
          }
        }
      }
    }
  }
  else if(NumFiles==1){
    if(type==4)
      new_id = id - s_tallies[0];
    if(type==0)
      new_id = id;
  }
  else{
    printf("%d detected!  original ID=%d, type=%d\n",new_id, id, type);
    fflush(stdout);
  }
  return new_id;
}

/*###################### GALDATAONLY NO INDEXES ########################*/
static PyObject *
galdataonly(PyObject *self, PyObject *args, PyObject *keywds)
{
  PyObject *array;
  //PyObject *array2;
  char *Directory;
  int Snap;
  int ndim = 2;
  FILE *fp_pos, *fp_index, *fp_type, *fp_cat, *fp_prop;
  int tot_gal, len;
  int k, i;
  int galnum;
  float x,y,z;
  unsigned int id;
  int Numfiles;

  static char *kwlist[]={"dir","snapnumber","galnum",NULL};
  if(!PyArg_ParseTupleAndKeywords(args,keywds,"sii",kwlist,&Directory,&Snap,&galnum)){
    PyErr_Format(PyExc_TypeError,"incorrect input!  galprop data file directory, snap number, and galaxy number - see readme.txt");
    //return NULL;
  }
  
  sprintf(infile,"%s/pos_%03d",Directory,Snap);
  if(!(fp_pos=fopen(infile,"rb"))){
    PyErr_Format(PyExc_IOError,"can't open file: '%s'",infile);
    //return NULL;
  }

  sprintf(infile,"%s/index_%03d",Directory,Snap);
  if(!(fp_index=fopen(infile,"rb"))){
    PyErr_Format(PyExc_IOError,"can't open file: '%s'",infile);
    //return NULL;
  }
  sprintf(infile,"%s/type_%03d",Directory,Snap);
  if(!(fp_type=fopen(infile,"rb"))){
    PyErr_Format(PyExc_IOError,"can't open file: '%s'",infile);
    //return NULL;
  }
  sprintf(infile,"%s/catalogue_%03d",Directory,Snap);
  if(!(fp_cat=fopen(infile,"rb"))){
    PyErr_Format(PyExc_IOError,"can't open file: '%s'",infile);
    //return NULL;
  }

  /*
  sprintf(infile,"%s/properties_%03d",Directory,Snap);
  if(!(fp_prop=fopen(infile,"rb"))){
    PyErr_Format(PyExc_IOError,"can't open file: '%s'",infile);
    return NULL;
  }
  */

  fseek( fp_pos,   sizeof(int), SEEK_CUR);
  fseek( fp_index, sizeof(int), SEEK_CUR);
  fseek( fp_type,  sizeof(int), SEEK_CUR);
  fread( &tot_gal, sizeof(int), 1, fp_cat);
  //fseek( fp_prop,  sizeof(int), SEEK_CUR);
  
  //printf("  total number of galaxies: %d\n", tot_gal);

  for(i=0; i<tot_gal; i++){
    fread( &len,         sizeof(int), 1, fp_cat);
    fseek( fp_cat,       sizeof(int), SEEK_CUR);
    //fseek( fp_prop, 12*sizeof(float), SEEK_CUR);

    if (i==galnum){
      //printf("particles in target galaxy %d: %d\n",galnum,len);
      npy_intp dims[2]={len,5};
      array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);

      for(k=0; k<len; k++){
        fread( &x,    sizeof(float), 1, fp_pos);
        fread( &y,    sizeof(float), 1, fp_pos);
        fread( &z,    sizeof(float), 1, fp_pos);
        fread( &id,   sizeof(int),   1, fp_index);
        fread( &type, sizeof(int),   1, fp_type);

        DATA(array, k, 0) = x;
        DATA(array, k, 1) = y;
        DATA(array, k, 2) = z;
        DATA(array, k, 3) = 0;
        DATA(array, k, 4) = (double)type;
      }
    }
    else{
      for(k=0; k<len; k++){
        fseek( fp_pos, 3*sizeof(float), SEEK_CUR);
        fseek( fp_index,   sizeof(int), SEEK_CUR);
        fseek( fp_type,    sizeof(int), SEEK_CUR);
      }
    }  
  }
  fclose(fp_pos);
  fclose(fp_index);
  fclose(fp_type);
  fclose(fp_cat);
  //fclose(fp_prop);

  return PyArray_Return(array);
}


static PyObject *
test(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args,"si",&filename,&NumFiles)){
    PyErr_Format(PyExc_TypeError,"wrong input");
    //return NULL;
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
