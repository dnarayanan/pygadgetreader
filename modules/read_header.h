#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

#ifdef ENABLE_HDF5
#include <hdf5.h>
#endif

int read_header()
{
  
  /*
  if(NumFiles>1) 
    sprintf(infile,"%s.%d",filename,j);
  else 
    sprintf(infile,"%s",filename);
  */

#ifdef ENABLE_HDF5
  if(HDF5_FILE){
    int status = 0;

    hid_t hdf5_file, hdf5_headergrp, hdf5_attribute;
    
    hdf5_file      = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    hdf5_headergrp = H5Gopen1(hdf5_file, "/Header");

    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Time");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
    H5Aclose(hdf5_attribute);
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Redshift");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.redshift);
    H5Aclose(hdf5_attribute);
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "BoxSize");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.BoxSize);
    H5Aclose(hdf5_attribute);
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Omega0");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.Omega0);
    H5Aclose(hdf5_attribute);
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "OmegaLambda");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.OmegaLambda);
    H5Aclose(hdf5_attribute);
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "HubbleParam");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.HubbleParam);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_Sfr");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_sfr);
    H5Aclose(hdf5_attribute);
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_Cooling");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_cooling);
    H5Aclose(hdf5_attribute);
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_StellarAge");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_stellarage);
    H5Aclose(hdf5_attribute);
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_Feedback");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_feedback);
    H5Aclose(hdf5_attribute);
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_Metals");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_metals);
    H5Aclose(hdf5_attribute);

    if(H5Aexists(hdf5_headergrp, "Flag_fH2")){
      hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_fH2");
      H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_fH2);
      H5Aclose(hdf5_attribute);
    }
    else
      header.flag_fH2 = 0;
    if(H5Aexists(hdf5_headergrp, "Flag_Potential")){
      hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_Potential");
      H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_potential);
      H5Aclose(hdf5_attribute);
    }
    else
      header.flag_potential = 0;
    if(H5Aexists(hdf5_headergrp, "Flag_DelayTime")){
      hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_DelayTime");
      H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_delaytime);
      H5Aclose(hdf5_attribute);
    }
    else
      header.flag_delaytime = 0;
    if(H5Aexists(hdf5_headergrp, "Flag_TMax")){
      hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_TMax");
      H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_tmax);
      H5Aclose(hdf5_attribute);
    }
    else
      header.flag_tmax = 0;

    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumFilesPerSnapshot");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.num_files);
    H5Aclose(hdf5_attribute);
    
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, header.npartTotal);
    H5Aclose(hdf5_attribute);
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_ThisFile");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, header.npart);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "MassTable");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass);
    H5Aclose(hdf5_attribute);

    H5Gclose(hdf5_headergrp);
    H5Fclose(hdf5_file);

    if(header.flag_metals == 1)
      METALFACTOR = 1.0;

    header.num_files = 1;

    // Assign TOTAL particle counts
    Ngas   = header.npartTotal[0];
    Ndm    = header.npartTotal[1];
    Ndisk  = header.npartTotal[2];
    Nbulge = header.npartTotal[3];
    Nstar  = header.npartTotal[4];
    Nbdry  = header.npartTotal[5];
    Ntotal = Ngas+Ndm+Ndisk+Nbulge+Nstar+Nbdry;
    return 1;
  }
#endif


  if(Debug) printf("reading file...ERR=%d\n",ERR);

  sprintf(infile,"%s",filename);
  if(!(infp=fopen(infile,"r"))){
    sprintf(infile,"%s.%d",filename,j);
    if(!(infp=fopen(infile,"r"))) {
      ERR = 1;
      PyErr_Format(PyExc_TypeError,"can't open file : '%s'",infile);
      //return NULL;
    }
  }
  /*
  sprintf(infile,"%s",filename);
  if( access( filename, F_OK ) != -1){
    infp=fopen(infile,"r");
  }
  else{
    sprintf(infile,"%s.%d",filename,j);
    if( access( filename, F_OK) != -1){
      infp=fopen(infile,"r");
    }
    else
      ERR = 1;    
  }
  */
  
  if(Debug) printf("after file open...ERR=%d\n",ERR);
  
  if(ERR){
    if(Debug) printf("lol in error statement...\n");
    return 0;
  }

  if(Debug) printf("after ERR statement...\n");

  /*
  if(!(infp=fopen(infile,"r"))){
    ERR = 1;
    PyErr_Format(PyExc_TypeError,"can't open file: '%s'",infile);
    //return NULL;
  }
  */

  if(Tipsy==1){
    fread(&t_header, sizeof(t_header),1,infp);
    MAXDIM = t_header.ndim;
    Ngas   = t_header.ngas;
    Ndm    = t_header.ndark;
    Nstar  = t_header.nstar;
    Ntotal = t_header.ntotal;
  }
  else{
    int skip1, skip2;
    char* blocklabel;

    blocklabel="Header";

    rewind(infp);
    // READ HEADER
    fread(&skip1,sizeof(int),1,infp);
    fread(&header,sizeof(header),1,infp);
    fread(&skip2,sizeof(int),1,infp);

    if(header.flag_metals == 1)
      METALFACTOR = 1.0;

    if(Debug)
      printf("HEADER: skip1=%d  skip2=%d\n",skip1,skip2);

    if(errorcheck(skip1,skip2,blocklabel)){
       PyErr_Format(PyExc_IndexError,"skips before and after %s don't match!  %d vs %d",blocklabel,skip1,skip2); 
       return 2;
    }
        
    //printf("numfiles=%d \t\t header.num_files=%d \n",NumFiles,header.num_files);
    
    NumFiles = header.num_files;
    
    if(NumFiles != header.num_files){
      PyErr_Format(PyExc_IndexError,"NumFiles(%d) != header.num_files(%d)!",NumFiles,header.num_files);
      //return NULL;
    }
    if(header.flag_metals > 0 && NMETALS != header.flag_metals)
      printf("WARNING: NMETALS(%d) != header.flag_metals(%d)!\n",NMETALS,header.flag_metals);
      //PyErr_Format(PyExc_IndexError,"NMETALS(%d) != header.flag_metals(%d)!",NMETALS,header.flag_metals);

    
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
  return 1;
}
