#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

void rockstar_halos(int values)
{
  
  int64_t *dmlist;
  int64_t *pids;
  int64_t i,j,k;
  int64_t pid;
  
  int ndim=2;
  int64_t n_particles = 0;
  int64_t n_halos = 0;

  int breaker = 0;

  int nread   = 5000;
  int n_files = 0;
  
  int hcounter = 0;
  int pcounter = 0;
  int hstart   = 0;

  //count number of files, halos, and particles
  for(j=0; j<nread; j++){    

    //see if we can open the jth file
    sprintf(infile,"%s.%d.bin",filename,j);
    if(!(infp=fopen(infile,"r"))){
      sprintf(infile,"%s.%d",filename,j);
      if(!(infp=fopen(infile,"r"))) {
	breaker = 1;
	break;
      }
    }
    
    n_files ++;
    fread(&rs_header, sizeof(rs_header),1,infp);

    n_particles += rs_header.num_particles;
    n_halos     += rs_header.num_halos;

    if(Debug)
      printf("found %d particles and %d halos in file %d\n",
	     rs_header.num_particles,rs_header.num_halos,j);

    fclose(infp);
  }
  /*
  if(breaker)
    fclose(infp);
  */
  nread = n_files;

  //allocate memory
  rs_halos = (struct rshalo*)malloc(sizeof(struct rshalo) * n_halos);
  dmlist   = (int64_t*)malloc(sizeof(int64_t) * n_particles);
  pids     = (int64_t*)malloc(sizeof(int64_t) * n_particles);

  //particles
  if(values==1){
    npy_intp dims[2]={n_particles,2};
    array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_UINT32);
  }
  //halos
  else if(values==0){
    npy_intp dims[2]={n_halos,10};
    array = (PyArrayObject *)PyArray_SimpleNew(ndim,dims,PyArray_DOUBLE);
  }


  //cycle through files
  for(j=0; j<nread; j++){

    if(Debug) printf("opening %s.%d.bin\n",filename,j);
    sprintf(infile,"%s.%d.bin",filename,j);
    if(!(infp=fopen(infile,"r"))){
      sprintf(infile,"%s.%d",filename,j);
      if(!(infp=fopen(infile,"r"))) {
	ERR = 1;
	PyErr_Format(PyExc_TypeError,"can't open file : '%s'",infile);
      }
    }
    
    if(Debug) printf("reading header...\n");
    fread(&rs_header, sizeof(rs_header),1,infp);

    if(Debug) printf("cycling through halos...\n");
    //read the next halo struct into the appropriate index
    hstart = hcounter;
    for(i=0; i<rs_header.num_halos; i++){
      fread(&rs_halos[hcounter],sizeof(struct rshalo),1,infp);
      hcounter += 1;
    }
    
    if(Debug) printf("cycling through halos to read in particle data\n");
    //cycle through halos and read in particle data
    for(i=0; i<rs_header.num_halos; i++){	       
      for(k=0; k<rs_halos[hstart+i].num_p; k++){
	//printf("index %d   nump=%e halo %d\n",k,(double)rs_halos[hstart+1].num_p,i);
	fread(&pid, sizeof(int64_t),1,infp);
	pids[pcounter]   = pid;
	dmlist[pcounter] = rs_halos[hstart+i].id;
	pcounter += 1;
      }
    }
    if(Debug) printf("closing file...\n");
    fclose(infp);
  }

  /*
  if(Debug){
    for(i=0; i<n_particles; i++)
      printf("ID:%d  HID:%d\n",pids[i],dmlist[i]);
  }
  */

  printf("int64_t: %d\n",sizeof(int64_t));

  printf("%d\n",rs_halos[0].id);
  for(i=0; i<6; i++)
    printf(" %f\n",rs_halos[0].pos[i]);
  for(i=0;i<3;i++)
    printf(" %f\n",rs_halos[0].corevel[i]);
  for(i=0;i<3;i++)
    printf(" %f\n",rs_halos[0].bulkvel[i]);
  
  printf("\n");

  printf("%e\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n",
	 rs_halos[0].m,
	 rs_halos[0].r,
	 rs_halos[0].child_r,
	 rs_halos[0].vmax_r,
	 rs_halos[0].mgrav,
	 rs_halos[0].vmax,
	 rs_halos[0].rvmax,
	 rs_halos[0].rs,
	 rs_halos[0].klypin_rs,
	 rs_halos[0].vrms);

  printf("\n");
	 
  for(i=0;i<3;i++)
    printf(" %f\n",rs_halos[0].J[i]);
  
  printf("%f\n%f\n",
	 rs_halos[0].energy,
	 rs_halos[0].spin);

  for(i=0;i<4;i++)
    printf(" %f\n",rs_halos[0].alt_m[i]);

  printf("\n");

  printf("%f\n%f\n%f\n%f\n",
	 rs_halos[0].Xoff,
	 rs_halos[0].Voff,
	 rs_halos[0].b_to_a,
	 rs_halos[0].c_to_a);

  for(i=0;i<3;i++)
    printf("%f\n",rs_halos[0].A[i]);
  printf("%f\n%f\n",
	 rs_halos[0].b_to_a2,
	 rs_halos[0].c_to_a2);
  for(i=0;i<3;i++)
    printf("%f\n",rs_halos[0].A2[i]);

  printf("%f\n%f\n%f\n%f\n",
	 rs_halos[0].bullock_spin,
	 rs_halos[0].kin_to_pot,
	 rs_halos[0].m_pe_b,
	 rs_halos[0].m_pe_d);

  printf("\n");

  printf("%d\n%d\n%d\n%d\n%d\n%d\n%f\n%f\n%f\n",
	 rs_halos[0].num_p,
	 rs_halos[0].num_child_particles,
	 rs_halos[0].p_start,
	 rs_halos[0].desc,
	 rs_halos[0].flags,
	 rs_halos[0].n_core,
	 rs_halos[0].min_pos_err,
	 rs_halos[0].min_vel_err,
	 rs_halos[0].min_bulkvel_err);

  //assign values to the correct retrun array
  int pc = 0;
  if(values==0){
    for(i=0; i<n_halos; i++){
      DATA(array,pc,0) = (double)rs_halos[i].id;
      DATA(array,pc,1) = (double)rs_halos[i].num_p;
      DATA(array,pc,2) = rs_halos[i].m;
      DATA(array,pc,3) = rs_halos[i].r;
      DATA(array,pc,4) = rs_halos[i].pos[0];
      DATA(array,pc,5) = rs_halos[i].pos[1];
      DATA(array,pc,6) = rs_halos[i].pos[2];
      DATA(array,pc,7) = rs_halos[i].J[0];
      DATA(array,pc,8) = rs_halos[i].J[1];
      DATA(array,pc,9) = rs_halos[i].J[2];
      pc++;
    }
  }
  else if(values==1){
    for(i=0; i<n_particles; i++){
      DMDATA(array,pc,0) = pids[i];
      DMDATA(array,pc,1) = dmlist[i];
      pc++;
    }
  }
  return;
}
