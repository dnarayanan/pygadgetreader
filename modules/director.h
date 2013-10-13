#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

void readpos(){
  if(Tipsy)
    tipsy_posvel();
  else
    gadget_posvel();
}
void readvel(){
  if(Tipsy)
    tipsy_posvel();
  else
    gadget_posvel();
}

void readmass(){
  if(Tipsy)
    tipsy_bin();
  else
    gadget_mass();
  return;
}

void readu(){
  if(Tipsy)
    tipsy_bin();
  else
    gas_props();
  return;
}

void readrho(){
  if(Tipsy)
    tipsy_bin();
  else
    gas_props();
 return;
}

void readNE(){
  if(Tipsy)
    tipsy_aux();
  else
    gas_props();
  return;
}

void readNH(){
  if(Tipsy)
    tipsy_aux();
  else
    gas_props();
  return;
}

void readHSML(){
  if(Tipsy)
    tipsy_bin();
  else
    gas_props();
  return;
}

void readSFR(){
  if(Tipsy)
    tipsy_aux();
  else
    gas_props();
  return;
}

void readdelaytime(){
  if(Tipsy)
    tipsy_aux();
  else
    gas_props();
  return;
}

void readfh2(){
  if(Tipsy==0)
    gas_props();
  return;
}
void readsigma(){
  if(Tipsy==0)
    gas_props();
  return;
}

void readage(){
  if(Tipsy)
    tipsy_bin();
  else
    gadget_readage();
  return;
}

void readZ(){
  if(Tipsy)
    tipsy_bin();
  else
    gadget_readZ();
  return;
}

void readmetals(){
  if(Tipsy)
    tipsy_auxmetals();
  else
    gadget_readmetals();
  return;
}

void readtmax(){
  if(Tipsy)
    tipsy_aux();
  else
    gadget_readtmax();
  return;
}

void readnspawn(){
  if(Tipsy)
    tipsy_aux();
  else
    gadget_readnspawn();
  return;
}

void readpotential(){
  if(Tipsy)
    tipsy_bin();
  else
    gadget_readpotential();
  return;
}

void readpid(){
  if(Tipsy)
    read_tipsy_ID();
  else
    gadget_readpid();
}
