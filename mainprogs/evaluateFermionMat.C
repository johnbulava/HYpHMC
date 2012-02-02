#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <iostream>


#include "Complex.h"



int main(int argc,char **argv) {
  int N = 6;
  int len = 8*N*N*N*N;
  int count = 0;
  Complex avg[len];
  Complex sigma[len];
  int I;
  for (I=0; I<len; I++) {
    avg[I].x=0;
    avg[I].y=0;
    sigma[I].x=0;
    sigma[I].y=0;    
  }
  

  char* inFileName = "dataBase/data/results/hybrid/FermionMatColN6Nf2Kap0.150Lam0.050Y0.100Rho1.000R0.500Job1.dat";
  char* outFileName = "data/FermiTestN6Y0.1.dat";
  
  printf("Reading data...\n");
  FILE* inFile = fopen(inFileName,"r");
  bool ok = true;
  while (ok) {
    double x,y;
    for (I=0; I<len; I++) {
      if (fscanf(inFile,"%lf %lf ",&x,&y)!=2) {
        ok = false;
	break;
      }
      avg[I].x += x;
      avg[I].y += y;
      sigma[I].x += x*x;
      sigma[I].y += y*y;
    }
    if (ok) count++;
  }
  fclose(inFile);
  
  printf("Calculating data, count = %d...\n",count);
  for (I=0; I<len; I++) {
    avg[I].x /= count;
    avg[I].y /= count;
    sigma[I].x /= count;
    sigma[I].y /= count;
    sigma[I].x = sqrt(sigma[I].x - avg[I].x*avg[I].x) / sqrt(count);
    sigma[I].y = sqrt(sigma[I].y - avg[I].y*avg[I].y) / sqrt(count);
  }
  
  
  printf("Writing data...\n");
  FILE* outFile = fopen(outFileName,"w");
  for (I=0; I<len; I++) {
    int d3 = (I/8) % N;
    if (d3>N/2) d3 -= N; 
    int d2 = (I/(8*N)) % N;
    if (d2>N/2) d2 -= N; 
    int d1 = (I/(8*N*N)) % N;
    if (d1>N/2) d1 -= N; 
    int d0 = (I/(8*N*N*N)) % N;
    if (d0>N/2) d0 -= N; 

    double dt = fabs(d0)+fabs(d1)+fabs(d2)+fabs(d3);
    double dd = sqrt(d0*d0+d1*d1+d2*d2+d3*d3);
    
    fprintf(outFile,"%1.15f %1.15f %1.15f %1.15f %1.15f %1.15f\n",dt,dd,avg[I].x,avg[I].y,sigma[I].x,sigma[I].y);
  }  
  
  fclose(outFile);
  printf("READY\n");
  
}
