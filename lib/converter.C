#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include <stddef.h>
#include <sys/types.h>
#include <dirent.h>

  
double zero = 0.0;  
double NaN = 0.0/zero;
  

void convert(char* fileName) {
  double dummy_phi, dummy_stagphi;
  double dummy_det;
  double dummy_logDet_1, dummy_logDet_2, dummy_logDet_3,dummy_logDet_4;
  double dummy_logDet_5, dummy_logDet_6, dummy_logDet_7,dummy_logDet_8;
  char* restLine = new char[100000];



  printf("Reading from data file: '%s'\n",fileName);	    
  FILE* infile = fopen(fileName,"r");
  char* outFileName = new char[500];
  snprintf(outFileName,500,"converted/%s",fileName);
  printf("Writing to data file: '%s'\n",outFileName);	    
  
  
  FILE* outfile = fopen(outFileName,"w");

  while (fscanf(infile,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &dummy_phi, &dummy_stagphi,
   &dummy_logDet_1,&dummy_logDet_2,&dummy_logDet_3,&dummy_logDet_4,
   &dummy_logDet_5,&dummy_logDet_6,&dummy_logDet_7,&dummy_logDet_8)==10) {
    fgets(restLine, 100000, infile);
    
    fprintf(outfile,"%f %f %f %f %f %f %f %f %f %f\n",dummy_phi, dummy_stagphi, dummy_logDet_1, dummy_logDet_2, dummy_logDet_3, dummy_logDet_4, 
    dummy_logDet_5, dummy_logDet_6, dummy_logDet_7, dummy_logDet_8);

  }

  fclose(infile);
  fclose(outfile);
}


void scanForFiles() {
  DIR *dp;
  struct dirent *ep;
  double dummy_kappa, dummy_lambda, dummy_Y, dummy_rho, dummy_r;
  int dummy_N, dummy_Nf;
     
  dp = opendir (".");
  if (dp != NULL) {
    while (ep = readdir (dp)) {
      if (sscanf(ep->d_name,"measureN%dNf%dKap%lfLam%lfY%lfRho%lfR%lf.dat",&(dummy_N),&(dummy_Nf),&(dummy_kappa),&(dummy_lambda),&(dummy_Y),&(dummy_rho),&(dummy_r))==7) {
        printf("Found data file: '%s'\n",ep->d_name);	    

        convert(ep->d_name);
      }
    }
    closedir (dp);
  } else {
    printf("Couldn't open the '.' directory!!!\n");
    exit(0);
  }
}



int main(int argc,char **argv) {
  printf("   ***   Converter started!!!   ***\n\n");
  printf("NaN value = %f\n",NaN);
  scanForFiles();
}
