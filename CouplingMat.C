#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <iostream>


#include "Complex.h"
#include "Quat.h"
#include "Global.h"
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include "Tools.C"


int N, Nhalf;
double rho = 1.0;
double r = 0.5;
double* summands = NULL;




int main(int argc,char **argv) {
  iniTools(1);
  N = 120;
  Nhalf = N / 2;
  LogLevel = 2;

  printf("Calculating eigenvalues alpha...\n");
  summands = new double[N*N*N*N];
  int I;
  for (I=0; I<N*N*N*N; I++) summands[I] = NaN;
  double fac = (2*pi)/N;
  vector4D p;
  vector4D sinp;
  vector4D sinphalf;
  double pTildeSquare = 0;
  double pDoubleTildeSquare = 0;
  Complex result;
  double nenner;
  int I0,I1,I2,I3;
  int pos = 0;
  for (I0=0; I0<N; I0++) {
    p[0] = fac * I0;
    sinp[0] = sin(p[0]);
    sinphalf[0] = 2*sin(p[0]/2);
    for (I1=0; I1<N; I1++) {
      p[1] = fac * I1;
      sinp[1] = sin(p[1]);
      sinphalf[1] = 2*sin(p[1]/2);
      for (I2=I1; I2<N; I2++) {
        p[2] = fac * I2;
        sinp[2] = sin(p[2]);
        sinphalf[2] = 2*sin(p[2]/2);
        for (I3=I2; I3<N; I3++) {
          p[3] = fac * I3;
          sinp[3] = sin(p[3]);
          sinphalf[3] = 2*sin(p[3]/2);
	  pos = I0*N*N*N + I1*N*N + I2*N + I3;
	  pTildeSquare = sinp[0]*sinp[0] + sinp[1]*sinp[1] + sinp[2]*sinp[2] + sinp[3]*sinp[3];
	  pDoubleTildeSquare = sinphalf[0]*sinphalf[0] +sinphalf[1]*sinphalf[1] +sinphalf[2]*sinphalf[2] +sinphalf[3]*sinphalf[3];
	  
          nenner = sqrt(pTildeSquare + sqr(r*pDoubleTildeSquare - rho));
          result.x = rho + rho*(r*pDoubleTildeSquare - rho)/nenner;
          result.y = rho*sqrt(pTildeSquare)/nenner;
	  
	  summands[pos] = sinp[0] * sqrt((result.x*result.x + result.y*result.y) / (((result.x-2*rho)*(result.x-2*rho) + result.y*result.y) * pTildeSquare));
	  
	  if (I2==1) {
//	    result.print();
	//    printf("%d\n",pos);
	  
	  }
	  
	  
	  if (((I0==0)||(I0==Nhalf)) &&((I1==0)||(I1==Nhalf)) &&((I2==0)||(I2==Nhalf)) &&((I3==0)||(I3==Nhalf))) {
	    summands[pos] = 0;
	  }

	}
      }
//      exit(0);
    }
  }



  
  printf("Copying eigenvalues alpha...\n");
  for (I0=0; I0<N; I0++) {
    for (I1=0; I1<N; I1++) {
      for (I2=I1; I2<N; I2++) {
        for (I3=I2; I3<N; I3++) {	
	  int p1 = I0*N*N*N + I1*N*N + I2*N + I3;
          int p2;
	  
	  if (I0>Nhalf) {
  	    p2 = I0*N*N*N + I1*N*N + I2*N + I3;
	    summands[p2] = -summands[p1];
  	    p2 = I0*N*N*N + I1*N*N + I3*N + I2;
	    summands[p2] = -summands[p1];
  	    p2 = I0*N*N*N + I2*N*N + I1*N + I3;
	    summands[p2] = -summands[p1];
  	    p2 = I0*N*N*N + I2*N*N + I3*N + I1;
	    summands[p2] = -summands[p1];
  	    p2 = I0*N*N*N + I3*N*N + I1*N + I2;
	    summands[p2] = -summands[p1];
  	    p2 = I0*N*N*N + I3*N*N + I2*N + I1;
	    summands[p2] = -summands[p1];
	  } else {
  	    p2 = I0*N*N*N + I1*N*N + I2*N + I3;
	    summands[p2] = summands[p1];
  	    p2 = I0*N*N*N + I1*N*N + I3*N + I2;
	    summands[p2] = summands[p1];
  	    p2 = I0*N*N*N + I2*N*N + I1*N + I3;
	    summands[p2] = summands[p1];
  	    p2 = I0*N*N*N + I2*N*N + I3*N + I1;
	    summands[p2] = summands[p1];
  	    p2 = I0*N*N*N + I3*N*N + I1*N + I2;
	    summands[p2] = summands[p1];
  	    p2 = I0*N*N*N + I3*N*N + I2*N + I1;
	    summands[p2] = summands[p1];
	  }
	}
      }
    }
  }

  for (I=0; I<N*N*N*N; I++) { 
    if (isNaN(summands[I])) {
      printf("Nan at %d\n",I);
      exit(0);
    }
    
  }


  char* fileName = new char[200];
  snprintf(fileName,200,"data/CouplingMatL%dg.dat",N);
  FILE* file = fopen(fileName,"w");
  int dx;
  for (dx=2; dx<Nhalf; dx+=2) {
    printf("Calculting Fourier-Trafo for x = %d...\n",dx);
    
    double res = 0;
    for (I0=0; I0<N; I0++) {
      p[0] = fac * I0;
      sinp[0] = sin(p[0]*dx);    
      for (I1=0; I1<N; I1++) {
        for (I2=0; I2<N; I2++) {
          for (I3=0; I3<N; I3++) {	
  	    pos = I0*N*N*N + I1*N*N + I2*N + I3;
     
            res += sinp[0] * summands[pos];
	  }
	}
      }
    }
    res /= N*N*N*N;
    printf("...result is: %1.5f\n",res);
    fprintf(file,"%d %1.15f\n",dx,sqrt(sqr(res)));
    
  }
  fclose(file);
  
  
  
  
  
  
  
  
  
}
