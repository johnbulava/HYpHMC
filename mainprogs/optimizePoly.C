#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include <stddef.h>
#include <sys/types.h>
#include <dirent.h>

#include "PowerPolynomialsOneOverX.h"
#include "ChebyshevOneOverXPlusOneRelError.h"


double beta, epsilon;
int N, invBeta, NegExpEps;
double* coeff;
double bestNorm;


double calcNorm(double prec) {
  double x = -1;
  long double s;
  long double tnm1;
  long double tn;
  long double tnp1;
  double rel;
  double maxRel = 0;
  
  int I;
  while (x<=1.0) {
    s = coeff[0];
    tnm1 = 1;
    tn = x;
    for (I=1; I<=N; I++) {
      s += coeff[I] * tn;
      tnp1 = 2*x*tn - tnm1;
      tnm1 = tn;
      tn = tnp1;
    }
    rel = s*exp(beta*log(1+epsilon+x)) - 1;
    if (fabs(rel)>maxRel) maxRel = fabs(rel);
  
    x += prec;
  }
  return maxRel;
}

/*double calcNorm(double prec) {
  double x = epsilon;
  long double s;
  long double d;
  double rel;
  double maxRel = 0;
  
  int I;
  while (x<=1.0) {
    s = 0;
    d = 1;
    for (I=0; I<=N; I++) {
      s += coeff[I] * d;
      d *= x;
    }
    rel = s*exp(beta*log(x)) - 1;
    if (fabs(rel)>maxRel) maxRel = fabs(rel);
  
    x += prec;
  }
  return maxRel;
}*/
void iniCoeff() {
  coeff = new double[N+1];

  if (invBeta>1) {
    printf("Initializing Coeff with Power Polynomials");
    PowerPolynomialsOneOverX  PowerPolynomialTest(N, epsilon, 0.0, beta);   
  
    int I;
    for (I=0; I<=N; I++) {  
      coeff[I] = (double) PowerPolynomialTest.AppCoeff[I];
    }
  } else {
    printf("Initializing Coeff with Chebyshev Polynomials\n");
    ChebyshevOneOverXPlusOneRelError  ChebyshevRelErrorTest(N, epsilon);   
    int I;
    for (I=0; I<=N; I++) {  
      coeff[I] = (double) ChebyshevRelErrorTest.coeff[I];
    }
  }
  bestNorm = calcNorm(1E-6);
  printf("Start-Norm: %f\n", bestNorm);
}


void simEnealing(int iter, double T) {
  int I,I2;
  double SOld = calcNorm(1E-3);;
  double SNew;
  double pFac = 1;
  for (I=0; I<iter; I++) {
    for (I2=0; I2<=N; I2++) {
      double oldC = coeff[I2];
      coeff[I2] += 0.0002*(AdvancedZufall(AdvancedSeed)-0.5);
      SNew = calcNorm(pFac*1E-3);
//      printf("%f %f\n",SNew,pFac);
      if (SNew<bestNorm) {
        printf("FirstFound\n");
        SNew = calcNorm(1E-6);
        if (SNew<bestNorm) {
          bestNorm = SNew;
	  printf("New best Norm: %f\n", bestNorm);
	  exit(0);
        } else {
	  printf("   ***  Increasing Precision !!!\n");
	  pFac = pFac*0.1;
	}
      }
      
      if (AdvancedZufall(AdvancedSeed) < exp(-(SNew-SOld)/T)) {
        SOld = SNew;
      } else {
        coeff[I2] = oldC;
      }
    }
  }
}




int main(int argc,char **argv) {
  if (argc != 4) {
    printf("Wrong number of parameters!\n");
    exit(0);
  }
  if ((sscanf(argv[1],"%d",&N)!=1) || (sscanf(argv[2],"%d",&invBeta)!=1) || (sscanf(argv[3],"%d",&NegExpEps)!=1)) {
    printf("Parameter error!\n");
    exit(0);
  }
  beta = 1.0/invBeta;
  epsilon = exp(-NegExpEps*log(10));
  printf("\n                           *** Optimizing Polynomial for N=%d, beta=%f, epsilon=%1.15f *** \n\n",N,beta,epsilon);
  
  iniCoeff();  
  
  double T = 1000;
  while (true) {
    printf("Sim Step\n");
    simEnealing(1000,T);
    T = T * 0.5;
  }
  
  
  
  delete[] coeff;




}
