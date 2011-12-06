#include "ChebyshevOneOverXPlusOne.h"

void ChebyshevOneOverXPlusOne::ini(int setN, int setn, double setEps) {
  if (setN <= -1) setN = -1;
  if (setN==0) setN = 1;
  if (setn <= 0) setn = 0;
  if ((setN>0) && (setn>=setN)) setn = setN-1;
  if (setEps<=0) setEps = 1E-6;
  N = setN;
  n = setn;
  epsilon = setEps;
  if (N>0) {
    epsilon -= (1-cos(M_PI*0.5/N));
  }
  coeff = new double[n+1];
  calcCoeff();
}


void ChebyshevOneOverXPlusOne::calcCoeff() {
  printf("ChebyshevOneOverXPlusOne: Calcing Coefficients for N = %d and n = %d...\n",N,n);
  if (N>0) {
    double* zeroes = new double[N];
    double* Tnm1 = new double[N];
    double* Tn = new double[N];
    double Tnp1Dummy = 0;

    int I;
    for (I=0; I<N; I++) {
      zeroes[I] = -cos(M_PI*(I+0.5)/N);
      Tnm1[I] = 1;
      Tn[I] = zeroes[I];
    }
    printf("Smallest Zero: %1.15f, Largest Zero %1.15f\n", zeroes[0], zeroes[N-1]);
    coeff[0] = 0;
    for (I=0; I<N; I++) {
      coeff[0] += Tnm1[I]/(1+epsilon+zeroes[I]);
    }
    coeff[0] /= N;
    
    int cn = 1;
    while (cn<=n) {    
      coeff[cn] = 0;
      for (I=0; I<N; I++) {
        coeff[cn] += Tn[I]/(1+epsilon+zeroes[I]);
	  
        Tnp1Dummy = 2*zeroes[I]*Tn[I] - Tnm1[I];
        Tnm1[I] = Tn[I];
        Tn[I] = Tnp1Dummy;
      }
      coeff[cn] /= 0.5*N;
      cn++;
    }
  
    delete[] zeroes;
    delete[] Tnm1;
    delete[] Tn;
  } else {
    coeff[0] = 1/sqrt(epsilon*(2+epsilon));
    if (n>0) {
      coeff[1] = 2 - 2*(epsilon+1)*coeff[0];
    }
    if (n>1) {
      coeff[2] = -2*coeff[0] - 2*(epsilon+1)*coeff[1];
    }
    int cn=3;
    while (cn<=n) {    
      coeff[cn] = -coeff[cn-2] - 2*(epsilon+1)*coeff[cn-1];
      cn++;
    }
  }

  printf("sucessfully.\n");
}


void ChebyshevOneOverXPlusOne::desini() {
  N = 0;
  n = 0;
  delete [] coeff;
  coeff = NULL;
}
 

ChebyshevOneOverXPlusOne::ChebyshevOneOverXPlusOne() {
  ini(1,1,1E-6);
}

 
ChebyshevOneOverXPlusOne::ChebyshevOneOverXPlusOne(int setN, int setn, double setEps) {
  ini(setN, setn, setEps);
}


ChebyshevOneOverXPlusOne::~ChebyshevOneOverXPlusOne() {
  desini();
}


void ChebyshevOneOverXPlusOne::recalcCoeff(int setN, int setn, double setEps) {
  desini();
  ini(setN, setn, setEps);
}


void ChebyshevOneOverXPlusOne::printToFile() {
  printf("Printing Chebyshev-coefficients to disk...");
  FILE* file = fopen("data/ChebyshevOneOverXPlusOneCoefficients.dat","w");
  int I;
  for (I=0; I<=n; I++) {
    fprintf(file,"%d %1.15f\n",I,coeff[I]);
  }
  fclose(file);
  printf("ready.\n");
}


double ChebyshevOneOverXPlusOne::calcControlPlot() {
  double maxRelDiff = 0;
  int p = N;
  double offset = 0.5;
  
  if (N==-1) {
    p = (int)(sqrt(1/epsilon));
    offset = 0;
  }

  int interP = 100;
  printf("Calcing Control plot for Chebyshey-Polynomial...");
  double* zeroes = new double[interP*p];
  double* Tnm1 = new double[interP*p];
  double* Tn = new double[interP*p];
  double* val = new double[interP*p];
  double Tnp1Dummy = 0;
    
  int I;
  for (I=0; I<interP*p; I++) {
    zeroes[I] = -cos(M_PI*(I+interP*offset)/(interP*p));
    val[I] = coeff[0] * 1;
    Tnm1[I] = 1;
    Tn[I] = zeroes[I];
  }

  int cn = 1;
  while (cn<=n) {    
    for (I=0; I<interP*p; I++) {
      val[I] += coeff[cn] * Tn[I];

      Tnp1Dummy = 2*zeroes[I]*Tn[I] - Tnm1[I];
      Tnm1[I] = Tn[I];
      Tn[I] = Tnp1Dummy;
    }
    cn++;
  }

  FILE* file = fopen("data/ChebyshevOneOverXPlusOneControlPlot.dat","w");
  for (I=0; I<interP*p; I++) {
    if (abs((val[I]-(1/(1+epsilon+zeroes[I]))) * (1+epsilon+zeroes[I])) > maxRelDiff) {
      maxRelDiff = abs((val[I]-(1/(1+epsilon+zeroes[I]))) * (1+epsilon+zeroes[I]));
    }
    fprintf(file,"%1.15f %1.15f %1.15f %1.15f %1.15f %1.15f\n",zeroes[I],val[I],1/(1+epsilon+zeroes[I]),(val[I]-1/(1+epsilon+zeroes[I])), (val[I]-1/(1+epsilon+zeroes[I]))*(1+epsilon+zeroes[I]), (val[I]-1/(1+epsilon+zeroes[I]))*sqr(1+epsilon+zeroes[I]));    
  }
  fclose(file);
  printf("ready\n");
  
  delete[] zeroes;
  delete[] val;
  delete[] Tnm1;
  delete[] Tn;
  return maxRelDiff;
}
