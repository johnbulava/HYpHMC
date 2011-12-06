#include "ChebyshevOneOverXPlusOneRelError.h"
#include <stdio.h>
#include <math.h>
#include "Tools.h"

void ChebyshevOneOverXPlusOneRelError::ini(int setn, double setEps) {
  if (setn <= 0) setn = 0;
  if (setEps<=0) setEps = 1E-6;
  n = setn;
  epsilon = setEps;
  coeff = new double[n+1];
  calcCoeff();
}


void ChebyshevOneOverXPlusOneRelError::calcCoeff() {
  printf("ChebyshevOneOverXPlusOneRelError: Calcing Coefficients for n = %d...\n",n);

  double t1 = 1;
  double t2 = -1 - epsilon;
  double rho = 1;
  double Tnp1Dummy = 0;
  int I;
  for (I=1; I<=n+1; I++) {
    rho = t2;
    Tnp1Dummy = 2*(-1 - epsilon)*t2 - t1;
    t1 = t2;
    t2 = Tnp1Dummy;
  }
  rho = -1/rho;
  
  double* zeroes = new double[n+1];
  double* target = new double[n+1];
  double* Tnm1 = new double[n+1];
  double* Tn = new double[n+1];

  int I2;
  for (I=0; I<n+1; I++) {
    zeroes[I] = -cos(M_PI*(I+0.5)/(n+1));
    Tnm1[I] = 1;
    Tn[I] = zeroes[I];
  }

  for (I2=0; I2<n; I2++) {
    for (I=0; I<n+1; I++) {
      Tnp1Dummy = 2*zeroes[I]*Tn[I] - Tnm1[I];
      Tnm1[I] = Tn[I];
      Tn[I] = Tnp1Dummy;
    }
  }

  for (I=0; I<n+1; I++) {
    target[I] = (1+rho*Tn[I]) / (1+epsilon+zeroes[I]);
  }
  
  for (I=0; I<n+1; I++) {
    Tnm1[I] = 1;
    Tn[I] = zeroes[I];
  }
  
  coeff[0] = 0;
  for (I=0; I<n+1; I++) {
    coeff[0] += Tnm1[I] * target[I];
  }
  coeff[0] /= (n+1);
    
  int cn = 1;
  while (cn<=n) {    
    coeff[cn] = 0;
    for (I=0; I<n+1; I++) {
      coeff[cn] += Tn[I] * target[I];
	  
      Tnp1Dummy = 2*zeroes[I]*Tn[I] - Tnm1[I];
      Tnm1[I] = Tn[I];
      Tn[I] = Tnp1Dummy;
    }
    coeff[cn] /= 0.5*(n+1);
    cn++;
  }
  
  delete[] zeroes;
  delete[] target;
  delete[] Tnm1;
  delete[] Tn;

  printf("sucessfully.\n");
}


void ChebyshevOneOverXPlusOneRelError::desini() {
  n = 0;
  delete [] coeff;
  coeff = NULL;
}
 

ChebyshevOneOverXPlusOneRelError::ChebyshevOneOverXPlusOneRelError() {
  ini(1,1E-6);
}

 
ChebyshevOneOverXPlusOneRelError::ChebyshevOneOverXPlusOneRelError(int setn, double setEps) {
  ini(setn, setEps);
}


ChebyshevOneOverXPlusOneRelError::~ChebyshevOneOverXPlusOneRelError() {
  desini();
}


void ChebyshevOneOverXPlusOneRelError::recalcCoeff(int setn, double setEps) {
  desini();
  ini(setn, setEps);
}


void ChebyshevOneOverXPlusOneRelError::printToFile() {
  printf("Printing Chebyshev-coefficients to disk...");
  FILE* file = fopen("data/ChebyshevOneOverXPlusOneRelErrorCoefficients.dat","w");
  int I;
  for (I=0; I<=n; I++) {
    fprintf(file,"%d %1.15f\n",I,coeff[I]);
  }
  fclose(file);
  printf("ready.\n");
}


double ChebyshevOneOverXPlusOneRelError::calcControlPlot() {
  double maxRelDiff = 0;
  int p = n+1;
  
  int interP = 100;
  printf("Calcing Control plot for Chebyshey-Polynomial...");
  double* zeroes = new double[interP*p];
  double* Tnm1 = new double[interP*p];
  double* Tn = new double[interP*p];
  double* val = new double[interP*p];
  double Tnp1Dummy = 0;
    
  int I;
  for (I=0; I<interP*p; I++) {
    zeroes[I] = -1.0 + (2.0*I)/(interP*p);
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

  FILE* file = fopen("data/ChebyshevOneOverXPlusOneRelErrorControlPlot.dat","w");
  for (I=0; I<interP*p; I++) {
    if (abs((val[I]-(1/(1+epsilon+zeroes[I]))) * (1+epsilon+zeroes[I])) > maxRelDiff) {
      maxRelDiff = abs((val[I]-(1/(1+epsilon+zeroes[I]))) * (1+epsilon+zeroes[I]));
    }
    fprintf(file,"%1.15f %1.15f %1.15f %1.15f %1.15f\n",zeroes[I],val[I],1/(1+epsilon+zeroes[I]),(val[I]-1/(1+epsilon+zeroes[I])), (val[I]-1/(1+epsilon+zeroes[I]))*(1+epsilon+zeroes[I]));    
  }
  fclose(file);
  printf("ready\n");
  
  delete[] zeroes;
  delete[] val;
  delete[] Tnm1;
  delete[] Tn;
  return maxRelDiff;
}
