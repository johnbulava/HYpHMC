#include "PowerPolynomialsOneOverX.h"

void PowerPolynomialsOneOverX::ini(int setN, double setEps, double setAl, double setBet) {
  if (setN <= 0) setN = 1;
  if (setEps<=0) setEps = 1E-6;
/*  if (setAl<=0) setAl = 0;
  if (setAl>=2) setAl = 1.99;
  if (setBet<=0) setBet = 0;
  if (setBet>=1) setBet = 1.0;*/
  
  if ((2*setBet-setAl) <=-1) {
    printf("ERROR!!! Wrong parameters...\n");
    exit(0);
  }
  
  N = setN;
  epsilon = setEps;
  alpha = setAl;
  beta = setBet;
  
  Pcoeff = new long double*[N+1];
  int I;
  for (I=0; I<=N; I++) {
    Pcoeff[I] = new long double[N+1];
  }  
  AppCoeff = new long double[N+1];
  calcPCoeff();
  calcAppCoeff();  
}


long double PowerPolynomialsOneOverX::PScalar(long double* p1, long double* p2) {
  int n,m;
  long double s = 0.0;
  
  for (n=0; n<=N; n++) {
    for (m=0; m<=N; m++) {
      long double f = epsilon;
      f = expl(logl(f)*(1+2*beta-alpha+n+m));
      f = (1 - f) / (1+2*beta-alpha+n+m);
      
      s = s + f * p1[n] * p2[m];
    }
  }

  return s;
}


long double PowerPolynomialsOneOverX::PScalarWithOneOverX(long double* p) {
  int n;
  long double s = 0.0;
  
  for (n=0; n<=N; n++) {
    long double f = epsilon;
    f = expl(logl(f)*(1+2*beta-alpha+n-beta));
    f = (1 - f) / (1+2*beta-alpha+n-beta);
    s = s + f * p[n] * 1;
  }

  return s;
}


void PowerPolynomialsOneOverX::calcPCoeff() {
  printf("PowerPolynomialsOneOverX: Calcing Coefficients of PowerPolynomials for N = %d...\n",N);

  int n,m,k;
  long double s;
  for (n=0; n<=N; n++) {
    for (m=0; m<=N; m++) {
      Pcoeff[n][m] = 0;
    }
  }
  
  for (n=0; n<=N; n++) {
    Pcoeff[n][n] = 1.0;
    
    for (m=0; m<n; m++) {
      s = PScalar(Pcoeff[n],Pcoeff[m]);
      for (k=0; k<=m; k++) {
        Pcoeff[n][k] = Pcoeff[n][k] - s * Pcoeff[m][k];
      }
    }
    
    s = expl(0.5*logl(PScalar(Pcoeff[n], Pcoeff[n])));
    for (k=0; k<=n; k++) {
      Pcoeff[n][k] = Pcoeff[n][k] / s;
    }
  }

  printf("sucessfully.\n");
}


void PowerPolynomialsOneOverX::calcAppCoeff() {
  printf("PowerPolynomialsOneOverX: Calcing Coefficients of Approximating Polynomial of 1/x for N = %d...\n",N);
  int n,m;
  long double s;
  
  for (n=0; n<=N; n++) {
    AppCoeff[n] = 0;
  }  
  
  for (n=0; n<=N; n++) {
    s = PScalarWithOneOverX(Pcoeff[n]);
    for (m=0; m<=n; m++) {
      AppCoeff[m] = AppCoeff[m] + s*Pcoeff[n][m];
    }
  }

  printf("sucessfully.\n");
}


void PowerPolynomialsOneOverX::desini() {
  N = 0;
  int I;
  for (I=0; I<=N; I++) {
    delete [] Pcoeff[I];
  }
  delete [] Pcoeff;  
  Pcoeff = NULL;
  delete [] AppCoeff;
  AppCoeff = NULL;
}
 

PowerPolynomialsOneOverX::PowerPolynomialsOneOverX() {
  ini(1,1,1E-4,1);
}

 
PowerPolynomialsOneOverX::PowerPolynomialsOneOverX(int setN, double setEps, double setAl, double setBet) {
  ini(setN, setEps, setAl, setBet);
}


PowerPolynomialsOneOverX::~PowerPolynomialsOneOverX() {
  desini();
}


void PowerPolynomialsOneOverX::recalcCoeff(int setN, double setEps, double setAl, double setBet) {
  desini();
  ini(setN, setEps, setAl, setBet);
}


void PowerPolynomialsOneOverX::printToFile() {
  printf("Printing PowerPolynomial-coefficients to disk...");
  FILE* file = fopen("data/PowerPolynomialsOneOverXAppCoefficients.dat","w");
  int n;
  for (n=0; n<=N; n++) {
    fprintf(file,"%d %1.15f\n", n, (double) (AppCoeff[n]));
  }
  fclose(file);
  file = fopen("data/PowerPolynomialsOneOverXPCoefficients.dat","w");
  int m;
  for (n=0; n<=N; n++) {
    for (m=0; m<=N; m++) {
      fprintf(file,"%1.15f ", (double) (Pcoeff[n][m]));
    }
    fprintf(file,"\n");
  }
  fclose(file);
  printf("ready.\n");
}


double PowerPolynomialsOneOverX::calcControlPlot() {
  double maxRelDiff = 0;
  long double x = epsilon;    
  long double s;
  int n;
  
  FILE* file = fopen("data/PowerPolynomialsOneOverXControlPlot.dat","w");
  while (x<=1.0) {
    s = 0;
    for (n=0; n<=N; n++) {
      s = s + AppCoeff[n] * expl(n*logl(x));
    }

    double v3 = 1;
    double v4 = 0;
    double v5 = 0;
    if (x>=epsilon) { 
      v3 = (double) (expl(-beta*logl(x)));
      v4 = (double) (s-v3);
      v5 = (double) (v4/v3);
      if (fabs(v5) > maxRelDiff) {
        maxRelDiff = fabs(v5);
      }
    }
    
    fprintf(file,"%1.15f %1.15f %1.15f %1.15f %1.15f\n",(double) x,(double) (s), v3, v4, v5);    
    x = x + epsilon/100;
  }
  fclose(file);
  printf("ready\n");

  return maxRelDiff;
}
