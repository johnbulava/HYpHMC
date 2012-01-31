#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include <pthread.h>

#include "Complex.h"
#include "Quat.h"
#include "Global.h"
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include "Tools.C"



double arithgeo(double z, double a, double b, double& agm) {
  if (a <= b) {
    agm = a;
    return sin( z * a);
  }
  double xi = arithgeo(z, 0.5*(a+b), sqrt(a*b), agm);
  return 2*a*xi / ((a+b) + (a-b)*xi*xi);
}


void sncndnK(double z, double k, double& sn, double& K) {
  double agm = 0;
  sn = arithgeo(z, 1.0, sqrt(1.0-k*k), agm);
  K = M_PI / (2.0 * agm);
}


double getCl(double k, int n, int l) {
  double K = 0;
  double sn = 0;
  sncndnK(0, k, sn, K);
  
  double arg = l*K / (2*n+1);
  sncndnK(arg, k, sn, K);
  
  return sn*sn / (1 - sn*sn);
}


double evalFunc(double x, double* cl, int n) {
  double val1 = 1.0;
  for (int l=1; l<=n; l++) {
    double c = cl[2*l];  
    val1 *= (x+c) / (1.0 + c);
  }
  double val2 = 1.0;
  for (int l=1; l<=n; l++) {
    double c = cl[2*l-1];
    val2 *= (x+c) / (1.0 + c);
  }
  return val1 / val2;
}


double findNormFac(double* cl, int n, double k) {
  double xmax = 1.0 / (1.0 - k*k);
  double fmax = 0;
  double fmin = 1E10;
  
  for (int I=0; I<10000; I++) {
    double x = 1+ I*(xmax-1)/10000.0;
    double f = sqrt(x)*evalFunc(x, cl, n) - 1.0;
    if (f>fmax) fmax = f;
    if (f<fmin) fmin = f;
  }
  
  return (1.0 + 0.5*(fmax-fmin)) / (1+fmax);
}


int main(int argc,char **argv) {
  iniTools(5517);

  int n = 9;
  double xmax = 1000;
  double k = sqrt(1.0 - 1.0/xmax);
  double* cl = new double[2*n+1];
  cl[0] = 0;
  for (int I=1; I<=2*n; I++) {
    cl[I] = getCl(k, n, I);
  }
  
  double fac = findNormFac(cl, n, k);
  
  
  FILE* file = fopen("ZolotarevApprox.dat","w");
  fprintf(file, "plot [1:%1.15f] sqrt(x)*%1.15f*(", 1.0/ (1.0-k*k), fac);
  for (int l=1; l<=n; l++) {
    double c = cl[2*l];
    fprintf(file, "((x + %1.15f)/(1+%1.15f))",c,c);
    if (l<n) fprintf(file, " * ");
  }
  fprintf(file, ") / (");
  for (int l=1; l<=n; l++) {
    double c = cl[2*l-1];
    fprintf(file, "((x + %1.15f)/(1+%1.15f))",c,c);
    if (l<n) fprintf(file, " * ");
  }
  fprintf(file, ")-1\n");
  fclose(file);
  delete[] cl;

}

