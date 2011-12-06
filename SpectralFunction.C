#include "EvaluateObservablePropagatorBase.h"
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
#include "Tools.h"
#include "FermionMatrixOperations.h"
#include "ComplexPolynom.h"
#include "HighPrecisionComplex.h"
#include "HighPrecisionComplexPolynom.h"
#include <cln/number.h>
#include <cln/float.h>
#include <cln/complex.h>



int main(int argc,char **argv) {
  iniTools(5517);

  double mH = 0.2;
  double mG=0.01;
  double lamRen=0.01;
  int n=1;

  double vren=sqrt((mH*mH - mG*mG)/(8*lamRen));
  Complex HiggsPole = findPoleOfBosonic1LoopPropagatorFromRenPT(mH, mG, lamRen, vren, n);
  double Gamma = 2*HiggsPole.x;
  printf("Vren: %1.4f, Width: %1.4e\n", vren, Gamma);

  char* fileName = new char[1000];
  snprintf(fileName, 1000, "SpectralFunctionMh%1.5fMg%1.5fLamRen%1.5fVRen%1.5fGamma%1.5fN%d.dat",mH,mG,lamRen,vren,Gamma,n);
  FILE* fileS = fopen(fileName, "w");
  delete[] fileName;
  int N=10000;
  for (int I=0; I<N; I++) {
    double E = 0.75*I/N;
    double RHO = calcSpectralFunctionOfBosonic1LoopPropagatorFromRenPT(E, HiggsPole, mG, lamRen, vren, n);
    fprintf(fileS, "%1.15f %1.15f\n", E, RHO);
  }
  fclose(fileS);
}

