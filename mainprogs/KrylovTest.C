#define KRYLOVTESTMODUS
double KRYLOVERRROR = -1;

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
#include "Tools.h"
#include "NeubergerMatrix.h"
#include "WilsonMatrix.h"
#include "FermionMatrixOperations.h"
#include "ChebyshevOneOverXPlusOne.h"
#include "PowerPolynomialsOneOverX.h"
#include "ChebyshevOneOverXPlusOneRelError.h"
#include "SSEroutines.h"
#include "HMCPropagator.h"
#include "ExtremeFFT4D.h"
#include "PolynomialApproximation.h"
#include "pHMCForce.h"
#include "pHMCPropagator.h"
#include "GeneralChebyshevApproximation.h"
#include "MultiThreadedOperations.h"
#include "xSSE.h"


FermionMatrixOperations* fOps = NULL;
HMCPropagator* HMCProp = NULL;
int L0,L1,L2,L3;
Complex* input;
Complex* output1;
Complex* output2;
Complex* output3;
Complex* interim;
Complex* interim2;
double* phiField;
double* dSdPhi;




void randomPhi(int mode) {
  if (mode==0) {
    int I;
    for (I=0; I<4*L0*L1*L2*L3; I++) {
      double z1,z2;
      AdvancedGaussZufall(AdvancedSeed, z1, z2);
      phiField[I] = z1;
    }
  }
}



void randomInput() {
  for (int I=0; I<fOps->getVectorLengthXtrSize(); I++) {
    input[I].x = NaN;
    input[I].y = NaN;    
  }
  for (int I=0; I<fOps->getVectorLength(); I++) {
    double z1,z2;
    AdvancedGaussZufall(AdvancedSeed, z1, z2);
    input[I] = Complex(z1, z2);
  }
  fOps->transformToXtraSizeArray(input, input);
}


double matrixFunction1(double x) {
  return exp(-log(x)/1.0);
}


double matrixFunction2(double x) {
  return exp(-log(x)/2.0);
}


double matrixFunction3(double x) {
  return exp(-log(x)/3.0);
}





void testKrylovMethod(double yN, double TOL) {
  printf("\nTesting MMdag-SQRT-Solver (Neuberger with Xi) with yN = %1.2f, TOL = %1.15f ...\n",yN,TOL);
  
  char* fileName = new char[1000];
  snprintf(fileName,1000,"KrylovTest_L%dx%dx%dx%dy%1.3f.dat",L0,L1,L2,L3,yN);
  FILE* file = fopen(fileName, "w");
  fprintf(file, "#Krylov approach test...\n");
  fclose(file);
  
  
  ComplexVector inputVec(fOps->getVectorLength(),input);
  ComplexVector interimVec(fOps->getVectorLength(),interim);
  
  fOps->setYukawaCoupling(yN);
  fOps->setPreconditioner(false, 1, 0);
  fOps->setQPreconditioner(false, 0.25, 0.25);
  fOps->setRPreconditioner(false, 1, 0);
  
  int maxIter = 1;
  while (true) {
    int neededIterDummy = 0;  
    int neededIter1 = 0;
    int neededIter2 = 0;
    int neededIter3 = 0;

    fOps->applyFermionMatrixMMDaggerFunction(input, output1, phiField, &matrixFunction1, TOL, maxIter, neededIter1, 0, NULL, false, false);
    double errorEst1 = KRYLOVERRROR;

    fOps->applyFermionMatrixMMDaggerFunction(input, interim, phiField, &matrixFunction2, TOL, maxIter, neededIter2, 0, NULL, false, false);
    double errorEst2 = KRYLOVERRROR;
    fOps->applyFermionMatrixMMDaggerFunction(interim, output2, phiField, &matrixFunction2, TOL, 0, neededIterDummy, 0, NULL, false, false);    

    fOps->applyFermionMatrixMMDaggerFunction(input, interim, phiField, &matrixFunction3, TOL, maxIter, neededIter3, 0, NULL, false, false);
    double errorEst3 = KRYLOVERRROR;
    fOps->applyFermionMatrixMMDaggerFunction(interim, interim2, phiField, &matrixFunction3, TOL, 0, neededIterDummy, 0, NULL, false, false);
    fOps->applyFermionMatrixMMDaggerFunction(interim2, output3, phiField, &matrixFunction3, TOL, 0, neededIterDummy, 0, NULL, false, false);
    
    fOps->transformFromXtraSizeArray(input, input);

    fOps->executeFermionMatrixFermionDaggerMatrixMultiplication(output1, interim, phiField, 0, 0, false);
    fOps->transformFromXtraSizeArray(interim, interim);
    ComplexVector diffVec1 = interimVec - inputVec;
    
    fOps->executeFermionMatrixFermionDaggerMatrixMultiplication(output2, interim, phiField, 0, 0, false);
    fOps->transformFromXtraSizeArray(interim, interim);
    ComplexVector diffVec2 = interimVec - inputVec;
    
    fOps->executeFermionMatrixFermionDaggerMatrixMultiplication(output3, interim, phiField, 0, 0, false);
    fOps->transformFromXtraSizeArray(interim, interim);
    ComplexVector diffVec3 = interimVec - inputVec;
    

    double normInput = inputVec.getNorm();
    double normDiffVec1 = diffVec1.getNorm();
    double normDiffVec2 = diffVec2.getNorm();
    double normDiffVec3 = diffVec3.getNorm();
    double relDiffVec1 = normDiffVec1 / normInput;
    double relDiffVec2 = normDiffVec2 / normInput;
    double relDiffVec3 = normDiffVec3 / normInput;

    fOps->transformToXtraSizeArray(input, input);
        
    printf("%d %e %e %e %e %e %e\n", maxIter, relDiffVec1,relDiffVec2,relDiffVec3, errorEst1, errorEst2, errorEst3); 
    file = fopen(fileName, "a");
    fprintf(file, "%d %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e\n", maxIter, relDiffVec1,relDiffVec2,relDiffVec3, errorEst1, errorEst2, errorEst3); 
    fclose(file);
    if ((neededIter1<maxIter) && (neededIter2<maxIter) && (neededIter3<maxIter)) break;
    maxIter++;
  }
  delete[] fileName;
}



void testKrylovMethod2(double yN, double TOL, double fac) {
  printf("\nTesting MMdag-SQRT-Solver (Neuberger with Xi) with yN = %1.2f, TOL = %1.15f and fac = %f...\n",yN,TOL, fac);
  
  char* fileName = new char[1000];
  snprintf(fileName,1000,"KrylovTestTolScan_L%dx%dx%dx%dy%1.3f.dat",L0,L1,L2,L3,yN);
  FILE* file = fopen(fileName, "w");
  fprintf(file, "#Krylov approach test...\n");
  fclose(file);
  
  
  ComplexVector inputVec(fOps->getVectorLength(),input);
  ComplexVector interimVec(fOps->getVectorLength(),interim);
  
  fOps->setYukawaCoupling(yN);
  fOps->setPreconditioner(false, 1, 0);
  fOps->setQPreconditioner(false, 0.25, 0.25);
  fOps->setRPreconditioner(false, 1, 0);
  
  double tol2 = 1.0;
  while (true) {
    int neededIterDummy = 0;  
    int neededIter1 = 0;
    int neededIter2 = 0;
    int neededIter3 = 0;

    fOps->applyFermionMatrixMMDaggerFunction(input, output1, phiField, &matrixFunction1, tol2, 0, neededIter1, 0, NULL, false, false);
    double errorEst1 = KRYLOVERRROR;

    fOps->applyFermionMatrixMMDaggerFunction(input, interim, phiField, &matrixFunction2, tol2, 0, neededIter2, 0, NULL, false, false);
    double errorEst2 = KRYLOVERRROR;
    fOps->applyFermionMatrixMMDaggerFunction(interim, output2, phiField, &matrixFunction2, TOL, 0, neededIterDummy, 0, NULL, false, false);    

    fOps->applyFermionMatrixMMDaggerFunction(input, interim, phiField, &matrixFunction3, tol2, 0, neededIter3, 0, NULL, false, false);
    double errorEst3 = KRYLOVERRROR;
    fOps->applyFermionMatrixMMDaggerFunction(interim, interim2, phiField, &matrixFunction3, TOL, 0, neededIterDummy, 0, NULL, false, false);
    fOps->applyFermionMatrixMMDaggerFunction(interim2, output3, phiField, &matrixFunction3, TOL, 0, neededIterDummy, 0, NULL, false, false);
    
    fOps->transformFromXtraSizeArray(input, input);

    fOps->executeFermionMatrixFermionDaggerMatrixMultiplication(output1, interim, phiField, 0, 0, false);
    fOps->transformFromXtraSizeArray(interim, interim);
    ComplexVector diffVec1 = interimVec - inputVec;
    
    fOps->executeFermionMatrixFermionDaggerMatrixMultiplication(output2, interim, phiField, 0, 0, false);
    fOps->transformFromXtraSizeArray(interim, interim);
    ComplexVector diffVec2 = interimVec - inputVec;
    
    fOps->executeFermionMatrixFermionDaggerMatrixMultiplication(output3, interim, phiField, 0, 0, false);
    fOps->transformFromXtraSizeArray(interim, interim);
    ComplexVector diffVec3 = interimVec - inputVec;
    

    double normInput = inputVec.getNorm();
    double normDiffVec1 = diffVec1.getNorm();
    double normDiffVec2 = diffVec2.getNorm();
    double normDiffVec3 = diffVec3.getNorm();
    double relDiffVec1 = normDiffVec1 / normInput;
    double relDiffVec2 = normDiffVec2 / normInput;
    double relDiffVec3 = normDiffVec3 / normInput;

    fOps->transformToXtraSizeArray(input, input);
        
    printf("%1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %d %d %d %1.15e\n", tol2, relDiffVec1,relDiffVec2,relDiffVec3, errorEst1, errorEst2, errorEst3, neededIter1, neededIter2, neededIter3, TOL); 
    file = fopen(fileName, "a");
    fprintf(file, "%1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %d %d %d %1.15e\n", tol2, relDiffVec1,relDiffVec2,relDiffVec3, errorEst1, errorEst2, errorEst3, neededIter1, neededIter2, neededIter3, TOL); 
    fclose(file);
    tol2 *= fac;
    if (tol2<=TOL) break;
  }
  delete[] fileName;
}



int main(int argc,char **argv) {
  LogLevel = 3;

  int FFTW_ThreadCount = 2;

  if (LogLevel>1) printf("Initialisiere FFTW with %d threads.\n",FFTW_ThreadCount);
  fftw_init_threads();
  fftw_plan_with_nthreads(FFTW_ThreadCount);  
  
  iniTools(1);
  L0 = 4;
  L1 = 4;
  L2 = 4;
  L3 = 8;
  bool usexFFT = false;
  int threadCountPerNode = 2;
  int ParaOpMode = 1;
  char* fftPlanDescriptor = NULL;
  readOptimalFermionVectorEmbeddingAndFFTPlanFromTuningDB(L0, L1, L2, L3, threadCountPerNode, ParaOpMode, usexFFT, 1, 1, 1, 1, fftPlanDescriptor);
  delete[] fftPlanDescriptor;
    
  fOps = new FermionMatrixOperations(L0, L1, L2, L3, 1.0, 0.5, 1.0);
  fOps->setxFFTusage(usexFFT);  
  input = fOps->createFermionVector();
  output1 = fOps->createFermionVector();
  output2 = fOps->createFermionVector(); 
  output3 = fOps->createFermionVector(); 
  interim = fOps->createFermionVector(); 
  interim2 = fOps->createFermionVector(); 

  
  HMCProp = new HMCPropagator(fOps, 0.01, 0.10, 10, 2.1);
  phiField = (double*) HMCProp->phiField;
  dSdPhi = (double*) fOps->createFermionVector(2);
  
  initializePerformanceProfiler("testOpsPerformanceProfile.dat");

  randomPhi(0);
  randomInput();

  testKrylovMethod2(1.0, 1E-10, 0.85);


  testKrylovMethod(1.0, 1E-10);


 
 
 
 
  fOps->destroyFermionVector(input);
  fOps->destroyFermionVector(output1);
  fOps->destroyFermionVector(output2);
  fOps->destroyFermionVector(output3);
  fOps->destroyFermionVector(interim);
  fOps->destroyFermionVector(interim2);
  Complex* c = (Complex*)dSdPhi;
  fOps->destroyFermionVector(c);
  dSdPhi = NULL;
  
  delete fOps;
  delete HMCProp;
  writePerformanceProfilingDataToDisk();
  desiniTools();
  fftw_cleanup_threads();
  if (LogLevel>1) printf("KrylovTest terminated correctly.\n");
}
