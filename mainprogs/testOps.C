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





double startTime;
double startCPUTime;
FermionMatrixOperations* fOps = NULL;
HMCPropagator* HMCProp = NULL;
int L0,L1,L2,L3;
Complex* input;
Complex* output;
Complex* output2;
Complex* interim;
Complex* interim2;
double* phiField;
double* dSdPhi;



void startTimer() {
  startTime = zeitwert();
  startCPUTime = cpuTime();
}


double timePassed() {
  double time = zeitwert()-startTime;
  return time;
}


double cpuTimePassed() {
  double time = cpuTime()-startCPUTime;
  return time;
}


void randomPhi(int mode) {
  if (mode==0) {
    int I;
    for (I=0; I<4*L0*L1*L2*L3; I++) {
      phiField[I] = 2*(zufall()-0.5);
    }
  }
}


void randomPhi(double m, double s, double z) {
  int z0,z1,z2,z3;
  int count = 0;
  double fac = 1;
  for (z0=0; z0<L0; z0++) {
    for (z1=0; z1<L1; z1++) {
      for (z2=0; z2<L2; z2++) {
        for (z3=0; z3<L3; z3++) {
	  fac = 1;
	  if (((z0+z1+z2+z3)%2)==1) fac = -1;
 
          phiField[count+0] = m + s*fac + 2*z*(zufall()-0.5);
          phiField[count+1] = 0;
          phiField[count+2] = 0;
          phiField[count+3] = 0;
          count+=4;
	}
      }
    }
  }
}


void randomInput() {
  for (int I=0; I<fOps->getVectorLengthXtrSize(); I++) {
    input[I].x = NaN;
    input[I].y = NaN;    
  }
  for (int I=0; I<fOps->getVectorLength(); I++) {
    input[I] = Complex(2*(zufall()-0.5),2*(zufall()-0.5));
  }
  fOps->transformToXtraSizeArray(input, input);
}


void randomOutput() {
  for (int I=0; I<fOps->getVectorLengthXtrSize(); I++) {
    output[I].x = NaN;
    output[I].y = NaN;    
  }
  for (int I=0; I<fOps->getVectorLength(); I++) {
    output[I] = Complex(2*(zufall()-0.5),2*(zufall()-0.5));
  }
  fOps->transformToXtraSizeArray(output, output);
}


void randomOutput2() {
  for (int I=0; I<fOps->getVectorLengthXtrSize(); I++) {
    output2[I].x = NaN;
    output2[I].y = NaN;    
  }
  for (int I=0; I<fOps->getVectorLength(); I++) {
    output2[I] = Complex(2*(zufall()-0.5),2*(zufall()-0.5));
  }
  fOps->transformToXtraSizeArray(output2, output2);
}


void randomInterim() {
  for (int I=0; I<fOps->getVectorLengthXtrSize(); I++) {
    interim[I].x = NaN;
    interim[I].y = NaN;    
  }
  for (int I=0; I<fOps->getVectorLength(); I++) {
    interim[I] = Complex(2*(zufall()-0.5),2*(zufall()-0.5));
  }
  fOps->transformToXtraSizeArray(interim, interim);
}


void randomInterim2() {
  for (int I=0; I<fOps->getVectorLengthXtrSize(); I++) {
    interim2[I].x = NaN;
    interim2[I].y = NaN;    
  }
  for (int I=0; I<fOps->getVectorLength(); I++) {
    interim2[I] = Complex(2*(zufall()-0.5),2*(zufall()-0.5));
  }
  fOps->transformToXtraSizeArray(interim2, interim2);
}


void testGaussGenerator() {
  printf("\nPerforming Gauss test...\n");
  fOps->setSize(20,20,20,20);
  Complex* output = fOps->createFermionVector();
  fOps->fillGaussRandomVector(output,-1);

  int res[1001];
  int I;
  for (I=0; I<1001; I++) res[I] = 0;
  for (I=0; I<fOps->getVectorLength(); I++) {
    int p = (int)(500+(100*output[I].x));
    if (p<0) p=0;
    if (p>1001) p = 1001;
    res[p]++;
    p = (int)(500+(100*output[I].y));
    if (p<0) p=0;
    if (p>1001) p = 1001;
    res[p]++;    
  }
  
  FILE* file;
  file = fopen("data/GaussTest.dat","w");
  for (I=0; I<1001; I++) {
    fprintf(file,"%f %d\n",(I-500.0)/100.0, res[I]);
  }
  fclose(file);
  fOps->destroyFermionVector(output);
  fOps->setSize(L0,L1,L2,L3);
  printf("ready.\n");
}


void testSSECopy() {
  printf("\nTesting SSE-Copying...\n");
  randomInput();

  int VL = fOps->getVectorLength();
  int count = 0;
  
  fOps->transformFromXtraSizeArray(input, input);

  startTimer();
  count = 0;
  while (timePassed()<10) {
    SSE_ZCopy(VL, input, 1, output2, 1);
    count++;
  }
  printf("C++: %d iterations within %1.2f seconds.\n",count,timePassed());
  
  
#ifdef useBLAS
  startTimer();
  count = 0;
  while (timePassed()<10) {
    cblas_zcopy(VL, input, 1, output, 1);
    count++;
  }
  printf("BLAS: %d iterations within %1.2f seconds.\n",count,timePassed());

  int I;
  double diff = 0;
  for (I=0; I<VL; I++) diff += sqr(output[I].x-output2[I].x) + sqr(output[I].y-output2[I].y);
  diff = sqrt(diff);
  printf("Difference between results: %1.15f\n",diff);
#endif
  
  fOps->transformToXtraSizeArray(output, output);
  fOps->transformToXtraSizeArray(output2, output2);  
}



void testSSEScalarProduct() {
  printf("\nTesting SSE-Scalar-Products...\n");
  randomInput();
  randomOutput();

  Complex* res = new Complex[2];

  int count = 0;
  
  startTimer();
  count = 0;
  while (timePassed()<10) {
    SSE_ComplexScalarProduct(L0,L1,L2,L3,xtraSize1,xtraSize2,xtraSize3, input, output, res[0]);
    count++;
  }
  printf("SSE-Improved: %d iterations (scalar product) within %1.2f seconds.\n",count,timePassed());

#ifdef useBLAS
  int VL = fOps->getVectorLength();

  startTimer();
  count = 0;
  fOps->transformFromXtraSizeArray(input, input);
  fOps->transformFromXtraSizeArray(output, output);
  while (timePassed()<10) {
    cblas_zdotc_sub(VL, input, 1, output, 1, &(res[1]));
    count++;
  }
  fOps->transformToXtraSizeArray(input, input);
  fOps->transformToXtraSizeArray(output, output);
  printf("BLAS: %d iterations (scalar product) within %1.2f seconds.\n",count,timePassed());
  printf("Difference between results = %1.15f \n",sqrt(sqr(res[0].x-res[1].x)+sqr(res[0].y-res[1].y)));
#endif
 
  startTimer();
  count = 0;
  while (timePassed()<10) {
    SSE_ComplexSquareNorm(L0,L1,L2,L3,xtraSize1,xtraSize2,xtraSize3, input, res[0].x);
    count++;
  }
  printf("SSE-Improved: %d iterations (norm) within %1.2f seconds.\n",count,timePassed());

#ifdef useBLAS
  startTimer();
  count = 0;
  fOps->transformFromXtraSizeArray(input, input);
  fOps->transformFromXtraSizeArray(output, output);
  while (timePassed()<10) {
    cblas_zdotc_sub(VL, input, 1, input, 1, &(res[1]));
    count++;
  }
  fOps->transformToXtraSizeArray(input, input);
  fOps->transformToXtraSizeArray(output, output);
  printf("BLAS: %d iterations (norm) within %1.2f seconds.\n",count,timePassed());
  printf("Difference between results = %1.15f \n",sqrt(sqr(res[0].x-res[1].x)));
#endif

  delete[] res;
}


void testDistributedVectorCopy(int opMode, bool enforceCoreTie) {
  printf("\nTesting Distributed Vector-Copy...\n");
  randomInput();
  randomOutput();
  int VL = fOps->getVectorLength();

  fOps->activateMultiThreadedOps(opMode, enforceCoreTie);
  MultiThreadedOperations* threadedOps = fOps->getMultiThreadedOps();
  DistributedMemoryObject* memObj1 = threadedOps->allocateDistibutedFermionVector(L0, L1, L2, L3);
  DistributedMemoryObject* memObj2 = threadedOps->allocateDistibutedFermionVector(L0, L1, L2, L3);
  threadedOps->copyFermionVectorToDistributedFermionVector(input, memObj1, L0, L1, L2, L3);
  threadedOps->copyFermionVectorToDistributedFermionVector(output, memObj2, L0, L1, L2, L3);

  int count = 0;
  count = 0;
  threadedOps->copyFermionVector(memObj1, memObj2, L0, L1, L2, L3);
  startTimer();
  while (timePassed()<10) {
    for (int I=0; I<10; I++) {
      threadedOps->copyFermionVector(memObj1, memObj2, L0, L1, L2, L3);
    }
    count+=10;
  }
  printf("Distributed Vector-Copy: %d iterations within %1.2f seconds ==> Bandwidth: %1.2f GB/sec.\n",count,timePassed(),2*16*VL*(count/timePassed())/(1024*1024*1024));
  fOps->copyDistributedFermionVectorToFermionVector(memObj2, output);
  
  fOps->transformFromXtraSizeArray(input, input);  
  fOps->transformFromXtraSizeArray(output, output);
  memcpy(output2, input, VL*16);  
  startTimer();
  count = 0;
  while (timePassed()<10) {
    for (int I=0; I<10; I++) {
      memcpy(output2, input, VL*16);  
    }
    count+=10;
  }
  printf("MEMCPY: %d iterations within %1.2f seconds. ==> Bandwidth: %1.2f GB/sec.\n",count,timePassed(),2*16*VL*(count/timePassed())/(1024*1024*1024));
  double dif = 0;  
  double norm = 0;
  for (int I=0; I<VL; I++) {
    dif += sqr(output[I].x-output2[I].x) + sqr(output[I].y-output2[I].y);
    norm += sqr(output[I].x) + sqr(output[I].y);
  }
  dif = sqrt(dif);
  norm = sqrt(norm);
  printf("Difference between results = %1.15f (relative: %1.3e) \n",dif, dif/norm);
  
#ifdef useBLAS
  cblas_zcopy(VL, input, 1, output2, 1);  
  startTimer();
  count = 0;
  while (timePassed()<10) {
    for (int I=0; I<10; I++) {
      cblas_zcopy(VL, input, 1, output2, 1);
    }
    count+=10;
  }
  printf("BLAS: %d iterations within %1.2f seconds. ==> Bandwidth: %1.2f GB/sec.\n",count,timePassed(),2*16*VL*(count/timePassed())/(1024*1024*1024));
  dif = 0;  
  norm = 0;
  for (int I=0; I<VL; I++) {
    dif += sqr(output[I].x-output2[I].x) + sqr(output[I].y-output2[I].y);
    norm += sqr(output[I].x) + sqr(output[I].y);
  }
  dif = sqrt(dif);
  norm = sqrt(norm);
  printf("Difference between results = %1.15f (relative: %1.3e) \n",dif, dif/norm);
#endif
 
  fOps->deactivateMultiThreadedOps();
  delete memObj1;
  delete memObj2;
}

void testDistributedVectorNorm(int opMode, bool enforceCoreTie) {
  printf("\nTesting Distributed Vector-Norm...\n");
  randomInput();
  randomOutput();

  Complex* res = new Complex[2];

  MultiThreadedOperations* threadedOps = new MultiThreadedOperations(opMode, enforceCoreTie);
system("numastat"); 
  DistributedMemoryObject* memObj1 = threadedOps->allocateDistibutedFermionVector(L0, L1, L2, L3);
system("numastat");  
  DistributedMemoryObject* memObj2 = threadedOps->allocateDistibutedFermionVector(L0, L1, L2, L3);
  threadedOps->copyFermionVectorToDistributedFermionVector(input, memObj1, L0, L1, L2, L3);
  threadedOps->copyFermionVectorToDistributedFermionVector(output, memObj2, L0, L1, L2, L3);

  int count = 0;
  count = 0;
  res[0].y=0;
  startTimer();
  while (timePassed()<10) {
    threadedOps->vectorNormOfFermionVector(memObj1, res[0].x, L0, L1, L2, L3);
    count++;
  }
  printf("Distributed Vector-Norm: %d iterations within %1.2f seconds.\n",count,timePassed());

#ifdef useBLAS
  int VL = fOps->getVectorLength();

  fOps->transformFromXtraSizeArray(input, input);
  fOps->transformFromXtraSizeArray(output, output);
  startTimer();
  count = 0;
  while (timePassed()<10) {
    cblas_zdotc_sub(VL, input, 1, input, 1, &(res[1]));
    count++;
  }
  printf("BLAS: %d iterations (scalar product) within %1.2f seconds.\n",count,timePassed());
  double dif = sqrt(sqr(res[0].x-res[1].x)+sqr(res[0].y-res[1].y));
  double norm = sqrt(sqr(res[0].x)+sqr(res[0].y));
  printf("Difference between results = %1.15f (relative: %1.3e) \n",dif, dif/norm);
#endif
 
  delete[] res;
  delete threadedOps;
  delete memObj1;
  delete memObj2;
}


void testDistributedScalarProduct(int opMode, bool enforceCoreTie) {
  printf("\nTesting Distributed Scalar-Product...\n");
  randomInput();
  randomOutput();

  Complex* res = new Complex[2];

  MultiThreadedOperations* threadedOps = new MultiThreadedOperations(opMode, enforceCoreTie);
  DistributedMemoryObject* memObj1 = threadedOps->allocateDistibutedFermionVector(L0, L1, L2, L3);
  DistributedMemoryObject* memObj2 = threadedOps->allocateDistibutedFermionVector(L0, L1, L2, L3);
  threadedOps->copyFermionVectorToDistributedFermionVector(input, memObj1, L0, L1, L2, L3);
  threadedOps->copyFermionVectorToDistributedFermionVector(output, memObj2, L0, L1, L2, L3);

  int count = 0;
  startTimer();
  count = 0;
  res[0].y=0;
  while (timePassed()<10) {
    threadedOps->scalarProductOfFermionVectors(memObj1, memObj2, res[0], L0, L1, L2, L3);
    count++;
  }
  printf("Distributed Scalar-Product: %d iterations within %1.2f seconds.\n",count,timePassed());

#ifdef useBLAS
  int VL = fOps->getVectorLength();

  fOps->transformFromXtraSizeArray(input, input);
  fOps->transformFromXtraSizeArray(output, output);
  startTimer();
  count = 0;
  while (timePassed()<10) {
    cblas_zdotc_sub(VL, input, 1, output, 1, &(res[1]));
    count++;
  }
  printf("BLAS: %d iterations (scalar product) within %1.2f seconds.\n",count,timePassed());
  double dif = sqrt(sqr(res[0].x-res[1].x)+sqr(res[0].y-res[1].y));
  double norm = sqrt(sqr(res[0].x)+sqr(res[0].y));
  printf("Difference between results = %1.15f (relative: %1.3e) \n",dif, dif/norm);
#endif
 
  delete[] res;
  delete threadedOps;
  delete memObj1;
  delete memObj2;
}	  


void testDistributedVectorAddition(int opMode, bool enforceCoreTie) {
  printf("\nTesting Distributed Vector-Addition...\n");
  randomInput();
  randomOutput();
  int VL = fOps->getVectorLength();
  int VLXtr = fOps->getVectorLengthXtrSize();
  for (int I=0; I<VLXtr; I++) output2[I] = output[I];

  Complex alpha(1.7,-2.3);

  MultiThreadedOperations* threadedOps = new MultiThreadedOperations(opMode, enforceCoreTie);
  DistributedMemoryObject* memObj1 = threadedOps->allocateDistibutedFermionVector(L0, L1, L2, L3);
  DistributedMemoryObject* memObj2 = threadedOps->allocateDistibutedFermionVector(L0, L1, L2, L3);
  threadedOps->copyFermionVectorToDistributedFermionVector(input, memObj1, L0, L1, L2, L3);
  threadedOps->copyFermionVectorToDistributedFermionVector(output, memObj2, L0, L1, L2, L3);

  threadedOps->vectorAdditionOfFermionVectors(memObj1, memObj2, alpha, L0, L1, L2, L3);
  threadedOps->copyDistributedFermionVectorToFermionVector(memObj2, output2, L0, L1, L2, L3);
  
#ifdef useBLAS
  fOps->transformFromXtraSizeArray(input, input);
  fOps->transformFromXtraSizeArray(output, output);
  fOps->transformFromXtraSizeArray(output2, output2);  

  cblas_zaxpy(VL, &alpha, input, 1, output, 1);
  double diff = 0;
  double norm = 0;
  for (int I=0; I<VL; I++) {
    diff += sqr(output[I].x - output2[I].x) + sqr(output[I].y - output2[I].y);
    norm += sqr(output[I].x) + sqr(output[I].y);
  }
  diff = sqrt(diff);
  norm = sqrt(norm);
  fOps->transformToXtraSizeArray(input, input);
  fOps->transformToXtraSizeArray(output, output);
  fOps->transformToXtraSizeArray(output2, output2);  
#endif


  int count = 0;
  startTimer();
  count = 0;
  while (timePassed()<10) {
    threadedOps->vectorAdditionOfFermionVectors(memObj1, memObj2, alpha, L0, L1, L2, L3);
    count++;
  }
  printf("Distributed Vector-Addition: %d iterations within %1.2f seconds.\n",count,timePassed());
  
  
#ifdef useBLAS
  fOps->transformFromXtraSizeArray(input, input);
  fOps->transformFromXtraSizeArray(output, output);
  startTimer();
  count = 0;
  while (timePassed()<10) {
    cblas_zaxpy(VL, &alpha, input, 1, output, 1);
    count++;
  }
  printf("BLAS: %d iterations (vector addition) within %1.2f seconds.\n",count,timePassed());
  printf("Difference between results = %1.15f (relative: %1.3e) \n",diff, diff/norm);
#endif
 
  delete threadedOps;
  delete memObj1;
  delete memObj2;
}	  


void testDistributedFourierTransformation(int opMode, bool enforceCoreTie, int threadCount) {
  printf("\nTesting Distributed Fourier-Transformation...\n");
  randomInput();
  randomOutput();
  int VL = fOps->getVectorLength();

  fOps->setxFFTusage(false);
  fOps->performFFT(input, output, ExtremeFFT4D_Forward);
  startTimer();
  int count = 0;
  while (timePassed()<10) {
    int I;
    for (I=0; I<10; I++) {
      fOps->performFFT(input, output, ExtremeFFT4D_Forward);
    }
    count+=10;
  }
  printf("%d iterations within %1.2f (CPU: %1.2f) seconds (FFTW) --> %d cycles.\n",count,timePassed(), cpuTimePassed(), (int)(timePassed()*2.4E9/count));


  MultiThreadedOperations* threadedOps = new MultiThreadedOperations(opMode, enforceCoreTie);
  DistributedMemoryObject* memObj1 = threadedOps->allocateDistibutedFermionVector(L0, L1, L2, L3);
  DistributedMemoryObject* memObj2 = threadedOps->allocateDistibutedFermionVector(L0, L1, L2, L3);
  threadedOps->copyFermionVectorToDistributedFermionVector(input, memObj1, L0, L1, L2, L3);
  threadedOps->copyFermionVectorToDistributedFermionVector(output, memObj2, L0, L1, L2, L3);

  threadedOps->perform_FFTWFourierTransformationOfFermionVector(memObj1, memObj2, ExtremeFFT4D_Forward, L0, L1, L2, L3);
  threadedOps->copyDistributedFermionVectorToFermionVector(memObj2, output2, L0, L1, L2, L3);  
  startTimer();
  count = 0;
  while (timePassed()<10) {
    int I;
    for (I=0; I<10; I++) {
      threadedOps->perform_FFTWFourierTransformationOfFermionVector(memObj1, memObj2, ExtremeFFT4D_Forward, L0, L1, L2, L3);
    }
    count+=10;
  }
  printf("%d iterations within %1.2f seconds (Distributed FFTW) --> %d cycles.\n",count,timePassed(), (int)(timePassed()*2.4E9/count));

  
  threadedOps->perform_xFFTFourierTransformationOfFermionVector(memObj1, memObj2, ExtremeFFT4D_Forward, threadCount, L0, L1, L2, L3);
  threadedOps->copyDistributedFermionVectorToFermionVector(memObj2, interim, L0, L1, L2, L3);  
  startTimer();
  count = 0;
  while (timePassed()<10) {
    int I;
    for (I=0; I<10; I++) {
      threadedOps->perform_xFFTFourierTransformationOfFermionVector(memObj1, memObj2, ExtremeFFT4D_Forward, threadCount, L0, L1, L2, L3);
    }
    count+=10;
  }
  printf("%d iterations within %1.2f seconds (Distributed xFFT) --> %d cycles.\n",count,timePassed(), (int)(timePassed()*2.4E9/count));
  
  fOps->transformFromXtraSizeArray(output,output);
  fOps->transformFromXtraSizeArray(output2,output2);
  fOps->transformFromXtraSizeArray(interim,interim);
  double diffFFTW = 0;
  double diffxFFT = 0;  
  for (int I=0; I<VL; I++) {
    diffFFTW += sqr(output[I].x-output2[I].x) + sqr(output[I].y-output2[I].y);
    diffxFFT += sqr(output[I].x-interim[I].x) + sqr(output[I].y-interim[I].y);    
  }
  printf("Difference FFTW: %1.15f, xFFT: %1.15f\n",diffFFTW, diffxFFT);
 
  delete threadedOps;
  delete memObj1;
  delete memObj2;
}	  


void testDistributedMulWithDerivativesOfMatB(int opMode, bool enforceCoreTie, double split) {
  printf("\nTesting Distributed Multiplication with Derivatives of Mat B with split=%f...\n", split);
  randomInput();
  randomOutput();
  int VL = fOps->getVectorLength();

  MultiThreadedOperations* threadedOps = new MultiThreadedOperations(opMode, enforceCoreTie);
  DistributedMemoryObject* memObj1 = threadedOps->allocateDistibutedFermionVector(L0, L1, L2, L3);
  DistributedMemoryObject* memObj2 = threadedOps->allocateDistibutedFermionVector(L0, L1, L2, L3);
  DistributedMemoryObject* memObj3 = threadedOps->allocateDistibutedMemory(4*L0*L1*L2*L3);
  threadedOps->copyFermionVectorToDistributedFermionVector(input, memObj1, L0, L1, L2, L3);
  threadedOps->copyFermionVectorToDistributedFermionVector(output, memObj2, L0, L1, L2, L3);

  startTimer();
  int count = 0;
  while (timePassed()<10) {
    int I;
    for (I=0; I<10; I++) {
      fOps->executeMultiplicationVectorWithDerivativesOfB(input, output, output2);
    }
    count+=10;
  }
  printf("%d iterations within %1.2f (CPU: %1.2f) seconds --> %d cycles.\n",count,timePassed(), cpuTimePassed(), (int)(timePassed()*2.4E9/count));

  startTimer();
  count = 0;
  while (timePassed()<10) {
    int I;
    for (I=0; I<10; I++) {
      threadedOps->perform_MultiplicationVectorWithDerivativesOfB(memObj1, memObj2, L0, L1, L2, L3, split, memObj3);
    }    
    count+=10;
  }
  printf("%d iterations within %1.2f seconds (Distributed) --> %d cycles.\n",count,timePassed(), (int)(timePassed()*2.4E9/count));

  threadedOps->calcSumOfUniformlyDistributedComplexVectors(memObj3, output, 4*L0*L1*L2*L3);
  
  double diff = 0;
  for (int I=0; I<VL/2; I++) {
    diff += sqr(output[I].x-output2[I].x) + sqr(output[I].y-output2[I].y);
  }
  printf("Difference %1.15f:\n",diff);
 
  delete threadedOps;
  delete memObj1;
  delete memObj2;
  delete memObj3;
}


void testDistributedMMDaggerxQuasiHermiteanNeubergerWithChi(int OpMode, bool enforceCoreTie, int threadCount, double yN, bool useR, bool inFourierSpace) {
  printf("\nTesting Distributed MMdag in quasi-hermitean mode with yN = %1.2f, useR = %d, and inFourierSpace=%d...\n",yN, useR,inFourierSpace);  
  randomPhi(0);
  randomInput();
  randomOutput();
  int VL = fOps->getVectorLength();
  
  fOps->setYukawaCoupling(yN);
  fOps->setxFFT_DistributedFFT_ThreadCount(threadCount);
  fOps->setRPreconditioner(true, 1.0, 1.0);
  double rho,r;
  fOps->getDiracParameters(rho, r);
  
  fOps->setxFFTusage(false);
  fOps->executeFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(input, interim, phiField, useR, inFourierSpace);           
  startTimer();
  int count = 0;
  while (timePassed()<10) {
    for (int I=0; I<10; I++) {
      fOps->executeFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(input, output, phiField, useR, inFourierSpace);   
    }
    count+=10;
  }
  printf("%d iterations (FFTW, serial) within %1.2f seconds.\n",count,timePassed());


  fOps->activateMultiThreadedOps(OpMode, enforceCoreTie);  
  DistributedMemoryObject* memObj1 = fOps->createDistributedFermionVector();
  DistributedMemoryObject* memObj2 = fOps->createDistributedFermionVector();
  DistributedMemoryObject* memObj3 = fOps->createDistributedUniformLatticeComplexVector(2); 
  fOps->copyFermionVectorToDistributedFermionVector(input, memObj1);
  fOps->copyFermionVectorToDistributedFermionVector(output, memObj2);
  fOps->copyPhiFieldUniformlyToDistributedPhiField(phiField, memObj3);
  
  fOps->executeDistributedFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(memObj1, memObj2, memObj3, useR, inFourierSpace);
  startTimer();
  count = 0;
  while (timePassed()<10) {
    for (int I=0; I<10; I++) {
      fOps->executeDistributedFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(memObj1, memObj2, memObj3, useR, inFourierSpace);
    }
    count+=10;
  }
  printf("%d iterations (FFTW, Distributed) within %1.2f seconds.\n",count,timePassed());
  fOps->copyDistributedFermionVectorToFermionVector(memObj2, output); 
  fOps->transformFromXtraSizeArray(output, interim2);
  
  fOps->setxFFTusage(true);
  fOps->executeFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(input, output2, phiField, useR, inFourierSpace);       
  fOps->executeDistributedFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(memObj1, memObj2, memObj3, useR, inFourierSpace);
  startTimer();
  count = 0;
  while (timePassed()<10) {
    for (int I=0; I<10; I++) {
      fOps->executeDistributedFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(memObj1, memObj2, memObj3, useR, inFourierSpace);
    }
    count+=10;
  }
  printf("%d iterations (xFFT, Distributed) within %1.2f seconds.\n",count,timePassed());
  fOps->copyDistributedFermionVectorToFermionVector(memObj2, output); 
  fOps->transformFromXtraSizeArray(output, output);
  fOps->transformFromXtraSizeArray(output2, output2); 
  fOps->transformFromXtraSizeArray(interim, interim); 
  
  double diffFFTWDis = 0;
  double diffxFFT = 0;
  double diffxFFTDis = 0;
  for (int I=0; I<VL; I++) {
    diffFFTWDis += sqr(interim2[I].x-interim[I].x) + sqr(interim2[I].y-interim[I].y);
    diffxFFT += sqr(output2[I].x-interim[I].x) + sqr(output2[I].y-interim[I].y);
    diffxFFTDis += sqr(output[I].x-interim[I].x) + sqr(output[I].y-interim[I].y);
  }
  printf("Differences relatice to FFTW-serial: FFTW-distributed: %1.15f, xFFT-serial: %1.15f, xFFT-distributed: %1.15f\n",diffFFTWDis,diffxFFT,diffxFFTDis);
 
  delete memObj1;
  delete memObj2;
  delete memObj3;
  fOps->deactivateMultiThreadedOps();
}


void testDistributedPiModeRemoverOperatorApplication(int opMode, bool enforceCoreTie) {
  printf("\nTesting Distributed Application of Pi-Mode-Remover Operator...\n");
  randomInput();
  randomOutput();
  int VL = fOps->getVectorLength();

  fOps->activateMultiThreadedOps(opMode, enforceCoreTie);  
  DistributedMemoryObject* memObj1 = fOps->createDistributedFermionVector();
  DistributedMemoryObject* memObj2 = fOps->createDistributedFermionVector();
  fOps->copyFermionVectorToDistributedFermionVector(input, memObj1);
  fOps->copyFermionVectorToDistributedFermionVector(output, memObj2);
  
  
  startTimer();
  int count = 0;
  while (timePassed()<10) {
    int I;
    for (I=0; I<10; I++) {
      fOps->executePiModeRemoverOperator(input, output2, true);
    }
    count+=10;
  }
  printf("%d iterations within %1.2f (CPU: %1.2f) seconds --> %d cycles.\n",count,timePassed(), cpuTimePassed(), (int)(timePassed()*2.4E9/count));

  startTimer();
  count = 0;
  while (timePassed()<10) {
    int I;
    for (I=0; I<10; I++) {
      fOps->executeDistributedPiModeRemoverOperator(memObj1, memObj2, true);
    }    
    count+=10;
  }
  printf("%d iterations within %1.2f seconds (Distributed) --> %d cycles.\n",count,timePassed(), (int)(timePassed()*2.4E9/count));

  fOps->copyDistributedFermionVectorToFermionVector(memObj2, output); 
  fOps->transformFromXtraSizeArray(output, output);
  fOps->transformFromXtraSizeArray(output2, output2); 
  
  double diff = 0;
  for (int I=0; I<VL; I++) {
    diff += sqr(output[I].x-output2[I].x) + sqr(output[I].y-output2[I].y);
  }
  printf("Difference %1.15f:\n",diff);
 
  delete memObj1;
  delete memObj2;
  fOps->deactivateMultiThreadedOps();
}



void testDistributedRPreconditionerOperatorApplication(int opMode, bool enforceCoreTie) {
  printf("\nTesting Distributed Application of R-Preconditioner Operator...\n");
  randomInput();
  randomOutput();
  int VL = fOps->getVectorLength();

  fOps->activateMultiThreadedOps(opMode, enforceCoreTie);  
  DistributedMemoryObject* memObj1 = fOps->createDistributedFermionVector();
  DistributedMemoryObject* memObj2 = fOps->createDistributedFermionVector();
  fOps->copyFermionVectorToDistributedFermionVector(input, memObj1);
  fOps->copyFermionVectorToDistributedFermionVector(output, memObj2);

  
  startTimer();
  int count = 0;
  while (timePassed()<10) {
    int I;
    for (I=0; I<10; I++) {
      fOps->executeRPreconditionerMatrixMultiplication(input, output2, false, false, true);
    }
    count+=10;
  }
  printf("%d iterations within %1.2f (CPU: %1.2f) seconds --> %d cycles.\n",count,timePassed(), cpuTimePassed(), (int)(timePassed()*2.4E9/count));

  startTimer();
  count = 0;
  while (timePassed()<10) {
    int I;
    for (I=0; I<10; I++) {
      fOps->executeDistributedRPreconditionerMatrixMultiplication(memObj1, memObj2, false, false, true);
    }    
    count+=10;
  }
  printf("%d iterations within %1.2f seconds (Distributed) --> %d cycles.\n",count,timePassed(), (int)(timePassed()*2.4E9/count));

  fOps->copyDistributedFermionVectorToFermionVector(memObj2, output); 
  fOps->transformFromXtraSizeArray(output, output);
  fOps->transformFromXtraSizeArray(output2, output2); 
  
  double diff = 0;
  for (int I=0; I<VL; I++) {
    diff += sqr(output[I].x-output2[I].x) + sqr(output[I].y-output2[I].y);
  }
  printf("Difference %1.15f:\n",diff);
 
  delete memObj1;
  delete memObj2;
  fOps->deactivateMultiThreadedOps();
}


void testDistributedDiracOperatorApplication(int opMode, bool enforceCoreTie) {
  printf("\nTesting Distributed Application of Neuberger-Dirac Operator...\n");
  randomInput();
  randomOutput();
  int VL = fOps->getVectorLength();

  MultiThreadedOperations* threadedOps = new MultiThreadedOperations(opMode, enforceCoreTie);
  DistributedMemoryObject* memObj1 = threadedOps->allocateDistibutedFermionVector(L0, L1, L2, L3);
  DistributedMemoryObject* memObj2 = threadedOps->allocateDistibutedFermionVector(L0, L1, L2, L3);
  DistributedMemoryObject* memObj3 = threadedOps->allocateDistibutedMemory(2*4*128);
  DistributedMemoryObject* memObj4 = threadedOps->allocateDistibutedMemory(L0*L1*L2*L3);
  threadedOps->copyFermionVectorToDistributedFermionVector(input, memObj1, L0, L1, L2, L3);
  threadedOps->copyFermionVectorToDistributedFermionVector(output, memObj2, L0, L1, L2, L3);
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(fOps->NeubergerDiracOpApplicationSinPStdData, memObj3, 2*4*128);
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(fOps->NeubergerDiracOpApplicationEWDataRescaled, memObj4, L0*L1*L2*L3);
  
  
  startTimer();
  int count = 0;
  while (timePassed()<10) {
    int I;
    for (I=0; I<10; I++) {
      SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace(L0, L1, L2, L3, input, output2, fOps->NeubergerDiracOpApplicationSinPStdData, fOps->NeubergerDiracOpApplicationEWDataRescaled);    
    }
    count+=10;
  }
  printf("%d iterations within %1.2f (CPU: %1.2f) seconds --> %d cycles.\n",count,timePassed(), cpuTimePassed(), (int)(timePassed()*2.4E9/count));

  startTimer();
  count = 0;
  while (timePassed()<10) {
    int I;
    for (I=0; I<10; I++) {
      threadedOps->perform_DiracTypeOperatorMultiplicationInFourierSpace(memObj1, memObj2, L0, L1, L2, L3, memObj3, memObj4);
    }    
    count+=10;
  }
  printf("%d iterations within %1.2f seconds (Distributed) --> %d cycles.\n",count,timePassed(), (int)(timePassed()*2.4E9/count));

  threadedOps->copyDistributedFermionVectorToFermionVector(memObj2, output, L0, L1, L2, L3);
  fOps->transformFromXtraSizeArray(output, output);
  fOps->transformFromXtraSizeArray(output2, output2); 
  
  double diff = 0;
  for (int I=0; I<VL; I++) {
    diff += sqr(output[I].x-output2[I].x) + sqr(output[I].y-output2[I].y);
  }
  printf("Difference %1.15f:\n",diff);
 
  delete threadedOps;
  delete memObj1;
  delete memObj2;
  delete memObj3;
  delete memObj4;
}


void testDistributedMultiCombinedOperator_VA1_RorCopy_VAc_RorCopy_D_MatrixMultiplication(int opMode, bool enforceCoreTie, int threadCount, bool Ddag, bool usexFFT, bool RprecUse) {
  printf("\nTesting Distributed Application of Combined Operator VA+R+VA+R+D...\n");
  randomInput();
  randomOutput();
  randomInterim();
  randomInterim2();
  randomOutput2();
  
  int VL = fOps->getVectorLength();
  int VLxtr = fOps->getVectorLengthXtrSize();



  fOps->setYukawaCoupling(1.7);
  fOps->setxFFT_DistributedFFT_ThreadCount(threadCount);
  fOps->setRPreconditioner(RprecUse, 1.34, 0.45);

  fOps->activateMultiThreadedOps(opMode, enforceCoreTie);  
  fOps->setxFFTusage(false);  
  MultiThreadedOperations* threadedOps = fOps->getMultiThreadedOps();
  DistributedMemoryObject* memObj1 = threadedOps->allocateDistibutedFermionVector(L0, L1, L2, L3);
  DistributedMemoryObject* memObj2 = threadedOps->allocateDistibutedFermionVector(L0, L1, L2, L3);
  DistributedMemoryObject* memObj3 = threadedOps->allocateDistibutedFermionVector(L0, L1, L2, L3);
  DistributedMemoryObject* memObj4 = threadedOps->allocateDistibutedFermionVector(L0, L1, L2, L3);
  threadedOps->copyFermionVectorToDistributedFermionVector(input, memObj1, L0, L1, L2, L3);
  threadedOps->copyFermionVectorToDistributedFermionVector(output, memObj2, L0, L1, L2, L3);
  threadedOps->copyFermionVectorToDistributedFermionVector(interim, memObj3, L0, L1, L2, L3);
  threadedOps->copyFermionVectorToDistributedFermionVector(interim2, memObj4, L0, L1, L2, L3);
  Complex alpha1(1.0,0.0);  
  Complex alpha2(1.5,3.7);
  


#ifdef useBLAS
  fOps->transformFromXtraSizeArray(interim, interim);
  fOps->transformFromXtraSizeArray(input, input);
  cblas_zaxpy(VL, &alpha1, input, 1, interim, 1);
  fOps->transformToXtraSizeArray(interim, interim);  
  fOps->transformToXtraSizeArray(input, input);  
#endif

  fOps->executeRPreconditionerMatrixMultiplication(interim, interim, false, false, true);

#ifdef useBLAS
  fOps->transformFromXtraSizeArray(interim, interim);
  fOps->transformFromXtraSizeArray(interim2, interim2);
  cblas_zaxpy(VL, &alpha2, interim2, 1, interim, 1);
  fOps->transformToXtraSizeArray(interim, interim);  
  fOps->transformToXtraSizeArray(interim2, interim2);  
#endif

  if (RprecUse) {
    fOps->executeRPreconditionerMatrixMultiplication(interim, input, false, false, true);
  } else {
    for (int I=0; I<VLxtr; I++) {
      input[I].x = interim[I].x;
      input[I].y = interim[I].y;      
    }
  }
  
  if (Ddag) {
    fOps->executeDiracDaggerMatrixMultiplication(input, output, true);
  } else {
    fOps->executeDiracMatrixMultiplication(input, output, true);
  }
  
  
  fOps->setxFFTusage(usexFFT);
  fOps->executeDistributedMultiCombinedOperator_VA1_RorCopy_VAc_RorCopy_D_MatrixMultiplication(memObj1, memObj2, memObj3, memObj4, alpha2, Ddag, true);

  threadedOps->copyDistributedFermionVectorToFermionVector(memObj1, output2, L0, L1, L2, L3);
  fOps->transformFromXtraSizeArray(output2, output2);
  fOps->transformFromXtraSizeArray(input, input);
  double diff1 = 0;
  for (int I=0; I<VL; I++) {
    diff1 += sqr(input[I].x-output2[I].x) + sqr(input[I].y-output2[I].y);
  }
  threadedOps->copyDistributedFermionVectorToFermionVector(memObj2, output2, L0, L1, L2, L3);
  fOps->transformFromXtraSizeArray(output2, output2);
  fOps->transformFromXtraSizeArray(output, output);
  double diff2 = 0;
  for (int I=0; I<VL; I++) {
    diff2 += sqr(output[I].x-output2[I].x) + sqr(output[I].y-output2[I].y);
  }
  threadedOps->copyDistributedFermionVectorToFermionVector(memObj3, output2, L0, L1, L2, L3);
  fOps->transformFromXtraSizeArray(output2, output2);
  fOps->transformFromXtraSizeArray(interim, interim);
  double diff3 = 0;
  for (int I=0; I<VL; I++) {
    diff3 += sqr(interim[I].x-output2[I].x) + sqr(interim[I].y-output2[I].y);
  }
  threadedOps->copyDistributedFermionVectorToFermionVector(memObj4, output2, L0, L1, L2, L3);
  fOps->transformFromXtraSizeArray(output2, output2);
  fOps->transformFromXtraSizeArray(interim2, interim2);
  double diff4 = 0;
  for (int I=0; I<VL; I++) {
    diff4 += sqr(interim2[I].x-output2[I].x) + sqr(interim2[I].y-output2[I].y);
  }

  startTimer();
  int count = 0;
  while (timePassed()<10) {
    int I;
    for (I=0; I<10; I++) {
      fOps->executeDistributedMultiCombinedOperator_VA1_RorCopy_VAc_RorCopy_D_MatrixMultiplication(memObj1, memObj2, memObj3, memObj4, alpha2, Ddag, true);
    }    
    count+=10;
  }
  printf("%d iterations within %1.2f seconds (Distributed) --> %d cycles.\n",count,timePassed(), (int)(timePassed()*2.4E9/count));
  printf("Differences v1: %1.15f, v2: %1.15f, v3: %1.15f, v4: %1.15f\n", diff1, diff2, diff3, diff4);

  delete memObj1;
  delete memObj2;
  delete memObj3;
  delete memObj4;
}


void testDistributedYukawaCouplingMatB(int opMode, bool enforceCoreTie, bool dag, double y, double split) {
  printf("\nTesting Distributed YukawaCouplingMatB with dag = %d, y=%f, split=%f...\n", dag, y, split);
  randomInput();
  randomOutput();
  int VL = fOps->getVectorLength();

  MultiThreadedOperations* threadedOps = new MultiThreadedOperations(opMode, enforceCoreTie);
  DistributedMemoryObject* memObj1 = threadedOps->allocateDistibutedFermionVector(L0, L1, L2, L3);
  DistributedMemoryObject* memObj2 = threadedOps->allocateDistibutedFermionVector(L0, L1, L2, L3);
  DistributedMemoryObject* memObj3 = threadedOps->allocateDistibutedMemory(2*L0*L1*L2*L3);
  threadedOps->copyFermionVectorToDistributedFermionVector(input, memObj1, L0, L1, L2, L3);
  threadedOps->copyFermionVectorToDistributedFermionVector(output, memObj2, L0, L1, L2, L3);
  threadedOps->copyComplexVectorToUniformlyDistributedMemory((Complex*) phiField, memObj3, 2*L0*L1*L2*L3);

  startTimer();
  int count = 0;
  while (timePassed()<10) {
    int I;
    if (dag) {
      for (I=0; I<10; I++) {
        perform_yBsplitDagger(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, y, split, phiField, input, output2);
      }    
    } else {
      for (I=0; I<10; I++) {
        perform_yBsplit(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, y, split, phiField, input, output2);
      }
    }
    count+=10;
  }
  printf("%d iterations within %1.2f (CPU: %1.2f) seconds --> %d cycles.\n",count,timePassed(), cpuTimePassed(), (int)(timePassed()*2.4E9/count));


  startTimer();
  count = 0;
  while (timePassed()<10) {
    int I;
    if (dag) {
      for (I=0; I<10; I++) {
        threadedOps->perform_yBDaggerOperatorMultiplication(memObj1, memObj2, L0, L1, L2, L3, y, split, memObj3);
      }
    } else {
      for (I=0; I<10; I++) {
        threadedOps->perform_yBOperatorMultiplication(memObj1, memObj2, L0, L1, L2, L3, y, split, memObj3);
      }    
    }
    count+=10;
  }
  printf("%d iterations within %1.2f seconds (Distributed) --> %d cycles.\n",count,timePassed(), (int)(timePassed()*2.4E9/count));

  threadedOps->copyDistributedFermionVectorToFermionVector(memObj2, output, L0, L1, L2, L3);
  
  fOps->transformFromXtraSizeArray(output,output);
  fOps->transformFromXtraSizeArray(output2,output2);
  double diff = 0;
  for (int I=0; I<VL; I++) {
    diff += sqr(output[I].x-output2[I].x) + sqr(output[I].y-output2[I].y);
  }
  printf("Difference %1.15f:\n",diff);
 
  delete threadedOps;
  delete memObj1;
  delete memObj2;
  delete memObj3;
}


void testSSEVectorAddition() {
  printf("\nTesting SSE-Vector-Additions...\n");
  int VL = fOps->getVectorLength();
  int VLXtr = fOps->getVectorLengthXtrSize();
  int count = 0;
  randomInput();
  randomOutput();
  Complex alpha(input[0].x, input[0].y);
  int I;
  xSSE* xSSEObj = new xSSE();
  
  for (I=0; I<VLXtr; I++) output2[I] = output[I];
  fOps->transformFromXtraSizeArray(input, input);
  fOps->transformFromXtraSizeArray(output, output);
#ifdef useBLAS
  cblas_zaxpy(VL, &alpha, input, 1, output, 1);
#endif
  fOps->transformToXtraSizeArray(input, input);
  xSSEObj->xSSE_ComplexVectorAddition_Wrapper(input, output2, alpha, 8, L0,L1,L2,L3,xtraSize1,xtraSize2,xtraSize3);  
//  SSE_ComplexVectorAddition(L0,L1,L2,L3,xtraSize1,xtraSize2,xtraSize3, alpha, input, output2);
  fOps->transformFromXtraSizeArray(output2, output2);
  double diff = 0;
#ifdef useBLAS
  for (I=0; I<VL; I++) {
    diff += sqr(output[I].x - output2[I].x) + sqr(output[I].y - output2[I].y);
  }
  diff = sqrt(diff);
#endif


  randomInput();
  randomOutput();
  for (I=0; I<VLXtr; I++) output2[I] = output[I];
  double n1 = 0;
  fOps->transformFromXtraSizeArray(input, input);
  fOps->transformFromXtraSizeArray(output, output);
#ifdef useBLAS
  cblas_zaxpy(VL, &alpha, input, 1, output, 1);
#endif
  fOps->transformToXtraSizeArray(input, input);
  SSE_ComplexVectorAdditionWithSquaredNorm(L0,L1,L2,L3,xtraSize1,xtraSize2,xtraSize3, alpha, input, output2, n1);
  fOps->transformFromXtraSizeArray(output2, output2);
  double diff2 = 0;
#ifdef useBLAS
  double n2 = 0;
  for (I=0; I<VL; I++) {
    diff2 += sqr(output[I].x - output2[I].x) + sqr(output[I].y - output2[I].y);
    n2 += sqr(output[I].x) + sqr(output[I].y);
  }
  diff2 = sqrt(diff2);
#endif


  randomInput();
  randomOutput();
  for (I=0; I<VLXtr; I++) output2[I] = output[I];
  fOps->transformFromXtraSizeArray(input, input);
  fOps->transformFromXtraSizeArray(output, output);
  for (I=0; I<VL; I++) {
    output[I].x = alpha.x*output[I].x + input[I].x;
    output[I].y = alpha.x*output[I].y + input[I].y;
  }
  fOps->transformToXtraSizeArray(input, input);
  SSE_ComplexVectorAdditionSPECIAL1(L0,L1,L2,L3,xtraSize1,xtraSize2,xtraSize3, alpha.x, input, output2);
  fOps->transformFromXtraSizeArray(output2, output2);
  double diff3 = 0;
  for (I=0; I<VL; I++) {
    diff3 += sqr(output[I].x - output2[I].x) + sqr(output[I].y - output2[I].y);
  }
  diff3 = sqrt(diff3);
  
  
  randomInput();
  randomOutput();
  for (I=0; I<VLXtr; I++) output2[I] = output[I];
  fOps->transformFromXtraSizeArray(input, input);
  fOps->transformFromXtraSizeArray(output, output);
  for (I=0; I<VL; I++) {
    output[I].x = alpha.x*(input[I].x - output[I].x);
    output[I].y = alpha.x*(input[I].y - output[I].y);
  }
  fOps->transformToXtraSizeArray(input, input);
  SSE_ComplexVectorAdditionSPECIAL2(L0,L1,L2,L3,xtraSize1,xtraSize2,xtraSize3, alpha.x, input, output2);
  fOps->transformFromXtraSizeArray(output2, output2);
  double diff4 = 0;
  for (I=0; I<VL; I++) {
    diff4 += sqr(output[I].x - output2[I].x) + sqr(output[I].y - output2[I].y);
  }
  diff4 = sqrt(diff4);

#ifdef useBLAS
  startTimer();
  count = 0;
  while (timePassed()<10) {
    cblas_zaxpy(VL, &alpha, input, 1, output, 1);
    count++;
  }
  printf("BLAS: %d iterations (vector add) within %1.2f seconds.\n",count,timePassed());
#endif 

  startTimer();
  count = 0;
  while (timePassed()<10) {
//SSE_ComplexVectorAddition(L0,L1,L2,L3,xtraSize1,xtraSize2,xtraSize3, alpha, input, output2);
    xSSEObj->xSSE_ComplexVectorAddition_Wrapper(input, output2, alpha, 8, L0,L1,L2,L3,xtraSize1,xtraSize2,xtraSize3);      
    count++;
  }
  printf("SSE-Improved: %d iterations (vector add) within %1.2f seconds.\n",count,timePassed());
  printf("Difference between results = %1.15f \n",diff);


#ifdef useBLAS
  startTimer();
  count = 0;
  Complex res;
  while (timePassed()<10) {
    cblas_zaxpy(VL, &alpha, input, 1, output, 1);
    cblas_zdotc_sub(VL, output, 1, output, 1, &res);
    count++;
  }
  printf("BLAS: %d iterations (vector add with scalar product) within %1.2f seconds.\n",count,timePassed());
#endif

  startTimer();
  count = 0;
  double n3;
  while (timePassed()<10) {
    SSE_ComplexVectorAdditionWithSquaredNorm(L0,L1,L2,L3,xtraSize1,xtraSize2,xtraSize3, alpha, input, output2, n3);
    count++;
  }
  printf("SSE-Improved: %d iterations (vector add with scalar product) within %1.2f seconds.\n",count,timePassed());
#ifdef useBLAS
  printf("Difference between results = %1.15f  %1.15f \n",diff2, fabs((double)(n2-n1)));
#else
  printf("Difference between results = %1.15f \n",diff2);
#endif  
  startTimer();
  count = 0;
  while (timePassed()<10) {
    for (I=0; I<VL; I++) {
      output[I].x = alpha.x*output[I].x + input[I].x;
      output[I].y = alpha.x*output[I].y + input[I].y;
    }
    count++;
  }
  printf("C++: %d iterations (vector add SPECIAL1) within %1.2f seconds.\n",count,timePassed());

  startTimer();
  count = 0;
  while (timePassed()<10) {
    SSE_ComplexVectorAdditionSPECIAL1(L0,L1,L2,L3,xtraSize1,xtraSize2,xtraSize3, alpha.x, input, output2);
    count++;
  }
  printf("SSE-Improved: %d iterations (vector add SPECIAL1) within %1.2f seconds.\n",count,timePassed());
  printf("Difference between results = %1.15f\n",diff3); 
  
  startTimer();
  count = 0;
  while (timePassed()<10) {
    for (I=0; I<VL; I++) {
      output[I].x = alpha.x*(input[I].x - output[I].x);
      output[I].y = alpha.x*(input[I].y - output[I].y);
    }
    count++;
  }
  printf("C++: %d iterations (vector add SPECIAL2) within %1.2f seconds.\n",count,timePassed());

  startTimer();
  count = 0;
  while (timePassed()<10) {
    SSE_ComplexVectorAdditionSPECIAL2(L0,L1,L2,L3,xtraSize1,xtraSize2,xtraSize3, alpha.x, input, output2);
    count++;
  }
  printf("SSE-Improved: %d iterations (vector add SPECIAL2) within %1.2f seconds.\n",count,timePassed());
  printf("Difference between results = %1.15f\n",diff4); 
  delete xSSEObj;
}
 
 
void testSSEDiracOperatorApplication() {
  printf("\nTesting SSE-Dirac-Operator (in Fourier-Space)...\n");
  randomInput();
  randomOutput();

  int VL = fOps->getVectorLength();
  int I;
  fOps->transformFromXtraSizeArray(input, interim);  
  executeNeubergerDiracMultiplicationInFourierSpace(L0, L1, L2, L3, interim, output, fOps->NeubergerDiracOpApplicationSinPStdData, fOps->NeubergerDiracOpApplicationEWData);
  SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace(L0, L1, L2, L3, input, output2, fOps->NeubergerDiracOpApplicationSinPStdData, fOps->NeubergerDiracOpApplicationEWData);
  fOps->transformFromXtraSizeArray(output2, output2);  
  double diff = 0;
  for (I=0; I<VL; I++) {
    diff += sqr(output[I].x - output2[I].x) + sqr(output[I].y - output2[I].y);
  }
  diff = sqrt(diff);

  startTimer();
  int count = 0;
  while (timePassed()<10) {
    executeNeubergerDiracMultiplicationInFourierSpace(L0, L1, L2, L3, input, output, fOps->NeubergerDiracOpApplicationSinPStdData, fOps->NeubergerDiracOpApplicationEWData);
    count++;
  }
  printf("C++: %d iterations within %1.2f seconds.\n",count,timePassed());

  startTimer();
  count = 0;
  while (timePassed()<10) {
    SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace(L0, L1, L2, L3, interim, output2, fOps->NeubergerDiracOpApplicationSinPStdData, fOps->NeubergerDiracOpApplicationEWData);
    count++;
  }
  printf("SSE-Improved: %d iterations within %1.2f seconds.\n",count,timePassed());
  printf("Difference between results = %1.15f \n",diff);
} 


void testSSEBOperatorApplication(double yN, double rho) {
  printf("\nTesting B-Operator with yN=%1.2f and rho=%1.2f...\n",yN, rho);
  randomInput();
  randomPhi(0);
  

  int VL = fOps->getVectorLength();
  int I;
  double twoRho = 2.0 * rho;
  fOps->transformFromXtraSizeArray(input, input);  
  performf_YB_2rho_AndCopyToOutput(L0,L1,L2,L3,xtraSize1,xtraSize2,xtraSize3, yN, twoRho, input, phiField, interim, output);
  fOps->transformToXtraSizeArray(input, input);    
  SSE_Performf_YB_2rho_AndCopyToOutput(L0,L1,L2,L3,xtraSize1,xtraSize2,xtraSize3, yN, twoRho, input, phiField, interim2, output2);
  fOps->transformFromXtraSizeArray(interim2, interim2);  
  fOps->transformFromXtraSizeArray(output2, output2);  
  double diff = 0;
  double diff2 = 0;
  for (I=0; I<VL; I++) {
    diff += sqr(output[I].x - output2[I].x) + sqr(output[I].y - output2[I].y);
    diff2 += sqr(interim[I].x - interim2[I].x) + sqr(interim[I].y - interim2[I].y);
  }
  diff = sqrt(diff);
  diff2 = sqrt(diff2);

  startTimer();
  int count = 0;
  while (timePassed()<10) {
    performf_YB_2rho_AndCopyToOutput(L0,L1,L2,L3,xtraSize1,xtraSize2,xtraSize3, yN, twoRho, input, phiField, interim, output);
    count++;
  }
  printf("C++: %d iterations within %1.2f seconds.\n",count,timePassed());

  startTimer();
  count = 0;
  while (timePassed()<10) {
    SSE_Performf_YB_2rho_AndCopyToOutput(L0,L1,L2,L3,xtraSize1,xtraSize2,xtraSize3, yN, twoRho, input, phiField, interim, output);
    count++;
  }
  printf("SSE-Improved: %d iterations within %1.2f seconds.\n",count,timePassed());
  printf("Difference between results = %1.15f, and diff2 = %1.15f \n",diff, diff2);
} 


void testMMDaggerxQuasiHermiteanNeubergerWithChi(double yN, bool useR, bool inFourierSpace) {
  printf("\nTesting MMdag in quasi-hermitean mode with yN = %1.2f...\n",yN);  

  randomPhi(0);
  randomInput();
  
  ComplexVector inputVec(fOps->getVectorLength(),input);
  ComplexVector outputVec(fOps->getVectorLength(),output);
  ComplexVector outputVec2(fOps->getVectorLength(),output2);
  
  fOps->setYukawaCoupling(yN);
  fOps->setRPreconditioner(true, 1.0, 1.0);
  double rho,r;
  fOps->getDiracParameters(rho, r);
    
  startTimer();
  int count = 0;
  while (timePassed()<10) {
    fOps->executeFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(input, output, phiField, useR, inFourierSpace);   
    count++;
  }
  printf("%d iterations within %1.2f seconds.\n",count,timePassed());
    
  fOps->executeFermionQuasiHermiteanMatrixMultiplication(input, output2, phiField, useR, true, inFourierSpace);
  fOps->executeFermionQuasiHermiteanMatrixMultiplication(output2, output2, phiField, useR, false, inFourierSpace);
  fOps->transformFromXtraSizeArray(output, output);
  fOps->transformFromXtraSizeArray(output2, output2);  
  ComplexVector diffVec = outputVec - outputVec2;
  printf("Comparison with two applications of single M: Difference Norm: %1.15f\n",diffVec.getNorm());
  
  fOps->executeDiracMatrixMultiplication(input, output, false);
  Complex alpha(-2*rho ,0);
  SSE_ComplexVectorAddition(L0,L1,L2,L3,xtraSize1,xtraSize2,xtraSize3, alpha, input, output);
  fOps->executeFermionQuasiHermiteanMatrixMultiplication(output, output2, phiField, false, false, false);     

  fOps->executeFermionMatrixMultiplication(input, output, phiField, false, NULL, NULL, 0, 0);
  fOps->executePiModeRemoverOperator(output, output, false);
  fOps->transformFromXtraSizeArray(output, output);
  fOps->transformFromXtraSizeArray(output2, output2);  
  ComplexVector diffVec2 = outputVec - outputVec2;
  printf("Comparison with non-hermitean MMdag-Matrix: Difference Norm: %1.15f\n",diffVec2.getNorm());
  
}


void testMxNeubergerWithChi(double yN, double massSplit, double explicitMass, bool compare, bool dag) {
  if (dag) {
    printf("\nTesting daggered Dirac-Operator (Neuberger with Xi) with yN = %1.2f...\n",yN);
  } else {
    printf("\nTesting Dirac-Operator (Neuberger with Xi) with yN = %1.2f...\n",yN);  
  }
  randomPhi(0);
  randomInput();
  
  ComplexVector inputVec(fOps->getVectorLength(),input);
  ComplexVector outputVec(fOps->getVectorLength(),output);
  ComplexVector outputVec2(fOps->getVectorLength(),output2);
  
  fOps->setYukawaCoupling(yN);
  fOps->setMassSplitRatio(massSplit);
  fOps->setExplicitMass(explicitMass);

  startTimer();
  int count = 0;
  while (timePassed()<10) {
    if (dag) {
      fOps->executeFermionDaggerMatrixMultiplication(input, output, phiField, 0);   
    } else {
      fOps->executeFermionMatrixMultiplication(input, output, phiField, false, NULL, NULL, 0, 0);   
    }
    count++;
  }
  printf("%d iterations within %1.2f seconds.\n",count,timePassed());

  if (compare) {
    fOps->transformFromXtraSizeArray(input, input);
    fOps->transformFromXtraSizeArray(output, output);
    ComplexMatrix mat(1);
    fOps->constructNeubergerWithXiFermionMatrix(mat, (vector4D*) phiField, massSplit);
    if (dag) mat.dagger();
    outputVec2 = mat * inputVec;
    ComplexVector diffVec = outputVec - outputVec2;
    printf("Comparison with full-matrix Mul: Difference Norm: %1.15f\n",diffVec.getNorm());
  }
}


void testMMDaggerxNeubergerWithChi(double yN, double massSplit, double explicitMass, bool compare) {
  printf("\nTesting MMdag-Operator (Neuberger with Xi) with yN = %1.2f...\n",yN);
  randomPhi(0);
  randomInput();

  ComplexVector inputVec(fOps->getVectorLength(),input);
  ComplexVector outputVec(fOps->getVectorLength(),output);
  ComplexVector outputVec2(fOps->getVectorLength(),output2);

  fOps->setYukawaCoupling(yN);
  fOps->setMassSplitRatio(massSplit);
  fOps->setExplicitMass(explicitMass);

  startTimer();
  int count = 0;
  while (timePassed()<10) {
    fOps->executeFermionDaggerMatrixMultiplication(input, output, phiField, 0);   
    fOps->executeFermionMatrixMultiplication(output, output, phiField, false, NULL, NULL, 0, 0);   
    count++;
  }
  printf("%d iterations within %1.2f seconds.\n",count,timePassed());

  if (compare) {
    fOps->transformFromXtraSizeArray(input, input);    
    fOps->transformFromXtraSizeArray(output, output);    
    ComplexMatrix mat(1);
    ComplexMatrix mat2(1);
    fOps->constructNeubergerWithXiFermionMatrix(mat, (vector4D*) phiField, massSplit);
    fOps->constructNeubergerWithXiFermionMatrix(mat2, (vector4D*) phiField, massSplit);
    mat2.dagger();
    ComplexMatrix mat3 = mat*mat2;
    outputVec2 = mat3 * inputVec;
    ComplexVector diffVec = outputVec - outputVec2;
    printf("Comparison with full-matrix Mul: Difference Norm: %1.15f\n",diffVec.getNorm());
  }
}


void testCompactMMDaggerxNeubergerWithChi(double yN, bool inFourierSpace, bool precUse, double precM, double precS, bool QprecUse, double QprecMu, double QprecBeta) {
  printf("\nTesting compact MMdag-Operator (Neuberger with Xi) with yN = %1.2f...\n",yN);
  randomPhi(0);
  randomInput();

  ComplexVector inputVec(fOps->getVectorLength(),input);
  ComplexVector outputVec(fOps->getVectorLength(),output);
  ComplexVector outputVec2(fOps->getVectorLength(),output2);

  fOps->setYukawaCoupling(yN);
  fOps->setPreconditioner(precUse, precM, precS);
  fOps->setQPreconditioner(QprecUse, QprecMu, QprecBeta);

  startTimer();
  int count = 0;
  while (timePassed()<10) {
    fOps->executeFermionDaggerMatrixMultiplication(input, output, phiField, 1);   
    fOps->executeFermionMatrixMultiplication(output, output, phiField, false, NULL, NULL, 2, 1);   
    count++;
  }
  printf("For single M-Applications: %d iterations within %1.2f seconds.\n",count,timePassed());

  if (inFourierSpace) {
    fOps->performFFT(input,output2,ExtremeFFT4D_Forward);
    int VLxtrSize = fOps->getVectorLengthXtrSize();
    double VolumeNorm = 1.0 / (L0*L1*L2*L3);
    for (int I=0; I<VLxtrSize; I++) {
      input[I].x = output2[I].x * VolumeNorm;
      input[I].y = output2[I].y * VolumeNorm;
    }
  }
  startTimer();
  count = 0;
  while (timePassed()<10) {
    fOps->executeFermionMatrixFermionDaggerMatrixMultiplication(input, output2, phiField, 2, 1, inFourierSpace);
    count++;
  }
  printf("For compact application: %d iterations within %1.2f seconds.\n",count,timePassed());
  if (inFourierSpace) {
    fOps->performFFT(output2,input,ExtremeFFT4D_Backward);
    int VLxtrSize = fOps->getVectorLengthXtrSize();
    for (int I=0; I<VLxtrSize; I++) {
      output2[I].x = input[I].x;
      output2[I].y = input[I].y;
    }
  }

  ComplexVector diffVec = outputVec - outputVec2;
  printf("Comparison: Difference Norm: %1.15f\n",diffVec.getNorm());
}


void testPreconditioner(double yN, double m, double s) {
  printf("\nTesting Preconditioners with m = %1.2f and s = %1.2f...\n",m,s);
  randomPhi(m,s,0);
  randomInput();
  
  ComplexVector inputVec(fOps->getVectorLength(),input);
  ComplexVector outputVec(fOps->getVectorLength(),output);
  ComplexVector outputVec2(fOps->getVectorLength(),output2);

  fOps->setYukawaCoupling(yN);
  fOps->setPreconditioner(true, m, s);


  startTimer();
  int count = 0;  
  while (timePassed()<10) {
    fOps->executeFermionMatrixStaticInverseMultiplication(input, output, true, true);
    count++;
  }
  printf("Double static inverse MdagM in Fourier Space: %d iterations within %1.2f seconds.\n",count,timePassed());

  startTimer();
  count = 0;
  while (timePassed()<10) {
    fOps->executeFermionMatrixStaticInverseMultiplication(input, output, false, true);
    count++;
  }
  printf("Single static inverse M in Fourier Space: %d iterations within %1.2f seconds.\n",count,timePassed());
  
  fOps->executeFermionMatrixStaticInverseMultiplication(input, output, false, false);
  fOps->executeFermionMatrixMultiplication(output, output2, phiField, false, NULL, NULL, 0, 0);
  fOps->transformFromXtraSizeArray(output2, output2);
  fOps->transformFromXtraSizeArray(input, input);
  ComplexVector diffVec = inputVec - outputVec2;
  fOps->transformToXtraSizeArray(input, input);  
  printf("Single static inverse M: Comparison with full-matrix Mul: Difference Norm: %1.15f\n",diffVec.getNorm());


  fOps->executeFermionDaggerMatrixMultiplication(input, output, phiField, 0);   
  fOps->executeFermionMatrixStaticInverseMultiplication(output, output, true, false);
  fOps->executeFermionMatrixMultiplication(output, output2, phiField, false, NULL, NULL, 0, 0);
  fOps->transformFromXtraSizeArray(output2, output2);
  fOps->transformFromXtraSizeArray(input, input);
  diffVec = inputVec - outputVec2;
  fOps->transformToXtraSizeArray(input, input);  
  printf("Double static inverse MdagM: Comparison with full-matrix Mul: Difference Norm: %1.15f\n",diffVec.getNorm());


  startTimer();
  count = 0;
  while (timePassed()<10) {
    fOps->executeFermionMatrixMultiplication(input, output, phiField, false, NULL, NULL, 1, 0);
    count++;
  }
  printf("M with Preconditioner: %d iterations within %1.2f seconds.\n",count,timePassed());
  fOps->transformFromXtraSizeArray(output, output);
  fOps->transformFromXtraSizeArray(input, input);
  diffVec = inputVec - outputVec;
  fOps->transformToXtraSizeArray(input, input);  
  printf("Comparison with full-matrix Mul: Difference Norm: %1.15f\n",diffVec.getNorm());

  startTimer();
  count = 0;
  while (timePassed()<10) {
    fOps->executeFermionDaggerMatrixMultiplication(input, output, phiField, 0);   
    fOps->executeFermionMatrixMultiplication(output, output, phiField, false, NULL, NULL, 2, 0);
    count++;
  }
  printf("MMdag with Preconditioner: %d iterations within %1.2f seconds.\n",count,timePassed());
  fOps->transformFromXtraSizeArray(output, output);
  fOps->transformFromXtraSizeArray(input, input);
  diffVec = inputVec - outputVec;
  fOps->transformToXtraSizeArray(input, input);  
  printf("Comparison with full-matrix Mul: Difference Norm: %1.15f\n",diffVec.getNorm());
  
  fOps->setPreconditioner(false, m, s);
}


void testSolverNeubergerWithChi(double yN, bool doubleM, bool QHMode, double TOL) {
  if (doubleM) {
    printf("\nTesting MMdag-Solver (Neuberger with Xi) with yN = %1.2f, QHMode=%d, TOL = %1.15f ...\n",yN,QHMode,TOL);
  } else {
    printf("\nTesting M-Solver (Neuberger with Xi) with yN = %1.2f, QHMode=%d, TOL = %1.15f ...\n",yN,QHMode,TOL);  
  }
  randomPhi(0);
  randomInput();
  
  ComplexVector inputVec(fOps->getVectorLength(),input);
  ComplexVector outputVec(fOps->getVectorLength(),output);
  ComplexVector outputVec2(fOps->getVectorLength(),output2);
  
  fOps->setYukawaCoupling(yN);

  if (doubleM) {
    if (QHMode) {
      fOps->executeFermionQuasiHermiteanMatrixMultiplication(input, output, phiField, false, true, false);
      fOps->executeFermionQuasiHermiteanMatrixMultiplication(output, output, phiField, false, false, false);
    } else {
      fOps->executeFermionDaggerMatrixMultiplication(input, output, phiField, 0);   
      fOps->executeFermionMatrixMultiplication(output, output, phiField, false, NULL, NULL, 0, 0);       
    }
  } else {
    if (QHMode) {
      fOps->executeFermionQuasiHermiteanMatrixMultiplication(input, output, phiField, false, false, false);    
    } else {
      fOps->executeFermionMatrixMultiplication(input, output, phiField, false, NULL, NULL, 0, 0); 
    }
  }

  startTimer();
  int count = 0;
  while (timePassed()<10) {
    int neededIter = 0;
    fOps->solveFermionMatrixLGS(output, output2, phiField, TOL, doubleM, QHMode, -1, neededIter);
    count++;
  }
  printf("%d iterations within %1.2f seconds.\n",count,timePassed());


  fOps->transformFromXtraSizeArray(output2, output2);
  fOps->transformFromXtraSizeArray(input, input);
  ComplexVector diffVec = inputVec - outputVec2;
  printf("Comparison with exact solution: Absolute difference Norm: %1.15f, Relative difference Norm: %1.15f\n",diffVec.getNorm(),diffVec.getNorm()/inputVec.getNorm());
  
  int I0, I1, I2, I3, I;
  double fac;
  double m = 0;
  double s = 0;
  count = 0;
  vector4D sum1; 
  sum1[0] = sum1[1] = sum1[2] = sum1[3] = 0;
  vector4D sum2;
  sum2[0] = sum2[1] = sum2[2] = sum2[3] = 0;
  for (I0 = 0; I0<L0; I0++) {
    for (I1 = 0; I1<L1; I1++) {
      for (I2 = 0; I2<L2; I2++) {
        for (I3 = 0; I3<L3; I3++) {
	  fac = 1;
	  if (((I0+I1+I2+I3) % 2) == 1) fac = -1;
	  for (I=0; I<4; I++) {
	    sum1[I] += phiField[count];
	    sum2[I] += fac*phiField[count];
	    count++;
	  }
	}
      }
    }
  }
  fOps->transformToXtraSizeArray(input, input);
  m = sqrt(sqr(sum1[0])+sqr(sum1[1])+sqr(sum1[2])+sqr(sum1[3])) / (L0*L1*L2*L3);
  s = sqrt(sqr(sum2[0])+sqr(sum2[1])+sqr(sum2[2])+sqr(sum2[3])) / (L0*L1*L2*L3);
  
  printf("Redoing Inversion with P-/R-preconditioner. Measured m = %1.3f, s = %1.3f\n",m,s);
  fOps->setPreconditioner(true, m, s);
  fOps->setRPreconditioner(true, m, 1.0);
  
  if (doubleM) {
    if (QHMode) {
      fOps->executeFermionQuasiHermiteanMatrixMultiplication(input, output, phiField, true, true, false);
      fOps->executeFermionQuasiHermiteanMatrixMultiplication(output, output, phiField, true, false, false);    
    } else {
      fOps->executeFermionDaggerMatrixMultiplication(input, output, phiField, 0);   
      fOps->executeFermionMatrixMultiplication(output, output, phiField, false, NULL, NULL, 2, 0);    
    }
  } else {
    if (QHMode) {
      fOps->executeFermionQuasiHermiteanMatrixMultiplication(input, output, phiField, true, false, false);        
    } else {
      fOps->executeFermionMatrixMultiplication(input, output, phiField, false, NULL, NULL, 1, 0); 
    }
  }
  
  startTimer();
  count = 0;
  while (timePassed()<10) {
    int neededIter = 0;
    fOps->solveFermionMatrixLGS(output, output2, phiField, TOL, doubleM, QHMode, -1, neededIter);
    count++;
  }
  printf("%d iterations within %1.2f seconds.\n",count,timePassed());

  fOps->transformFromXtraSizeArray(output2, output2);
  fOps->transformFromXtraSizeArray(input, input);
  diffVec = inputVec - outputVec2;
  fOps->transformToXtraSizeArray(input, input);  
  printf("Comparison with exact solution: Absolute difference Norm: %1.15f, Relative difference Norm: %1.15f\n",diffVec.getNorm(),diffVec.getNorm()/inputVec.getNorm());

  if (QHMode) {
    fOps->setPreconditioner(false, m, s);
    fOps->setQPreconditioner(false, 0.5, 0.5);
    fOps->setRPreconditioner(false, 1.0, 1.0);  
    return;
  }

  printf("Redoing Inversion with Q-preconditioner. \n");
  fOps->setQPreconditioner(true, 0.5, 0.5);
  
  if (doubleM) {
    fOps->executeFermionDaggerMatrixMultiplication(input, output, phiField, 1);   
    fOps->executeFermionMatrixMultiplication(output, output, phiField, false, NULL, NULL, 2, 1);
  } else {
    fOps->executeFermionMatrixMultiplication(input, output, phiField, false, NULL, NULL, 1, 1); 
  }
  
  startTimer();
  count = 0;
  while (timePassed()<10) {
    int neededIter = 0;
    fOps->solveFermionMatrixLGS(output, output2, phiField, TOL, doubleM, QHMode, -1, neededIter);
    count++;
  }
  printf("%d iterations within %1.2f seconds.\n",count,timePassed());

  fOps->transformFromXtraSizeArray(output2, output2);
  fOps->transformFromXtraSizeArray(input, input);
  diffVec = inputVec - outputVec2;
  fOps->transformToXtraSizeArray(input, input);  
  printf("Comparison with exact solution: Absolute difference Norm: %1.15f, Relative difference Norm: %1.15f\n",diffVec.getNorm(),diffVec.getNorm()/inputVec.getNorm());

  fOps->setPreconditioner(false, m, s);
  fOps->setQPreconditioner(false, 0.5, 0.5);
  fOps->setRPreconditioner(false, 1.0, 1.0);
}


double matrixFunction(double x) {
  return 1.0 / sqrt(x);
}


void testMMDaggerSQTRSolverNeubergerWithChi(double yN, double TOL, bool precUse, double precM, double precS, bool QprecUse, double QprecMu, double QprecBeta, bool RprecUse, double RprecM, double RprecF, int compareWithAuxVecCount) {
  printf("\nTesting MMdag-SQRT-Solver (Neuberger with Xi) with yN = %1.2f, TOL = %1.15f ...\n",yN,TOL);
  randomPhi(0);
  randomInput();
  
  ComplexVector inputVec(fOps->getVectorLength(),input);
  ComplexVector outputVec(fOps->getVectorLength(),output);
  ComplexVector outputVec2(fOps->getVectorLength(),output2);
  
  fOps->setYukawaCoupling(yN);
  fOps->setPreconditioner(precUse, precM, precS);
  fOps->setQPreconditioner(QprecUse, QprecMu, QprecBeta);
  fOps->setRPreconditioner(RprecUse, RprecM, RprecF);
  

//  fOps->executeFermionMatrixFermionDaggerMatrixMultiplication(input, output, phiField, 2, 1, true);
//  fOps->executeFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(input, output, phiField, true, true);



  int neededIter = 0;
  fOps->applyFermionMatrixMMDaggerFunction(input, output, phiField, &matrixFunction, TOL, 0, neededIter, 0, NULL, true, true);
  startTimer();
  int count = 0;
  while (timePassed()<10) {
    int neededIter = 0;
    fOps->applyFermionMatrixMMDaggerFunction(input, output, phiField, &matrixFunction, TOL, 0, neededIter, 0, NULL, true, true);
    count++;
  }
  printf("%d iterations within %1.2f seconds.\n",count,timePassed());
  neededIter = 0;
  fOps->applyFermionMatrixMMDaggerFunction(output, output, phiField, &matrixFunction, TOL, 0, neededIter, 0, NULL, true, true);
  
  fOps->executeFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(output, output2, phiField, true, true);

  fOps->transformFromXtraSizeArray(input, input);
  fOps->transformFromXtraSizeArray(output2, output2);
  ComplexVector diffVec = outputVec2 - inputVec;
  printf("Comparison with exact solution: Absolute difference Norm: %1.15f, Relative difference Norm: %1.15f\n",diffVec.getNorm(),diffVec.getNorm()/inputVec.getNorm());


  Complex** auxVecs= new Complex*[compareWithAuxVecCount];
  for (int I=0; I<compareWithAuxVecCount; I++) {
    auxVecs[I] = fOps->createFermionVector();
  }
  
  neededIter = 0;
  fOps->transformToXtraSizeArray(input, input);
  fOps->applyFermionMatrixMMDaggerFunction(input, output, phiField, &matrixFunction, TOL, 0, neededIter, compareWithAuxVecCount, auxVecs, true, true);
  startTimer();
  count = 0;
  while (timePassed()<10) {
    int neededIter = 0;
    fOps->applyFermionMatrixMMDaggerFunction(input, output, phiField, &matrixFunction, TOL, 0, neededIter, compareWithAuxVecCount, auxVecs, true, true);
    count++;
  }
  printf("%d iterations within %1.2f seconds.\n",count,timePassed());
  neededIter = 0;
  fOps->applyFermionMatrixMMDaggerFunction(output, output, phiField, &matrixFunction, TOL, 0, neededIter, compareWithAuxVecCount, auxVecs, true, true);

  fOps->executeFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(output, output2, phiField, true, true);

  for (int I=0; I<compareWithAuxVecCount; I++) {
    fOps->destroyFermionVector(auxVecs[I]);
  }
  delete[] auxVecs;

  fOps->transformFromXtraSizeArray(input, input);
  fOps->transformFromXtraSizeArray(output2, output2);
  ComplexVector diffVec2 = outputVec2 - inputVec;
  printf("Comparison with solution from auxVectors (%d): Absolute difference Norm: %1.15f, Relative difference Norm: %1.15f\n",compareWithAuxVecCount,diffVec2.getNorm(),diffVec2.getNorm()/inputVec.getNorm());
}


double matrixFunction2_Parameter = 0;
double matrixFunction2(double x) {
  return 1.0/sqrt(x+matrixFunction2_Parameter);
}


void testMMDaggerChebyshevPolynomialApplication(double yN, double TOL, bool precUse, double precM, double precS, bool QprecUse, double QprecMu, double QprecBeta) {
  printf("\nTesting MMdag-Chebyshev-Application (Neuberger with Xi) with yN = %1.2f, TOL = %1.15f ...\n",yN,TOL);
  randomPhi(0);
  randomInput();
  
  ComplexVector inputVec(fOps->getVectorLength(),input);
  ComplexVector outputVec(fOps->getVectorLength(),output);
  ComplexVector outputVec2(fOps->getVectorLength(),output2);
  
  fOps->setYukawaCoupling(yN);
  fOps->setPreconditioner(precUse, precM, precS);
  fOps->setQPreconditioner(QprecUse, QprecMu, QprecBeta);
  fOps->setRPreconditioner(true, 1.0, 1.0);

  matrixFunction2_Parameter = 1.0;
  GeneralChebyshevApproximation chebyTest; 
  chebyTest.calcApproximation(&matrixFunction2, 0.0, 2.0, TOL, 10000);

//  fOps->executeFermionMatrixFermionDaggerMatrixMultiplication(input, output, phiField, 2, 1, true);
  fOps->executeFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(input, output, phiField, true, true);
  Complex compFac(matrixFunction2_Parameter,0);
  SSE_ComplexVectorAddition(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, compFac, input, output); 
  
  startTimer();
  int count = 0;
  while (timePassed()<10) {
    fOps->applyFermionMatrixMMDaggerChebyshevPolynomial(output, output2, phiField, &chebyTest, true, true);
    count++;
  }
  printf("%d iterations within %1.2f seconds.\n",count,timePassed());
  fOps->applyFermionMatrixMMDaggerChebyshevPolynomial(output2, output2, phiField, &chebyTest, true, true);


  fOps->transformFromXtraSizeArray(input, input);
  fOps->transformFromXtraSizeArray(output2, output2);
  ComplexVector diffVec = outputVec2 - inputVec;
  printf("Comparison with exact solution: Absolute difference Norm: %1.15f, Relative difference Norm: %1.15f\n",diffVec.getNorm(),diffVec.getNorm()/inputVec.getNorm());
}


void testdSdPhiNeubergerWithChi(double yN, double h, double TOL) {
  printf("\nTesting dS/dPhi (Neuberger with Xi) with yN = %1.2f, h = %f, TOL = %1.15f ...\n",yN,h,TOL);
  randomPhi(0);
  randomInput();
  
  ComplexVector inputVec(fOps->getVectorLength(),input);
  ComplexVector outputVec(fOps->getVectorLength(),output);
  ComplexVector outputVec2(fOps->getVectorLength(),output2);
  
  fOps->setYukawaCoupling(yN);

  int testIndex = ((int)(zufall()*L0*L1*L2*L3)) % (L0*L1*L2*L3);
  int I;
  double dSdPhi_Comp[4];
  int neededIter = 0;
  fOps->solveFermionMatrixLGS(input, output, phiField, TOL, true, false, -1, neededIter);
  Complex dummy;
  SSE_ComplexScalarProduct(L0,L1,L2,L3,xtraSize1,xtraSize2,xtraSize3, input, output, dummy);
  double SOld = dummy.x;  
  for (I=0; I<4; I++) {
    phiField[4*testIndex+I] += h;
    fOps->solveFermionMatrixLGS(input, output, phiField, TOL, true, false, -1, neededIter);  
    phiField[4*testIndex+I] -= h;
    SSE_ComplexScalarProduct(L0,L1,L2,L3,xtraSize1,xtraSize2,xtraSize3, input, output, dummy);

    double SNew = dummy.x;  
    dSdPhi_Comp[I] = (SNew-SOld) / h;
  }

  startTimer();
  int count = 0;
  while (timePassed()<10) {
    fOps->executeMultiplicationVectorWithDerivativesOfMMdaggerInverse(input, output, dSdPhi, (double*) phiField, TOL);
    count++;
  }
  printf("%d iterations within %1.2f seconds.\n",count,timePassed());

  printf("Comparison with numerical solution:\n");
  for (I=0; I<4; I++) {  
    printf("  Relative difference Norm for Phi_%d: %1.15f\n",I,(dSdPhi_Comp[I]-dSdPhi[4*testIndex+I]) / (dSdPhi_Comp[I]));
  }
}


void testConditionNumberEstimation(double yN) {
  printf("\nTesting Condition Number Estimation (Neuberger with Xi) with yN = %1.2f ...\n",yN);
  randomPhi(0.5,5,0);
  fOps->setYukawaCoupling(yN);
  double eigMin, eigMax, cond;
//  fOps->estimateFermionMatrixConditionNumber(phiField,eigMin,eigMax,cond);
  printf("Estimated Min. EV %f, Max. EV. %f, Condition-Number: %f\n",eigMin,eigMax,cond);
  int I;

//  fOps->exactFermionMatrixConditionNumber(phiField,min,max,cond);
  
  ComplexMatrix mat(1);
  ComplexMatrix mat2(1);
  fOps->constructNeubergerWithXiFermionMatrix(mat, (vector4D*) phiField);
  fOps->constructNeubergerWithXiFermionMatrix(mat2, (vector4D*) phiField);
  mat2.dagger();
  ComplexMatrix mat3 = mat * mat2;  
  
  printEigenvalues("data/MatrixMMdag", mat3, Complex(1,0));
  eigMin = sqrt(sqr(mat3.eigenvalues[0].x) + sqr(mat3.eigenvalues[0].y));
  eigMax = eigMin;  
  for (I=0; I<mat3.matrixSize; I++) {
    double n = sqrt(sqr(mat3.eigenvalues[I].x) + sqr(mat3.eigenvalues[I].y));
    if (n>eigMax) eigMax = n;
    if (n<eigMin) eigMin = n;
  }
  cond = eigMin/eigMax;  
  printf("Real Min. EV. %f, Max. EV. %f, Condition-Number: %f\n",eigMin, eigMax, cond);
}


void testEigenvaluesNeubergerWithChi(double yN) {
  printf("\nTesting Eigenvalues (Neuberger with Xi) with yN = %1.2f ...\n",yN);
  randomPhi(0);
  fOps->setYukawaCoupling(yN);

  ComplexMatrix mat(1);
  fOps->constructNeubergerWithXiFermionMatrix(mat, (vector4D*) phiField);
  printEigenvalues("data/MatrixM3", mat, Complex(-0.5,0));
}


void testHMCPropagation(double yN, double kappa, double lambda) {
  printf("\nTesting Integrators with yN = %1.2f, kappa = %f, lambda = %1.15f...\n",yN, kappa, lambda);
  double S1, S2;
  randomPhi(0);
  randomInput();
  
  HMCProp->setKappa(kappa);
  HMCProp->setLambda(lambda);
  fOps->setYukawaCoupling(yN);
 
  HMCProp->samplePhiMomentumField();
  HMCProp->sampleOmegaFields();

  int iter = 10;
  while (iter<=280) {
    S1 = HMCProp->calcTotalAction(1E-12);
    HMCProp->LeapFrogMarkovStep(iter,  1.0 / iter, 1E-10, 1E-12);
    S2 = HMCProp->calcTotalAction(1E-12);
    printf("Leap-Frog integration: With %d Inversions achieved dS = %1.15f (=%1.15f if accepted)\n", iter, (HMCProp->SafterProp-HMCProp->SbeforeProp), S2-S1);
    HMCProp->restorePhiFields();
    iter *= 2;
  }
  printf("\n");
/*  iter = 5;
  while (2*iter<=280) {
    S1 = HMCProp->calcTotalAction(1E-12);
    HMCProp->OmelyanO2MarkovStepn(iter,  1.0 / iter, 1E-10, 1E-12);
    S2 = HMCProp->calcTotalAction(1E-12);
    printf("Omelyan Order 2 integration: With %d Inversions achieved dS = %1.15f (=%1.15f if accepted)\n", 2*iter, (HMCProp->SafterProp-HMCProp->SbeforeProp), S2-S1);
    HMCProp->restorePhiFields();
    iter *= 2;
  }
  printf("\n");*/
  iter = 2;
  while (5*iter<=280) {
    S1 = HMCProp->calcTotalAction(1E-12);
    HMCProp->OmelyanO4MarkovStep(iter,  1.0 / iter, 1E-10, 1E-12);
    S2 = HMCProp->calcTotalAction(1E-12);
    printf("Omelyan Order 4 integration: With %d Inversions achieved dS = %1.15f (=%1.15f if accepted)\n", 5*iter, (HMCProp->SafterProp-HMCProp->SbeforeProp), S2-S1);
    HMCProp->restorePhiFields();
    iter *= 2;
  }
}


/*void testFourierTrafo() {
  printf("\nTesting Fast-Fourier-Transformation ...\n");
  int VL = fOps->getVectorLength();

  fftw_plan InputVectorForwardTrafoPlan = fOps->getFFTPlan(input, output2, ExtremeFFT4D_Forward);						       
  randomInput();
							       

  fftw_execute(InputVectorForwardTrafoPlan); 


//  fOps->transformFromXtraSizeArray(input,input);
  ExtremeFFT xFFT;
//  xFFT.SlowRearrange(N, input, output);
//  xFFT.SlowFourierStep(N, 0, output, interim);
//  xFFT.SlowFourierStep(N, 1, interim, output);
//  xFFT.SlowFourierStep(N, 2, output, interim);
//  xFFT.SlowFourierStep(N, 3, interim, output2);  
//  fOps->transformToXtraSizeArray(output2,output2);
//  fOps->transformToXtraSizeArray(input,input);



  xFFT.FastFourierTrafo(L0,L1,L2,L3, input, output, ExtremeFFT4Drrrrrrrrr_Forward);

  fOps->transformFromXtraSizeArray(output,output);
  fOps->transformFromXtraSizeArray(output2,output2);
  int I;
  double diff = 0;
  int correct = 0;
  for (I=0; I<VL; I++) {
    if ((I>=0) && (I<=8)) {
      printf("own: ");output[I].print();
      output2[I].print();
    }
    if (sqrt(sqr(output[I].x-output2[I].x)+sqr(output[I].y-output2[I].y))<1E-10) {
      correct++;
      if (I%8==0) printf("Correct Entry Nr.: %d\n",I);
    } else {
      if (I%8<4) {
        printf("Incorrect index: %d\n",I);
        printf("own: ");output[I].print();
        output2[I].print();
        exit(0);
      }
    }
    diff += sqr(output[I].x-output2[I].x) + sqr(output[I].y-output2[I].y);
  }
  printf("Difference %1.15f:\n",diff);
  printf("Number of correct entries %d:\n",correct);
  bool usexFFT = fOps->getxFFTusage();

  fOps->setxFFTusage(false);
  startTimer();
  int count = 0;
  while (timePassed()<10) {
    int I;
    for (I=0; I<10; I++) {
      fOps->performFFT(input, output, ExtremeFFT4D_Forward);
    }
    count+=10;
  }
  printf("%d iterations within %1.2f (CPU: %1.2f) seconds (FFTW) --> %d cycles.\n",count,timePassed(), cpuTimePassed(), (int)(timePassed()*2.4E9/count));

  fOps->setxFFTusage(true);
  startTimer();
  count = 0;
  while (timePassed()<10) {
    int I;
    for (I=0; I<10; I++) {
      fOps->performFFT(input, output, ExtremeFFT4D_Forward);
    }
    count+=10;
  }
  printf("%d iterations within %1.2f seconds (xFFT) --> %d cycles.\n",count,timePassed(), (int)(timePassed()*2.4E9/count));

  fOps->setxFFTusage(usexFFT);
}*/


void testFourierTrafo4D(bool forw) {
  printf("\nTesting 4D-Fast-Fourier-Transformation ...\n");
  int VL = fOps->getVectorLength();
  randomInput();
  randomOutput();
  randomOutput2();
  bool usexFFT = fOps->getxFFTusage();
  
  
  fOps->setxFFTusage(false);
  fOps->performFFT(input, output, forw);
  startTimer();
  int count = 0;
  while (timePassed()<10) {
    int I;
    for (I=0; I<10; I++) {
      fOps->performFFT(input, output, forw);
    }
    count+=10;
  }
  printf("%d iterations within %1.2f (CPU: %1.2f) seconds (FFTW) --> %d cycles.\n",count,timePassed(), cpuTimePassed(), (int)(timePassed()*2.4E9/count));


  ExtremeFFT4D xFFT(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, 8, 1, NULL, true);
  xFFT.tune(ExtremeFFT4D_TuneLevel_Low, input, output);



  xFFT.FastFourierTrafo(input, output, forw);
  xFFT.stopThreadsAndWaitForTermination();


  fftw_plan InputVectorForwardTrafoPlan = fOps->getFFTPlan(input, output2, forw);
  fftw_execute(InputVectorForwardTrafoPlan); 





  fOps->transformFromXtraSizeArray(output,output);
  fOps->transformFromXtraSizeArray(output2,output2);
  int I;
  double diff = 0;
  int correct = 0;
  for (I=0; I<VL; I++) {
/*    if ((I>=0) && (I<=8)) {
      printf("own: ");output[I].print();
      output2[I].print();
    }*/
    if (sqrt(sqr(output[I].x-output2[I].x)+sqr(output[I].y-output2[I].y))<1E-10) {
      correct++;
//      if (I%8==0) printf("Correct Entry Nr.: %d\n",I);
    } else {
      if (I%8==0) {
//        printf("Incorrect index: %d\n",I);
//        printf("own: ");output[I].print();
//        output2[I].print();
//      exit(0);
      }
    }
    diff += sqr(output[I].x-output2[I].x) + sqr(output[I].y-output2[I].y);
  }
  printf("Difference %1.15f:\n",diff);
  printf("Number of correct entries %d:\n",correct);

  fOps->setxFFTusage(true);
  xFFT.FastFourierTrafo(input, output, forw);
  startTimer();
  count = 0;
  while (timePassed()<10) {
    int I;
    for (I=0; I<10; I++) {
//      fOps->performFFT(input, output, forw);
  xFFT.FastFourierTrafo(input, output, forw);
    }
    count+=10;
  }
  printf("%d iterations within %1.2f seconds (xFFT4D) --> %d cycles.\n",count,timePassed(), (int)(timePassed()*2.4E9/count));

  fOps->setxFFTusage(usexFFT);
}


int measureL1Ways(int L1CacheSize) {
  printf("Measuring L1-Ways...\n");
  int* Dummy1;
  int* Dummy2;
  int* Dummy3;
  int* Dummy4;
  int* Dummy5;
  int* Dummy6;
  int* Dummy7;
  int* Dummy8;
  int* Dummy9;
  int* Dummy10;
  long int Increment = 2*L1CacheSize;
  Complex* xxx = createSuperAlignedComplex(L1CacheSize*100);  //192 MB
  int I;
  for (I=0; I<L1CacheSize*100; I++) {
    xxx[I].x = zufall();
    xxx[I].y = zufall();
  }
  double T[9];
  
  printf("Testing 1 way\n");
  T[0] = zeitwert();
  __asm__ volatile (
          "mov %10, %0\n\t"
	  
          "mov %11, %1\n\t" "mov %1, %2\n\t" "add %1, %2\n\t" "mov %2, %3\n\t" "add %1, %3\n\t"	"mov %3, %4\n\t"
          "add %1, %4\n\t"  "mov %4, %5\n\t" "add %1, %5\n\t" "mov %5, %6\n\t" "add %1, %6\n\t" "mov %6, %7\n\t"
          "add %1, %7\n\t"  "mov %7, %8\n\t" "add %1, %8\n\t"	  
	  
          "movapd (%0), %%xmm0\n\t"
          "mov $100000000, %9\n\t"
          "L1WaysMeasureW1Loop:\n\t"

            "movapd %%xmm0, (%0)\n\t"

            "dec %9\n\t"
          "jg L1WaysMeasureW1Loop\n\t"
	  : "=r" (Dummy1), "=r" (Dummy2), "=r" (Dummy3), "=r" (Dummy4), "=r" (Dummy5), "=r" (Dummy6), "=r" (Dummy7), "=r" (Dummy8), "=r" (Dummy9), "=r" (Dummy10)
	  : "m" (xxx), "m" (Increment));
  T[0] = zeitwert() - T[0];
  
  printf("Testing 2 ways\n");
  T[1] = zeitwert();
  __asm__ volatile (
          "mov %10, %0\n\t"

          "mov %11, %1\n\t" "mov %1, %2\n\t" "add %1, %2\n\t" "mov %2, %3\n\t" "add %1, %3\n\t"	"mov %3, %4\n\t"
          "add %1, %4\n\t"  "mov %4, %5\n\t" "add %1, %5\n\t" "mov %5, %6\n\t" "add %1, %6\n\t" "mov %6, %7\n\t"
          "add %1, %7\n\t"  "mov %7, %8\n\t" "add %1, %8\n\t"	  

          "movapd (%0), %%xmm0\n\t"
          "mov $100000000, %9\n\t"
          "L1WaysMeasureW2Loop:\n\t"

            "movapd %%xmm0, (%0)\n\t"
            "movapd %%xmm0, (%0,%1)\n\t"

            "dec %9\n\t"
          "jg L1WaysMeasureW2Loop\n\t"
	  : "=r" (Dummy1), "=r" (Dummy2), "=r" (Dummy3), "=r" (Dummy4), "=r" (Dummy5), "=r" (Dummy6), "=r" (Dummy7), "=r" (Dummy8), "=r" (Dummy9), "=r" (Dummy10)
	  : "m" (xxx), "m" (Increment));
  T[1] = zeitwert() - T[1];
  
  printf("Testing 3 ways\n");
  T[2] = zeitwert();
  __asm__ volatile (
          "mov %10, %0\n\t"

          "mov %11, %1\n\t" "mov %1, %2\n\t" "add %1, %2\n\t" "mov %2, %3\n\t" "add %1, %3\n\t"	"mov %3, %4\n\t"
          "add %1, %4\n\t"  "mov %4, %5\n\t" "add %1, %5\n\t" "mov %5, %6\n\t" "add %1, %6\n\t" "mov %6, %7\n\t"
          "add %1, %7\n\t"  "mov %7, %8\n\t" "add %1, %8\n\t"	  

          "movapd (%0), %%xmm0\n\t"
          "mov $100000000, %9\n\t"
          "L1WaysMeasureW3Loop:\n\t"

            "movapd %%xmm0, (%0)\n\t"
            "movapd %%xmm0, (%0,%1)\n\t"
            "movapd %%xmm0, (%0,%2)\n\t"

            "dec %9\n\t"
          "jg L1WaysMeasureW3Loop\n\t"
	  : "=r" (Dummy1), "=r" (Dummy2), "=r" (Dummy3), "=r" (Dummy4), "=r" (Dummy5), "=r" (Dummy6), "=r" (Dummy7), "=r" (Dummy8), "=r" (Dummy9), "=r" (Dummy10)
	  : "m" (xxx), "m" (Increment));
  T[2] = zeitwert() - T[2];  
  
  printf("Testing 4 ways\n");
  T[3] = zeitwert();
  __asm__ volatile (
          "mov %10, %0\n\t"

          "mov %11, %1\n\t" "mov %1, %2\n\t" "add %1, %2\n\t" "mov %2, %3\n\t" "add %1, %3\n\t"	"mov %3, %4\n\t"
          "add %1, %4\n\t"  "mov %4, %5\n\t" "add %1, %5\n\t" "mov %5, %6\n\t" "add %1, %6\n\t" "mov %6, %7\n\t"
          "add %1, %7\n\t"  "mov %7, %8\n\t" "add %1, %8\n\t"	  

          "movapd (%0), %%xmm0\n\t"
          "mov $100000000, %9\n\t"
          "L1WaysMeasureW4Loop:\n\t"

            "movapd %%xmm0, (%0)\n\t"
            "movapd %%xmm0, (%0,%1)\n\t"
            "movapd %%xmm0, (%0,%2)\n\t"
            "movapd %%xmm0, (%0,%3)\n\t"

            "dec %9\n\t"
          "jg L1WaysMeasureW4Loop\n\t"
	  : "=r" (Dummy1), "=r" (Dummy2), "=r" (Dummy3), "=r" (Dummy4), "=r" (Dummy5), "=r" (Dummy6), "=r" (Dummy7), "=r" (Dummy8), "=r" (Dummy9), "=r" (Dummy10)
	  : "m" (xxx), "m" (Increment));
  T[3] = zeitwert() - T[3];  

  printf("Testing 5 ways\n");
  T[4] = zeitwert();
  __asm__ volatile (
          "mov %10, %0\n\t"

          "mov %11, %1\n\t" "mov %1, %2\n\t" "add %1, %2\n\t" "mov %2, %3\n\t" "add %1, %3\n\t"	"mov %3, %4\n\t"
          "add %1, %4\n\t"  "mov %4, %5\n\t" "add %1, %5\n\t" "mov %5, %6\n\t" "add %1, %6\n\t" "mov %6, %7\n\t"
          "add %1, %7\n\t"  "mov %7, %8\n\t" "add %1, %8\n\t"	  

          "movapd (%0), %%xmm0\n\t"
          "mov $100000000, %9\n\t"
          "L1WaysMeasureW5Loop:\n\t"

            "movapd %%xmm0, (%0)\n\t"
            "movapd %%xmm0, (%0,%1)\n\t"
            "movapd %%xmm0, (%0,%2)\n\t"
            "movapd %%xmm0, (%0,%3)\n\t"
            "movapd %%xmm0, (%0,%4)\n\t"

            "dec %9\n\t"
          "jg L1WaysMeasureW5Loop\n\t"
	  : "=r" (Dummy1), "=r" (Dummy2), "=r" (Dummy3), "=r" (Dummy4), "=r" (Dummy5), "=r" (Dummy6), "=r" (Dummy7), "=r" (Dummy8), "=r" (Dummy9), "=r" (Dummy10)
	  : "m" (xxx), "m" (Increment));
  T[4] = zeitwert() - T[4];  

  printf("Testing 6 ways\n");
  T[5] = zeitwert();
  __asm__ volatile (
          "mov %10, %0\n\t"

          "mov %11, %1\n\t" "mov %1, %2\n\t" "add %1, %2\n\t" "mov %2, %3\n\t" "add %1, %3\n\t"	"mov %3, %4\n\t"
          "add %1, %4\n\t"  "mov %4, %5\n\t" "add %1, %5\n\t" "mov %5, %6\n\t" "add %1, %6\n\t" "mov %6, %7\n\t"
          "add %1, %7\n\t"  "mov %7, %8\n\t" "add %1, %8\n\t"	  

          "movapd (%0), %%xmm0\n\t"
          "mov $100000000, %9\n\t"
          "L1WaysMeasureW6Loop:\n\t"

            "movapd %%xmm0, (%0)\n\t"
            "movapd %%xmm0, (%0,%1)\n\t"
            "movapd %%xmm0, (%0,%2)\n\t"
            "movapd %%xmm0, (%0,%3)\n\t"
            "movapd %%xmm0, (%0,%4)\n\t"
            "movapd %%xmm0, (%0,%5)\n\t"

            "dec %9\n\t"
          "jg L1WaysMeasureW6Loop\n\t"
	  : "=r" (Dummy1), "=r" (Dummy2), "=r" (Dummy3), "=r" (Dummy4), "=r" (Dummy5), "=r" (Dummy6), "=r" (Dummy7), "=r" (Dummy8), "=r" (Dummy9), "=r" (Dummy10)
	  : "m" (xxx), "m" (Increment));
  T[5] = zeitwert() - T[5];  

  printf("Testing 7 ways\n");
  T[6] = zeitwert();
  __asm__ volatile (
          "mov %10, %0\n\t"

          "mov %11, %1\n\t" "mov %1, %2\n\t" "add %1, %2\n\t" "mov %2, %3\n\t" "add %1, %3\n\t"	"mov %3, %4\n\t"
          "add %1, %4\n\t"  "mov %4, %5\n\t" "add %1, %5\n\t" "mov %5, %6\n\t" "add %1, %6\n\t" "mov %6, %7\n\t"
          "add %1, %7\n\t"  "mov %7, %8\n\t" "add %1, %8\n\t"	  

          "movapd (%0), %%xmm0\n\t"
          "mov $100000000, %9\n\t"
          "L1WaysMeasureW7Loop:\n\t"

            "movapd %%xmm0, (%0)\n\t"
            "movapd %%xmm0, (%0,%1)\n\t"
            "movapd %%xmm0, (%0,%2)\n\t"
            "movapd %%xmm0, (%0,%3)\n\t"
            "movapd %%xmm0, (%0,%4)\n\t"
            "movapd %%xmm0, (%0,%5)\n\t"
            "movapd %%xmm0, (%0,%6)\n\t"

            "dec %9\n\t"
          "jg L1WaysMeasureW7Loop\n\t"
	  : "=r" (Dummy1), "=r" (Dummy2), "=r" (Dummy3), "=r" (Dummy4), "=r" (Dummy5), "=r" (Dummy6), "=r" (Dummy7), "=r" (Dummy8), "=r" (Dummy9), "=r" (Dummy10)
	  : "m" (xxx), "m" (Increment));
  T[6] = zeitwert() - T[6];  
  
  printf("Testing 8 ways\n");
  T[7] = zeitwert();
  __asm__ volatile (
          "mov %10, %0\n\t"

          "mov %11, %1\n\t" "mov %1, %2\n\t" "add %1, %2\n\t" "mov %2, %3\n\t" "add %1, %3\n\t"	"mov %3, %4\n\t"
          "add %1, %4\n\t"  "mov %4, %5\n\t" "add %1, %5\n\t" "mov %5, %6\n\t" "add %1, %6\n\t" "mov %6, %7\n\t"
          "add %1, %7\n\t"  "mov %7, %8\n\t" "add %1, %8\n\t"	  

          "movapd (%0), %%xmm0\n\t"
          "mov $100000000, %9\n\t"
          "L1WaysMeasureW8Loop:\n\t"

            "movapd %%xmm0, (%0)\n\t"
            "movapd %%xmm0, (%0,%1)\n\t"
            "movapd %%xmm0, (%0,%2)\n\t"
            "movapd %%xmm0, (%0,%3)\n\t"
            "movapd %%xmm0, (%0,%4)\n\t"
            "movapd %%xmm0, (%0,%5)\n\t"
            "movapd %%xmm0, (%0,%6)\n\t"
            "movapd %%xmm0, (%0,%7)\n\t"

            "dec %9\n\t"
          "jg L1WaysMeasureW8Loop\n\t"
	  : "=r" (Dummy1), "=r" (Dummy2), "=r" (Dummy3), "=r" (Dummy4), "=r" (Dummy5), "=r" (Dummy6), "=r" (Dummy7), "=r" (Dummy8), "=r" (Dummy9), "=r" (Dummy10)
	  : "m" (xxx), "m" (Increment));
  T[7] = zeitwert() - T[7];  
  
  printf("Testing 9 ways\n");
  T[8] = zeitwert();
  __asm__ volatile (
          "mov %10, %0\n\t"

          "mov %11, %1\n\t" "mov %1, %2\n\t" "add %1, %2\n\t" "mov %2, %3\n\t" "add %1, %3\n\t"	"mov %3, %4\n\t"
          "add %1, %4\n\t"  "mov %4, %5\n\t" "add %1, %5\n\t" "mov %5, %6\n\t" "add %1, %6\n\t" "mov %6, %7\n\t"
          "add %1, %7\n\t"  "mov %7, %8\n\t" "add %1, %8\n\t"	  

          "movapd (%0), %%xmm0\n\t"
          "mov $100000000, %9\n\t"
          "L1WaysMeasureW9Loop:\n\t"

            "movapd %%xmm0, (%0)\n\t"
            "movapd %%xmm0, (%0,%1)\n\t"
            "movapd %%xmm0, (%0,%2)\n\t"
            "movapd %%xmm0, (%0,%3)\n\t"
            "movapd %%xmm0, (%0,%4)\n\t"
            "movapd %%xmm0, (%0,%5)\n\t"
            "movapd %%xmm0, (%0,%6)\n\t"
            "movapd %%xmm0, (%0,%7)\n\t"
            "movapd %%xmm0, (%0,%8)\n\t"

            "dec %9\n\t"
          "jg L1WaysMeasureW9Loop\n\t"
	  : "=r" (Dummy1), "=r" (Dummy2), "=r" (Dummy3), "=r" (Dummy4), "=r" (Dummy5), "=r" (Dummy6), "=r" (Dummy7), "=r" (Dummy8), "=r" (Dummy9), "=r" (Dummy10)
	  : "m" (xxx), "m" (Increment));
  T[8] = zeitwert() - T[8];  

  int ways = 1;
  double cmp = (T[8]/(8+1)) / 5;
  for (I=0; I<9; I++) {
    printf("Anzahl Writes pro Cache-Line: %d ==> %1.3f seconds, average per write: %1.3f\n",I+1,T[I],T[I]/(I+1));
    if ((T[I]/(I+1)) < cmp) ways = I+1;
  }
  printf("Number of detected L1-ways: %d\n", ways);
  return ways;
}


void testFFTembedding(int L1CacheSize, int ways, bool print, int* conflicts) {
  int pos[16];
  int count = 0;
  int I0,I1,I2,I3;
  int modulo = L1CacheSize / ways;
  int N0half = L0 / 2;
  int N1half = L1 / 2;
  int N2half = L2 / 2;
  int N3half = L3 / 2;
  int* hits = new int[modulo];
  
  for (I0=0; I0<2; I0++) {
    for (I1=0; I1<2; I1++) {
      for (I2=0; I2<2; I2++) {
        for (I3=0; I3<2; I3++) {
	  pos[count] = I3*N3half + I2*N2half*(L3+xtraSize3) + I1*N1half*(L3+xtraSize3)*(L2+xtraSize2) + I0*N0half*(L3+xtraSize3)*(L2+xtraSize2)*(L1+xtraSize1);
  	  pos[count] *= 128;
	  
  	  pos[count] = pos[count] % modulo;
	  
	  if (print) printf("FFT-Embedded Position %d is: %d\n",count,pos[count]);
          count++;	
	}
      }
    }
  }

  for (I0=0; I0<modulo; I0++) hits[I0] = 0;
  for (I0=0; I0<16; I0++) {
    for (I1=0; I1<N3half*128; I1++) {
      hits[(pos[I0]+I1) % modulo]++;
      hits[(pos[I0]+I1+(L3+xtraSize3)*128) % modulo]++;
    }
  }
  for (I0=0; I0<32; I0++) {
    conflicts[I0] = 0;
  }
  for (I0=0; I0<modulo; I0++) {
    if (hits[I0]>0) {
      conflicts[hits[I0]-1]++;
    }
  }

  int severeHits = 0;
  for (I0=1; I0<32; I0++) {
    if (print) printf("Number of hit-%d-conflicts: %d\n",I0+1,conflicts[I0]);
    if ((I0+1) > ways) severeHits += conflicts[I0];
  }
  if (print) printf("Number of severe Hits: %d (Number of L1-ways is %d)\n",severeHits, ways);
  delete[] hits;
}


void findOptimalEmbedding(int L1CacheSize, int ways, int max) {
  int e1, e2, e3;
  int* bestConf = new int[32];
  int bestSize = L0*L1*L2*L3;
  int bestE1 = 0;
  int bestE2 = 0;
  int bestE3 = 0;
  int* conflicts = new int[32];
  int I,I2;
  for (I=0; I<32; I++) bestConf[I] = 2*L0*L1*L2*L3*128;
  
  for (e3=0; e3<max; e3++) {
    for (e2=0; e2<max; e2++) {
      for (e1=0; e1<max; e1++) {
        xtraSize1 = e1;
        xtraSize2 = e2;
        xtraSize3 = e3;
	
	testFFTembedding(L1CacheSize, ways, false, conflicts);
	bool better = false;
        bool precEqual = true;
	for (I=0; I<32; I++) {
	  precEqual = true;
	  for (I2=I+1; I2<32; I2++) {
	    if (bestConf[I2] != conflicts[I2]) precEqual = false;
	  }
	  if ((precEqual) && (conflicts[I]<bestConf[I])) {
	    better = true;
	  }
	}
	precEqual = true;
	for (I=0; I<32; I++) {
	  if (bestConf[I] != conflicts[I]) precEqual = false;	  
	}
	if (precEqual) {
          if ((L0*(L3+e3)*(L2+e2)*(L1+e1)) < bestSize) {
	    better = true;
	  }
	}

        if (better) {	
	  for (I=0; I<32; I++) {
  	    bestConf[I] = conflicts[I];
	  }
	  bestE1 = e1;
	  bestE2 = e2;
	  bestE3 = e3;
	  bestSize = L0*(L3+e3)*(L2+e2)*(L1+e1);
	}
      }
    }
  }
  xtraSize1 = bestE1;
  xtraSize2 = bestE2;
  xtraSize3 = bestE3;
  printf("Best embedding: e3=%d, e2=%d, e1=%d\n",xtraSize3,xtraSize2,xtraSize1);
  
  testFFTembedding(L1CacheSize, ways, true, conflicts);
  delete[] conflicts;
}


void testArpackEWfinder(int what, double yN, bool compare, int nev, double TOL, bool doubleM) {
  printf("\nTesting Arpack EW-Finder with yN = %1.2f, nev = %d, TOL = %1.15f, and what=%d...\n",yN,nev, TOL,what);
  if (doubleM) {
    printf("Applying double M\n");
  } else {
    printf("Applying single M\n");  
  }
  randomPhi(0);
  randomInput();
  
  ComplexVector inputVec(fOps->getVectorLength(),input);
  ComplexVector outputVec(fOps->getVectorLength(),output);
  ComplexVector outputVec2(fOps->getVectorLength(),output2);
  
  fOps->setYukawaCoupling(yN);

  int VL = fOps->getVectorLength();
  int I;
  Complex* EigenV = fOps->calcFermionMatrixARPACKEigenValues(what,nev, phiField, TOL, doubleM, NULL,false,true);   
  bool changed=true;
  while (changed) {
    changed = false;
    for (I=0; I<nev-1; I++) {
      if ((what / 2) == 0) { 
        if (EigenV[I].x < EigenV[I+1].x) {
	  Complex dummy = EigenV[I];
  	  EigenV[I] = EigenV[I+1];
	  EigenV[I+1] = dummy;
          changed = true;	
        }
      }
      if ((what / 2) == 1) { 
        if (EigenV[I].y < EigenV[I+1].y) {
	  Complex dummy = EigenV[I];
  	  EigenV[I] = EigenV[I+1];
	  EigenV[I+1] = dummy;
          changed = true;	
        }
      }
      if ((what / 2) == 2) { 
        if (sqr(EigenV[I].x)+sqr(EigenV[I].y) < sqr(EigenV[I+1].x)+sqr(EigenV[I+1].y)) {
	  Complex dummy = EigenV[I];
  	  EigenV[I] = EigenV[I+1];
	  EigenV[I+1] = dummy;
          changed = true;	
        }
      }
    }
  }

  printf("Searched for the first %d eigenvalues.\n",nev);
  for (I=0; I<nev; I++) {
    printf("  -> EV %d: %1.15f %1.15fi\n",I,EigenV[I].x,EigenV[I].y);
  }

  if (compare) {    
    ComplexMatrix mat(1);
    if (doubleM) {
      ComplexMatrix mat2(1);    
      ComplexMatrix mat3(1);    
      fOps->constructNeubergerWithXiFermionMatrix(mat2, (vector4D*) phiField);
      fOps->constructNeubergerWithXiFermionMatrix(mat3, (vector4D*) phiField);
      mat3.dagger();
      mat = mat2*mat3;
    } else {
      fOps->constructNeubergerWithXiFermionMatrix(mat, (vector4D*) phiField);
    }
    mat.calcEigenvalues();
    for (I=0; I<nev; I++) {
      int I2;
      double bestNorm = 1E10;
      int bestInd = I;
      for (I2=I; I2<VL; I2++) {
        double norm = sqr(EigenV[I].x-mat.eigenvalues[I2].x) + sqr(EigenV[I].y-mat.eigenvalues[I2].y);
        if (norm<bestNorm) {
	  bestNorm = norm;
	  bestInd = I2;
	}
      }
      Complex dummy = mat.eigenvalues[I];
      mat.eigenvalues[I] = mat.eigenvalues[bestInd];
      mat.eigenvalues[bestInd] = dummy;
    }
    
    printf("Exact first %d eigenvalues.\n",nev);
    for (I=0; I<nev; I++) {
      printf("  -> EV %d: %1.15f %1.15fi\n",I,mat.eigenvalues[I].x,mat.eigenvalues[I].y);
    }
    double diff = 0;
    for (I=0; I<nev-1; I++) {
      diff += sqr(EigenV[I].x-mat.eigenvalues[I].x) + sqr(EigenV[I].y-mat.eigenvalues[I].y);
    }
    diff = sqrt(diff);
    printf("Difference between results: %1.15f\n",diff);
  }
  printf("Test of Arpacks Eigenvalues READY!\n\n");
}


void testPHMCdSdOmegaNeubergerWithChi(double yN, double massSplit, double explicitMass, double h, bool precUse, double precM, double precS, bool QprecUse, double QprecMu, double QprecBeta, bool RprecUse, double RprecM, double RprecF) {
  printf("\nTesting pHMC dS/dOmega (Neuberger with Xi) with yN = %1.2f, h = %f ...\n",yN,h);
  randomPhi(0);
  
  fOps->setYukawaCoupling(yN);
  fOps->setPreconditioner(precUse, precM, precS);
  fOps->setQPreconditioner(QprecUse, QprecMu, QprecBeta);
  fOps->setRPreconditioner(RprecUse, RprecM, RprecF);
  fOps->setMassSplitRatio(massSplit);
  fOps->setExplicitMass(explicitMass);
  pHMCForce* pHMCforce = new pHMCForce(fOps, true, 1,0, 0); 
  
  int rootCount = 1;
  Complex* roots = new Complex[2*rootCount];
  int I;
  for (I=0; I<rootCount; I++) {
    roots[2*I] = Complex(zufall(),zufall());
    roots[2*I+1] = roots[2*I];
    roots[2*I+1].y = -roots[2*I+1].y;    
  }
  pHMCforce->setApproxPolyRoots(0, roots, 2*rootCount, 1.6, 1.0);
	pHMCforce->randomOmegaField();  
  pHMCforce->setQuasiHermiteanMode(true);
  double S = pHMCforce->calcOmegaAction(0, phiField);
  printf(" -> Omega Action: %f\n",S);

  int t0 = (int) (zufall()*L0);
  int t1 = (int) (zufall()*L1);
  int t2 = (int) (zufall()*L2);
  int t3 = (int) (zufall()*L3);
  int tx = (int) (zufall()*8);
  int testIndex = tx + 8*t3 + 8*(L3+xtraSize3)*t2 + 8*(L3+xtraSize3)*(L2+xtraSize2)*t1 + 8*(L3+xtraSize3)*(L2+xtraSize2)*(L1+xtraSize1)*t0;
  Complex* omega = pHMCforce->getOmega();
  omega[testIndex].x += 0.5*h;
  Complex dSdOmega;
  dSdOmega.x = pHMCforce->calcOmegaAction(0, phiField);
  omega[testIndex].x -= h;
  dSdOmega.x -= pHMCforce->calcOmegaAction(0, phiField);
  omega[testIndex].x += 0.5*h;
  dSdOmega.x /= h;
  omega[testIndex].y += 0.5*h;
  dSdOmega.y = pHMCforce->calcOmegaAction(0, phiField);
  omega[testIndex].y -= h;
  dSdOmega.y -= pHMCforce->calcOmegaAction(0, phiField);
  omega[testIndex].y += 0.5*h;
  dSdOmega.y /= h;  
  printf(" ->Numerical result for derivative: ");
  dSdOmega.print();
    
  startTimer();
  Complex* dSdOmegaExact = NULL;
  int count = 0;
  while (timePassed()<10) {
    dSdOmegaExact = pHMCforce->calcAllForces(0, phiField, S);
    count++;
  }
  printf(" ==> For %d roots: %d iterations within %1.2f seconds.\n",2*rootCount, count,timePassed());

  printf(" ->Exact result for derivative: ");
  dSdOmegaExact[testIndex].print();
  
  printf(" ->Comparison with numerical solution:\n");
  printf("   Relative difference Norm for dSdOmega_%d = %1.15f+%1.15f i \n",testIndex,
   (dSdOmega.x-dSdOmegaExact[testIndex].x) / (dSdOmegaExact[testIndex].x),
   (dSdOmega.y-dSdOmegaExact[testIndex].y) / (dSdOmegaExact[testIndex].y));


  delete pHMCforce;
}


void testDistributedPHMCPolynomialMMdagInverseSQRTomegaAction(int OpMode, bool enforceCoreTie, double yN, double massSplit, double explicitMass, bool precUse, double precM, double precS, bool QprecUse, double QprecMu, double QprecBeta, bool RprecUse, double RprecM, double RprecF) {
  printf("\nTesting Distributed pHMC Polynomial MMdag Inverse SQRT Omega Action with yN = %1.2f, explicitMass = %1.2f, and massSplit = %1.2f...\n",yN, explicitMass, massSplit);
  randomPhi(0);
  randomInput();
  int VLxtr = fOps->getVectorLengthXtrSize();
  
  fOps->setYukawaCoupling(yN);
  fOps->setPreconditioner(precUse, precM, precS);
  fOps->setQPreconditioner(QprecUse, QprecMu, QprecBeta);
  fOps->setRPreconditioner(RprecUse, RprecM, RprecF);
  fOps->setMassSplitRatio(massSplit);
  fOps->setExplicitMass(explicitMass);
  
  int rootCount = 10;
  Complex* roots = new Complex[2*rootCount];
  int I;
  for (I=0; I<rootCount; I++) {
    roots[2*I] = Complex(1*zufall(),1*zufall());
    roots[2*I+1] = roots[2*I];
    roots[2*I+1].y = -roots[2*I+1].y;    
  }

  fOps->deactivateMultiThreadedOps();  
  pHMCForce* pHMCforce1 = new pHMCForce(fOps, true, 1, 0, 30); 
  pHMCforce1->setApproxPolyRoots(0, roots, 2*rootCount, 1.5E-4, 1.0);  
  Complex* omegaField1 = pHMCforce1->getOmega();
  SSE_ZCopy(VLxtr, input, 1, omegaField1, 1);
  pHMCforce1->setQuasiHermiteanMode(true);

  pHMCforce1->calcOmegaAction(0, phiField);
  startTimer();
  int count = 0;
  while (timePassed()<10) {
    pHMCforce1->calcOmegaAction(0, phiField);
    count++;
  }
  printf(" ==> For %d roots: %d iterations within %1.2f seconds.\n",2*rootCount, count,timePassed());
  double polS1 = pHMCforce1->getActOmegaAction(0);

  fOps->activateMultiThreadedOps(OpMode, enforceCoreTie);  
  pHMCForce* pHMCforce2 = new pHMCForce(fOps, true, 1, 0, 40); 
  pHMCforce2->setApproxPolyRoots(0, roots, 2*rootCount, 1.5E-4, 1.0);  
  Complex* omegaField2 = pHMCforce2->getOmega();
  SSE_ZCopy(VLxtr, input, 1, omegaField2, 1);
  pHMCforce2->setQuasiHermiteanMode(true);

  pHMCforce2->calcOmegaAction(0, phiField);
  startTimer();
  count = 0;
  while (timePassed()<10) {
    pHMCforce2->calcOmegaAction(0, phiField);
    count++;
  }
  printf(" ==> For %d roots (Distributed): %d iterations within %1.2f seconds.\n",2*rootCount, count,timePassed());
  double polS2 = pHMCforce2->getActOmegaAction(0);


  double diff = fabs(polS2 - polS1);
  printf("Difference between exact actions: %1.15f (relative %1.15e)\n", diff, diff/polS1);

  delete pHMCforce2;
  fOps->deactivateMultiThreadedOps();  
  delete pHMCforce1;
  delete[] roots;  
}


void testDistributedPHMCExactMMdagInverseSQRTomegaAction(int OpMode, bool enforceCoreTie, double yN, double massSplit, double explicitMass, bool precUse, double precM, double precS, bool QprecUse, double QprecMu, double QprecBeta, bool RprecUse, double RprecM, double RprecF) {
  printf("\nTesting Distributed pHMC Exact MMdag Inverse SQRT Omega Action with yN = %1.2f, explicitMass = %1.2f, and massSplit = %1.2f...\n",yN, explicitMass, massSplit);
  randomPhi(0);
  randomInput();
  int VLxtr = fOps->getVectorLengthXtrSize();
  
  fOps->setYukawaCoupling(yN);
  fOps->setPreconditioner(precUse, precM, precS);
  fOps->setQPreconditioner(QprecUse, QprecMu, QprecBeta);
  fOps->setRPreconditioner(RprecUse, RprecM, RprecF);
  fOps->setMassSplitRatio(massSplit);
  fOps->setExplicitMass(explicitMass);
  
  int rootCount = 10;
  Complex* roots = new Complex[2*rootCount];
  int I;
  for (I=0; I<rootCount; I++) {
    roots[2*I] = Complex(1*zufall(),1*zufall());
    roots[2*I+1] = roots[2*I];
    roots[2*I+1].y = -roots[2*I+1].y;    
  }

  fOps->deactivateMultiThreadedOps();  
  pHMCForce* pHMCforce1 = new pHMCForce(fOps, true, 1, 0, 30); 
  pHMCforce1->setApproxPolyRoots(0, roots, 2*rootCount, 1.5E-4, 1.0);  
  Complex* omegaField1 = pHMCforce1->getOmega();
  SSE_ZCopy(VLxtr, input, 1, omegaField1, 1);
  pHMCforce1->setQuasiHermiteanMode(true);

  int neededIter = 0;
  pHMCforce1->calcExactMMdagInverseSQRTomegaAction(phiField, 0.5, neededIter);
  startTimer();
  int count = 0;
  while (timePassed()<10) {
    pHMCforce1->calcExactMMdagInverseSQRTomegaAction(phiField, 0.5, neededIter);
    count++;
  }
  printf(" ==> For %d roots: %d iterations within %1.2f seconds.\n",2*rootCount, count,timePassed());
  double exactS1 = pHMCforce1->getExactOmegaMMdagInverseSQRTAction();

  fOps->activateMultiThreadedOps(OpMode, enforceCoreTie);  
  pHMCForce* pHMCforce2 = new pHMCForce(fOps, true, 1, 0, 40); 
  pHMCforce2->setApproxPolyRoots(0, roots, 2*rootCount, 1.5E-4, 1.0);  
  Complex* omegaField2 = pHMCforce2->getOmega();
  SSE_ZCopy(VLxtr, input, 1, omegaField2, 1);
  pHMCforce2->setQuasiHermiteanMode(true);

  pHMCforce2->calcExactMMdagInverseSQRTomegaAction(phiField, 0.5, neededIter);
  startTimer();
  count = 0;
  while (timePassed()<10) {
    pHMCforce2->calcExactMMdagInverseSQRTomegaAction(phiField, 0.5, neededIter);
    count++;
  }
  printf(" ==> For %d roots (Distributed): %d iterations within %1.2f seconds.\n",2*rootCount, count,timePassed());
  double exactS2 = pHMCforce2->getExactOmegaMMdagInverseSQRTAction();


  double diff = fabs(exactS2 - exactS1);
  printf("Difference between exact actions: %1.15f (relative %1.15e)\n", diff, diff/exactS1);

  delete pHMCforce2;
  fOps->deactivateMultiThreadedOps();  
  delete pHMCforce1;
  delete[] roots;
}


void testDistributedPHMCOmegaSampling(int OpMode, bool enforceCoreTie, double yN, double massSplit, double explicitMass, bool precUse, double precM, double precS, bool QprecUse, double QprecMu, double QprecBeta, bool RprecUse, double RprecM, double RprecF) {
  printf("\nTesting Distributed pHMC Omega sampling with yN = %1.2f, explicitMass = %1.2f, and massSplit = %1.2f...\n",yN, explicitMass, massSplit);
  randomPhi(0);
  randomInput();
  int VL = fOps->getVectorLength();
  
  fOps->setYukawaCoupling(yN);
  fOps->setPreconditioner(precUse, precM, precS);
  fOps->setQPreconditioner(QprecUse, QprecMu, QprecBeta);
  fOps->setRPreconditioner(RprecUse, RprecM, RprecF);
  fOps->setMassSplitRatio(massSplit);
  fOps->setExplicitMass(explicitMass);
  
  int rootCount = 10;
  Complex* roots = new Complex[2*rootCount];
  int I;
  for (I=0; I<rootCount; I++) {
    roots[2*I] = Complex(1*zufall(),1*zufall());
    roots[2*I+1] = roots[2*I];
    roots[2*I+1].y = -roots[2*I+1].y;    
  }

  fOps->deactivateMultiThreadedOps();  
  pHMCForce* pHMCforce1 = new pHMCForce(fOps, true, 1, 0, 0); 
  pHMCforce1->setApproxPolyRoots(0, roots, 2*rootCount, 1.5E-4, 1.0);  
  Complex* omegaField1 = pHMCforce1->getOmega();
  pHMCforce1->setQuasiHermiteanMode(true);

  pHMCforce1->sampleOmegaFields(phiField);
  startTimer();
  int count = 0;
  while (timePassed()<10) {
    pHMCforce1->sampleOmegaFields(phiField);
    count++;
  }
  printf(" ==> For %d roots: %d iterations within %1.2f seconds.\n",2*rootCount, count,timePassed());

  AdvancedSeed = -1234;
  AdvancedZufall(AdvancedSeed);
  pHMCforce1->sampleOmegaFields(phiField);

  fOps->activateMultiThreadedOps(OpMode, enforceCoreTie);  
  pHMCForce* pHMCforce2 = new pHMCForce(fOps, true, 1, 0, 0); 
  pHMCforce2->setApproxPolyRoots(0, roots, 2*rootCount, 1.5E-4, 1.0);  
  Complex* omegaField2 = pHMCforce2->getOmega();
  pHMCforce2->setQuasiHermiteanMode(true);

  pHMCforce2->sampleOmegaFields(phiField);
  startTimer();
  count = 0;
  while (timePassed()<10) {
    pHMCforce2->sampleOmegaFields(phiField);
    count++;
  }
  printf(" ==> For %d roots (Distributed): %d iterations within %1.2f seconds.\n",2*rootCount, count,timePassed());

  AdvancedSeed = -1234;
  AdvancedZufall(AdvancedSeed);
  pHMCforce2->sampleOmegaFields(phiField);

  double diff = 0;
  fOps->transformFromXtraSizeArray(omegaField1, omegaField1);
  fOps-> transformFromXtraSizeArray(omegaField2, omegaField2); 
  for (int I=0; I<VL; I++) {
    diff += sqr(omegaField1[I].x-omegaField2[I].x) + sqr(omegaField1[I].y-omegaField2[I].y);
  }
  printf("Difference of Omega-Fields: %1.15f\n", diff);

  delete pHMCforce2;
  fOps->deactivateMultiThreadedOps();  
  delete pHMCforce1;
  delete[] roots;
}


void testDistributedPHMCForcesNeubergerWithChi(int OpMode, bool enforceCoreTie, int threadCount, bool usexFFT, double yN, double massSplit, double explicitMass, bool precUse, double precM, double precS, bool QprecUse, double QprecMu, double QprecBeta, bool RprecUse, double RprecM, double RprecF) {
  printf("\nTesting Distributed pHMC Forces (Neuberger with Xi) with yN = %1.2f, explicitMass = %1.2f, and massSplit = %1.2f...\n",yN, explicitMass, massSplit);
  randomPhi(0);
  randomInput();
  int VL = fOps->getVectorLength();
  int VLxtr = fOps->getVectorLengthXtrSize();
  
  fOps->setYukawaCoupling(yN);
  fOps->setxFFT_DistributedFFT_ThreadCount(threadCount);
  fOps->setPreconditioner(precUse, precM, precS);
  fOps->setQPreconditioner(QprecUse, QprecMu, QprecBeta);
  fOps->setRPreconditioner(RprecUse, RprecM, RprecF);
  fOps->setMassSplitRatio(massSplit);
  fOps->setExplicitMass(explicitMass);
  
  int rootCount = 10;
  Complex* roots = new Complex[2*rootCount];
  int I;
  for (I=0; I<rootCount; I++) {
    roots[2*I] = Complex(1*zufall(),1*zufall());
    roots[2*I+1] = roots[2*I];
    roots[2*I+1].y = -roots[2*I+1].y;    
  }

  fOps->deactivateMultiThreadedOps();  
  fOps->setxFFTusage(false);  
  pHMCForce* pHMCforce1 = new pHMCForce(fOps, true, 1, 0, 0); 
  pHMCforce1->setApproxPolyRoots(0, roots, 2*rootCount, 1.5E-4, 1.0);  
  Complex* omegaField1 = pHMCforce1->getOmega();
  SSE_ZCopy(VLxtr, input, 1,omegaField1 , 1);
  pHMCforce1->setQuasiHermiteanMode(true);
  double S1 = pHMCforce1->calcOmegaAction(0, phiField);
  double* dSdPhiExact1 = (double*)pHMCforce1->getdSdPhi();
  double S1b = 0;
  Complex* dSdOmega1 = NULL;

  dSdOmega1 = pHMCforce1->calcAllForces(0, phiField, S1b);
  startTimer();
  int count = 0;
  while (timePassed()<10) {
    dSdOmega1 = pHMCforce1->calcAllForces(0, phiField, S1b);
    count++;
  }
  printf(" ==> For %d roots: %d iterations within %1.2f seconds.\n",2*rootCount, count,timePassed());


  fOps->activateMultiThreadedOps(OpMode, enforceCoreTie);  
  fOps->setxFFTusage(usexFFT);
  pHMCForce* pHMCforce2 = new pHMCForce(fOps, true, 1, 0, 0); 
  pHMCforce2->setApproxPolyRoots(0, roots, 2*rootCount, 1.5E-4, 1.0);  
  Complex* omegaField2 = pHMCforce2->getOmega();
  SSE_ZCopy(VLxtr, input, 1,omegaField2 , 1);
  pHMCforce2->setQuasiHermiteanMode(true);
  double S2 = pHMCforce2->calcOmegaAction(0, phiField);
  double* dSdPhiExact2 = (double*)pHMCforce2->getdSdPhi();
  double S2b = 0;
  Complex* dSdOmega2 = NULL;

  dSdOmega2 = pHMCforce2->calcAllForces(0, phiField, S2b);
  startTimer();
  count = 0;
  while (timePassed()<10) {
    dSdOmega2 = pHMCforce2->calcAllForces(0, phiField, S2b);
    count++;
  }
  printf(" ==> For %d roots (Distributed): %d iterations within %1.2f seconds.\n",2*rootCount, count,timePassed());

  double diff = 0;
  diff = fabs(S1-S2);
  printf("Difference of actions (a): %1.15f\n", diff);
  diff = 0;
  diff = fabs(S1b-S2b);
  printf("Difference of actions (b): %1.15f\n", diff);
  diff = 0;
  fOps->transformFromXtraSizeArray(dSdOmega1, dSdOmega1);
  fOps->transformFromXtraSizeArray(dSdOmega2, dSdOmega2);
  for (int I=0; I<VL; I++) {
    diff += sqr(dSdOmega1[I].x-dSdOmega2[I].x) + sqr(dSdOmega1[I].y-dSdOmega2[I].y);
  }
  printf("Difference of Omega-Forces: %1.15f\n", diff);
  diff = 0;
  for (int I=0; I<VL/2; I++) {
    diff += sqr(dSdPhiExact1[I]-dSdPhiExact2[I]);
  }
  printf("Difference of Phi-Forces: %1.15f\n", diff);

  delete pHMCforce2;
  fOps->deactivateMultiThreadedOps();  
  delete pHMCforce1;
  delete[] roots;
}


void testPHMCdSdPhiNeubergerWithChi(double yN, double massSplit, double explicitMass, double h, bool precUse, double precM, double precS, bool QprecUse, double QprecMu, double QprecBeta, bool RprecUse, double RprecM, double RprecF) {
  printf("\nTesting pHMC dS/dPhi (Neuberger with Xi) with yN = %1.2f, h = %f ...\n",yN,h);
  randomPhi(0);
  
  fOps->setYukawaCoupling(yN);
  fOps->setPreconditioner(precUse, precM, precS);
  fOps->setQPreconditioner(QprecUse, QprecMu, QprecBeta);
  fOps->setRPreconditioner(RprecUse, RprecM, RprecF);
  fOps->setMassSplitRatio(massSplit);
  fOps->setExplicitMass(explicitMass);
  pHMCForce* pHMCforce = new pHMCForce(fOps, true, 1, 0, 0); 
  
  int rootCount = 10;
  Complex* roots = new Complex[2*rootCount];
  int I;
  for (I=0; I<rootCount; I++) {
    roots[2*I] = Complex(1*zufall(),1*zufall());
    roots[2*I+1] = roots[2*I];
    roots[2*I+1].y = -roots[2*I+1].y;    
  }
  pHMCforce->setApproxPolyRoots(0, roots, 2*rootCount, 1.5E-4, 1.0);
  pHMCforce->randomOmegaField();
  pHMCforce->setQuasiHermiteanMode(true);
  double S = pHMCforce->calcOmegaAction(0, phiField);
  printf(" -> Omega Action: %f\n",S);

  int t0 = (int) (zufall()*L0);
  int t1 = (int) (zufall()*L1);
  int t2 = (int) (zufall()*L2);
  int t3 = (int) (zufall()*L3);
  int tx = (int) (zufall()*4);
  int testIndex = tx + 4*t3 + 4*L3*t2 + 4*L3*L2*t1 + 4*L3*L2*L1*t0;
  phiField[testIndex] += 0.5*h;
  double dSdPhi;
  dSdPhi = pHMCforce->calcOmegaAction(0, phiField);    
  phiField[testIndex] -= h;
  dSdPhi -= pHMCforce->calcOmegaAction(0, phiField);
  phiField[testIndex] += 0.5*h;
  dSdPhi /= h;

  printf(" ->Numerical result for derivative: %f\n",dSdPhi);
    
  startTimer();
  int count = 0;
  while (timePassed()<10) {
    pHMCforce->calcAllForces(0, phiField, S);
    count++;
  }
  printf(" ==> For %d roots: %d iterations within %1.2f seconds.\n",2*rootCount, count,timePassed());


  double* dSdPhiExact = (double*)pHMCforce->getdSdPhi();
  printf(" ->Exact result for derivative: %f\n",0.5*dSdPhiExact[testIndex]);
  
  printf(" ->Comparison with numerical solution:\n");
  printf("   Relative difference Norm for dSdPhi_%d = %1.15f \n",testIndex, (dSdPhi-0.5*dSdPhiExact[testIndex]) / (0.5*dSdPhiExact[testIndex]));
  

  delete pHMCforce;
}


void testPHMCSampleOmegaNeubergerWithChi(double yN, bool precUse, double precM, double precS, bool QprecUse, double QprecMu, double QprecBeta, bool RprecUse, double RprecM, double RprecF, bool useKrylov) {
  printf("\nTesting pHMC OmegaField-Sampling (Neuberger with Xi) with yN = %1.2f ...\n",yN);
  randomPhi(0);
  
  fOps->setYukawaCoupling(yN);
  fOps->setPreconditioner(precUse, precM, precS);
  fOps->setQPreconditioner(QprecUse, QprecMu, QprecBeta);
  fOps->setRPreconditioner(RprecUse, RprecM, RprecF);
  pHMCForce* pHMCforce = new pHMCForce(fOps, true, 1,0, 0); 
  
  int rootCount = 10;
  Complex* roots = new Complex[2*rootCount];
  int I;
  for (I=0; I<rootCount; I++) {
    roots[2*I] = Complex(zufall(),zufall());
    roots[2*I+1] = roots[2*I];
    roots[2*I+1].y = -roots[2*I+1].y;    
  }
  pHMCforce->setApproxPolyRoots(0, roots, 2*rootCount, 1.6E-2, 1.00);
  pHMCforce->randomOmegaField();  
  pHMCforce->setQuasiHermiteanMode(true);

  Complex* omega = pHMCforce->getOmega();
  Complex* dSdOmegaExact = NULL;
  double S = 0;
  dSdOmegaExact = pHMCforce->calcAllForces(0, phiField, S);

  startTimer();
  int count = 0;
  while (timePassed()<10) {
    if (useKrylov) {
      pHMCforce->applyInverseSQRTPolynomialKRYLOV(dSdOmegaExact, output, phiField, 0);
    } else {
      pHMCforce->applyInverseSQRTPolynomialCHEBYSHEV(dSdOmegaExact, output, phiField, 0);    
    }
    count++;
  }
  printf(" ==> For inversion of Polynomial: %d iterations within %1.2f seconds.\n", count,timePassed());
  if (useKrylov) {
    pHMCforce->applyInverseSQRTPolynomialKRYLOV(output, output, phiField, 0);
  } else {
    pHMCforce->applyInverseSQRTPolynomialCHEBYSHEV(output, output, phiField, 0);  
  }
  
  
  ComplexVector inputVec(fOps->getVectorLength(),omega);
  ComplexVector outputVec(fOps->getVectorLength(),output);
  
  fOps->transformFromXtraSizeArray(omega, omega);
  fOps->transformFromXtraSizeArray(output, output);
  ComplexVector diffVec = inputVec - outputVec;
  printf("Comparison with Polynomial-Application: Difference Norm: %1.15f\n",diffVec.getNorm());

  delete pHMCforce;
}


void testpHMCPropagation(int OpMode, bool enforceCoreTie, bool FACC, double yN, double kappa, double lambda, bool sphMode, double sphZeta, double massSplit, double explicitMass, double current, double c6, double c8, double c10, double lam6, double lam8, double lam10, bool precUse, double precM, double precS, bool QprecUse, double QprecMu, double QprecBeta, bool RprecUse, double RprecM, double RprecF) {
  printf("\nTesting Integrators for pHMC with yN = %1.2f, kappa = %f, lambda = %1.15f, SphericalMode = %d, SphericalZeta = %1.3e, current = %1.15f, c6 = %1.15f, c8 = %1.15f, c10 = %1.15f, lam6 = %1.15f, lam8 = %1.15f, lam10 = %1.15f,...\n",yN, kappa, lambda, sphMode, sphZeta,current, c6, c8, c10, lam6, lam8, lam10);
  double S1, S2;
  randomPhi(0);
  
  int nf = 1;
  double gamma = 0.15;
  double theta = 0.05;
  int subPolCnt = 1;
  double* polEps = new double[1+subPolCnt];
  polEps[0] = 0.2; 
  polEps[1] = 0.4;   
  double* polLam = new double[1+subPolCnt]; 
  polLam[0] = 50.0;
  polLam[1] = 50.0;
  int* polDeg = new int[1+subPolCnt];
  polDeg[0] = 20;
  polDeg[1] = 12;  
  int digit = 1000;
  double alpha = 0.5; 
  int maxPolDegPerNod = 24;
  int precMCnt = 0;
  double* precMss = NULL;

  int LevelCount = 2 + subPolCnt;
  int* iterations = new int[LevelCount];
  int* integrators = new int[LevelCount];
  int* subPolNr = new int[LevelCount];
  int I;
  for (I=2; I<LevelCount; I++) subPolNr[I] = 2*subPolCnt-2*I+3;
  subPolNr[0] = 0;
  subPolNr[1] = 2*subPolCnt;  
  double* propTOL = NULL;
  double finalTOL = 0;

  
  fOps->setYukawaCoupling(yN);
  fOps->setPreconditioner(precUse, precM, precS);
  fOps->setQPreconditioner(QprecUse, QprecMu, QprecBeta);
  fOps->setRPreconditioner(RprecUse, RprecM, RprecF);
  fOps->setMassSplitRatio(massSplit);
  fOps->setExplicitMass(explicitMass);
  
  if (OpMode>0) fOps->activateMultiThreadedOps(OpMode, enforceCoreTie);  

  pHMCPropagator* pHMCProp = new pHMCPropagator(fOps, lambda, kappa, current, c6, c8, c10, lam6, lam8, lam10, nf, gamma, sphMode, sphZeta, theta, subPolCnt, polEps, polLam, polDeg, precMCnt, precMss, digit, alpha, maxPolDegPerNod, 0); 
  if (FACC) {
    pHMCProp->setPhiForceFourierType(2, 1.00);
    pHMCProp->calcPhiMomentumMasses(1.0);    
  }
  pHMCProp->getNodesReady();
  pHMCProp->sampleALLMomenta();
  pHMCProp->synchronizedChangeOfQuasiHermiteanMode(true);
 
  int iter = 10;
  while (iter<=280) {
    S1 = pHMCProp->calcTotalAction(1E-12);
    pHMCProp->LeapFrogMarkovStep(iter,  1.0 / iter, 0, 0);
    S2 = pHMCProp->calcTotalAction(1E-12);
    printf("Leap-Frog integration: With %d Inversions achieved dS = %1.15f (=%1.15f if accepted)\n", iter, (pHMCProp->SafterProp-pHMCProp->SbeforeProp), S2-S1);
    pHMCProp->restoreALLfields(true);
    iter *= 2;
  }
  printf("\n");
  iter = 10;
  while (iter<=280) {
    iterations[2] = iter;
    iterations[1] = 10;
    iterations[0] = 10;
    integrators[2] = 0;
    integrators[1] = 0;    
    integrators[0] = 0;
    
    S1 = pHMCProp->calcTotalAction(1E-12);
    pHMCProp->multiTimeScaleMarkovStep(LevelCount, iterations, 1.0/iter, integrators, subPolNr, propTOL, finalTOL);
    S2 = pHMCProp->calcTotalAction(1E-12);
    printf("Nested Leap-Frog/Leap-Frog/Leap-Frog integration (%d/%d nested Steps): With %d main-Inversions achieved dS = %1.15f (=%1.15f if accepted)\n", iterations[1], iterations[0], iter, (pHMCProp->SafterProp-pHMCProp->SbeforeProp), S2-S1);
    pHMCProp->restoreALLfields(true);
    iter *= 2;
  }
  printf("\n");
  iter = 10;
  while (iter<=280) {
    iterations[2] = iter;
    iterations[1] = 10;
    iterations[0] = 2;
    integrators[2] = 0;
    integrators[1] = 0;    
    integrators[0] = 2;

    S1 = pHMCProp->calcTotalAction(1E-12);
    pHMCProp->multiTimeScaleMarkovStep(LevelCount, iterations, 1.0/iter, integrators, subPolNr, propTOL, finalTOL);
    S2 = pHMCProp->calcTotalAction(1E-12);
    printf("Nested Leap-Frog/Leap-Frog/OmelyanO4 integration (%d/%d nested Steps): With %d main-Inversions achieved dS = %1.15f (=%1.15f if accepted)\n", iterations[1], iterations[0], iter, (pHMCProp->SafterProp-pHMCProp->SbeforeProp), S2-S1);
    pHMCProp->restoreALLfields(true);
    iter *= 2;
  }
  printf("\n");
  iter = 5;
  while (2*iter<=280) {
    S1 = pHMCProp->calcTotalAction(1E-12);
    pHMCProp->OmelyanO2MarkovStep(iter,  1.0 / iter, 0, 0);
    S2 = pHMCProp->calcTotalAction(1E-12);
    printf("Omelyan Order 2 integration: With %d Inversions achieved dS = %1.15f (=%1.15f if accepted)\n", 2*iter, (pHMCProp->SafterProp-pHMCProp->SbeforeProp), S2-S1);
    pHMCProp->restoreALLfields(true);
    iter *= 2;
  }
  printf("\n");
  iter = 5;
  while (2*iter<=280) {
    iterations[2] = iter;
    iterations[1] = 10;
    iterations[0] = 10;
    integrators[2] = 1;
    integrators[1] = 0;    
    integrators[0] = 0;

    S1 = pHMCProp->calcTotalAction(1E-12);
    pHMCProp->multiTimeScaleMarkovStep(LevelCount, iterations, 1.0/iter, integrators, subPolNr, propTOL, finalTOL);
    S2 = pHMCProp->calcTotalAction(1E-12);
    printf("Nested OmelyanO2/Leap-Frog/Leap-Frog integration (%d/%d nested Steps): With %d main-Inversions achieved dS = %1.15f (=%1.15f if accepted)\n", iterations[1], iterations[0], 2*iter, (pHMCProp->SafterProp-pHMCProp->SbeforeProp), S2-S1);
    pHMCProp->restoreALLfields(true);
    iter *= 2;
  }
  printf("\n");
  iter = 5;
  while (2*iter<=280) {
    iterations[2] = iter;
    iterations[1] = 5;
    iterations[0] = 2;
    integrators[2] = 1;
    integrators[1] = 1;    
    integrators[0] = 2;

    S1 = pHMCProp->calcTotalAction(1E-12);
    pHMCProp->multiTimeScaleMarkovStep(LevelCount, iterations, 1.0/iter, integrators, subPolNr, propTOL, finalTOL);
    S2 = pHMCProp->calcTotalAction(1E-12);
    printf("Nested OmelyanO2/OmelyanO2/OmelyanO4 integration (%d/%d nested Steps): With %d main-Inversions achieved dS = %1.15f (=%1.15f if accepted)\n", iterations[1], iterations[0], 2*iter, (pHMCProp->SafterProp-pHMCProp->SbeforeProp), S2-S1);
    pHMCProp->restoreALLfields(true);
    iter *= 2;
  }
  printf("\n");
  iter = 2;
  while (5*iter<=280) {
    S1 = pHMCProp->calcTotalAction(1E-12);
    pHMCProp->OmelyanO4MarkovStep(iter,  1.0 / iter, 0, 0);
    S2 = pHMCProp->calcTotalAction(1E-12);
    printf("Omelyan Order 4 integration: With %d Inversions achieved dS = %1.15f (=%1.15f if accepted)\n", 5*iter, (pHMCProp->SafterProp-pHMCProp->SbeforeProp), S2-S1);
    pHMCProp->restoreALLfields(true);
    iter *= 2;
  }
  printf("\n");
  iter = 2;
  while (5*iter<=280) {
    iterations[2] = iter;
    iterations[1] = 10;
    iterations[0] = 10;
    integrators[2] = 2;
    integrators[1] = 0;    
    integrators[0] = 0;
    
    S1 = pHMCProp->calcTotalAction(1E-12);
    pHMCProp->multiTimeScaleMarkovStep(LevelCount, iterations, 1.0/iter, integrators, subPolNr, propTOL, finalTOL);
    S2 = pHMCProp->calcTotalAction(1E-12);
    printf("Nested OmelyanO4/Leap-Frog/Leap-Frog integration (%d/%d nested Steps): With %d main-Inversions achieved dS = %1.15f (=%1.15f if accepted)\n", iterations[1], iterations[0], 5*iter, (pHMCProp->SafterProp-pHMCProp->SbeforeProp), S2-S1);
    pHMCProp->restoreALLfields(true);
    iter *= 2;
  }
  printf("\n");
  iter = 2;
  while (5*iter<=280) {
    iterations[2] = iter;
    iterations[1] = 2;
    iterations[0] = 2;
    integrators[2] = 2;
    integrators[1] = 2;    
    integrators[0] = 2;
    
    S1 = pHMCProp->calcTotalAction(1E-12);
    pHMCProp->multiTimeScaleMarkovStep(LevelCount, iterations, 1.0/iter, integrators, subPolNr, propTOL, finalTOL);
    S2 = pHMCProp->calcTotalAction(1E-12);
    printf("Nested OmelyanO4/OmelyanO4/OmelyanO4 integration (%d/%d nested Steps): With %d main-Inversions achieved dS = %1.15f (=%1.15f if accepted)\n", iterations[1], iterations[0], 5*iter, (pHMCProp->SafterProp-pHMCProp->SbeforeProp), S2-S1);
    pHMCProp->restoreALLfields(true);
    iter *= 2;
  }
  
  delete pHMCProp;  
  if (OpMode>0) fOps->deactivateMultiThreadedOps();    
}


void testDerivativeB(double yN, double massSplit, double explicitMass, double h) {
  printf("\nTesting derivative of B with yN = %1.2f, h = %f ...\n",yN,h);
  randomPhi(0);
  randomInput();
  
  fOps->setYukawaCoupling(yN);
  fOps->setMassSplitRatio(massSplit);
  fOps->setExplicitMass(explicitMass);
  
/*  phiField[0] = 1.0;
  phiField[1] = 2.0;
  phiField[2] = 3.0;
  phiField[3] = 4.0;
  ComplexMatrix* mat = createPhiMatB(phiField, false, massSplit);
  mat->print();
  printf("\n");
  mat = createPhiMatB(phiField, true, massSplit);
  mat->print();
  
  delete mat;*/
  
  
  for (int testIndex=0; testIndex<4; testIndex++) {
    phiField[testIndex] += 0.5*h;
    mulWithPhiMatB(phiField, input, output, massSplit, yN);
    Complex x2(0,0);
    for (int I=0; I<8; I++) x2 = x2 + adj(input[I]) * output[I];
    phiField[testIndex] -= h;
    mulWithPhiMatB(phiField, input, output, massSplit, yN);
    phiField[testIndex] += 0.5*h;
    Complex x1(0,0);
    for (int I=0; I<8; I++) x1 = x1 + adj(input[I]) * output[I];
    Complex dx = (1.0/h) * (x2-x1);

    fOps->executeMultiplicationVectorWithDerivativesOfB(input, input, output);
    printf("Numerical derivative: %1.15e + %1.15e i \n",dx.x,dx.y);
    printf("Exact derivative: %1.15e + %1.15e i \n",yN*output[testIndex].x,yN*output[testIndex].y);
  }  
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
  
xtraSize1 = 1;
xtraSize2 = 1;
xtraSize3 = 1;
  
  fOps = new FermionMatrixOperations(L0, L1, L2, L3, 1.0, 0.5, 1.0);
  fOps->setxFFTusage(usexFFT);  
  input = fOps->createFermionVector();
  output = fOps->createFermionVector();
  output2 = fOps->createFermionVector(); 
  interim = fOps->createFermionVector(); 
  interim2 = fOps->createFermionVector(); 
  
  

//  testFourierTrafo4D(ExtremeFFT4D_Backward);
//  exit(0);

  
  
  HMCProp = new HMCPropagator(fOps, 0.01, 0.10, 10, 2.1);
  phiField = (double*) HMCProp->phiField;
  dSdPhi = (double*) fOps->createFermionVector(2);
  
 /* 
  PolynomialApproximation poly(128, 1000, 0.5, 1E-4);
  poly.plotApproxPolynomials("data/polyApproxCheck.dat");
  PowerPolynomialsOneOverX powerPoly(12, 1E-2, 0, 0.5);
  powerPoly.printToFile();
  powerPoly.calcControlPlot();
	*/
  
  initializePerformanceProfiler("testOpsPerformanceProfile.dat");


/*  int ways = 2; //measureL1Ways(65536);
  int* conflicts = new int[32];
  testFFTembedding(65536, ways, true, conflicts);
  delete[] conflicts;
  findOptimalEmbedding(65536, ways,16);
exit(0);*/

//  testFourierTrafo();

  //testConditionNumberEstimation(10);
//exit(0);
 

  testDerivativeB(1.2, 0.15, 0.0, 1E-5);
  testPHMCdSdOmegaNeubergerWithChi(1.2, 0.15, 0.0, 1E-5, false, 1, 0, false, 0.5, 0.5, false, 1.0, 1.0);
  testPHMCdSdPhiNeubergerWithChi(1.2, 0.15, 0.0, 1E-5, true, 1, 0, true, 0.5, 0.5, true, 1.0, 1.0);
  testPHMCSampleOmegaNeubergerWithChi(1.2, true, 1, 0, true, 0.5, 0.5, true, 1.0, 0.30, false);
  testpHMCPropagation(1, false, true, 1.2, 0.1, NaN, true, 0*0.23, 0.15, 0.0, 1.4, 0.0051, 0.0052, 0.0053, 0.054, 0.055, 0.056, true, 1, 0, true, 1, 0.5, true, 1.0, 1.0); 


//  testDistributedVectorCopy(1, false);
//  testDistributedVectorNorm(0, false);
//  testDistributedScalarProduct(2, false);
//  testDistributedVectorAddition(1, false);
//  testDistributedFourierTransformation(1,false, 2);
//  testDistributedYukawaCouplingMatB(1, false, false, 1.7, 0.15);
//  testDistributedMulWithDerivativesOfMatB(1, false, 1.0);
//  testDistributedDiracOperatorApplication(1, true);
//  testDistributedMultiCombinedOperator_VA1_RorCopy_VAc_RorCopy_D_MatrixMultiplication(1, false, 1, true, false, true);
//  testDistributedRPreconditionerOperatorApplication(1, false);
//  testDistributedPiModeRemoverOperatorApplication(2);
//  testDistributedMMDaggerxQuasiHermiteanNeubergerWithChi(1, false, 2, 1.5, true, true);
//  testDistributedPHMCForcesNeubergerWithChi(1, false, 1, false, 1.2, 0.15, 0.0, true, 1, 0, true, 0.5, 0.5, true, 1.0, 1.0);
//  testDistributedPHMCOmegaSampling(2, 1.2, 0.15, 0.0,  true, 1, 0, true, 0.5, 0.5, true, 1.0, 0.30);
//    testDistributedPHMCExactMMdagInverseSQRTomegaAction(2, 1.2, 0.15, 0.0,  true, 1, 0, true, 0.5, 0.5, true, 1.0, 0.30);
//    testDistributedPHMCPolynomialMMdagInverseSQRTomegaAction(2, 1.2, 0.15, 0.0,  true, 1, 0, true, 0.5, 0.5, true, 1.0, 0.30);
    
    
//  testSSECopy();
//  testSSEScalarProduct();
//  testSSEVectorAddition();
/*  testSSEDiracOperatorApplication();
  testSSEBOperatorApplication(1.5, 1.2);*/
  
//  testArpackEWfinder(1, 0, false, 1, 1E-1, false); 

//  testGaussGenerator();

//  testMxNeubergerWithChi(1.5, 0.15, 0.0, true, false);
//  testMxNeubergerWithChi(1.5, 0.15, 0.0, true, true);
//  testMMDaggerxNeubergerWithChi(1.5, 0.15, 0.0, true);
//  testMMDaggerxQuasiHermiteanNeubergerWithChi(1.5, true, true);
//  testCompactMMDaggerxNeubergerWithChi(1.5, true, true, 1, 0, true, 0.5, 0.5);

  testPreconditioner(0.7, 0.24, 0.03);
  testSolverNeubergerWithChi(5, true, true, 1E-8);
  testMMDaggerSQTRSolverNeubergerWithChi(1.1, 1E-10, true, 1.0, 0, true, 0.5, 0.5, true, 1.0, 1.0, 20);
  testMMDaggerChebyshevPolynomialApplication(1.1, 1E-10, true, 1.0, 0, true, 0.5, 0.5);
  testdSdPhiNeubergerWithChi(1.5, 1E-5, 1E-12);
//  testHMCPropagation(1.5, 0.10, 0.01);


//  testEigenvaluesNeubergerWithChi(0.0);
 
  fOps->destroyFermionVector(input);
  fOps->destroyFermionVector(output);
  fOps->destroyFermionVector(output2);
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
  if (LogLevel>1) printf("TestOps terminated correctly.\n");
}
