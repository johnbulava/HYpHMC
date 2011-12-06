#ifndef MultiThreadedOperations_included
#define MultiThreadedOperations_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <numa.h>
#include <unistd.h>


#include "Global.h"
#include "Complex.h"
#include "DistributedMemoryObject.h"
#include "ThreadControl.h"
#include "Tools.h"
#include "xSSE.h"


#define MultiThreadedOperationsThreadMainLoopWaitCycles 100
#define MultiThreadedOperationsThreadMainLoopSecurityWaitTimeInSecs 4
#define MultiThreadedOperationsGeneralWorkSpaceSize 100



class MultiThreadedOperations {
private:  
  int ThreadCount;
  int OperationMode;
  int threadCompleteMask;
  ThreadControl** Threads;
  xSSE* xSSEObj;
  inline void setParameterOnAllThreads(int nr, void* para);
  inline void setExecutionCommandOnAllThreads(int cmd);
  inline Complex getComplexNumberSumFromAllThreads(int nr); 
  void threadedExecute(int threadMask);
  Complex* generalWorkSpace;

  pthread_mutex_t pThreadMutex;
  
public:    
  MultiThreadedOperations(int opMode, bool enforceCoreTie); 
  ~MultiThreadedOperations();  
  
  DistributedMemoryObject* allocateDistibutedMemory(int sizeInComplexNumbers);
  DistributedMemoryObject* allocateDistibutedFermionVector(int L0, int L1, int L2, int L3);
  void terminateThreads();
  void copyComplexVectorToUniformlyDistributedMemory(Complex* vec, DistributedMemoryObject* memObj, int size);
  void calcSumOfUniformlyDistributedComplexVectors(DistributedMemoryObject* memObj, Complex* vec, int size);
  void copyFermionVectorToDistributedFermionVector(Complex* vec, DistributedMemoryObject* memObj, int L0, int L1, int L2, int L3);
  void copyDistributedFermionVectorToFermionVector(DistributedMemoryObject* memObj, Complex* vec, int L0, int L1, int L2, int L3);
  
  void zeroFermionVector(DistributedMemoryObject* memObj, int L0, int L1, int L2, int L3);
  void zeroUniformlyDistributedComplexVector(DistributedMemoryObject* memObj, int size);

  void vectorAdditionOfFermionVectors(DistributedMemoryObject* x, DistributedMemoryObject* y, Complex& alpha, int L0, int L1, int L2, int L3);
  void vectorAdditionOfUniformlyDistributedComplexVectors(DistributedMemoryObject* x, DistributedMemoryObject* y, Complex& alpha, int size);
  void scalarProductOfFermionVectors(DistributedMemoryObject* memObj1, DistributedMemoryObject* memObj2, Complex& res, int L0, int L1, int L2, int L3);
  void vectorNormOfFermionVector(DistributedMemoryObject* memObj, double& res, int L0, int L1, int L2, int L3);
  void multiplyFermionVectorWithComplexNumber(DistributedMemoryObject* input, DistributedMemoryObject* output, Complex alpha, int L0, int L1, int L2, int L3);
  void multiplyFermionVectorWithTimeIndexedComplexScalars(DistributedMemoryObject* input, DistributedMemoryObject* output, DistributedMemoryObject* scalars, int L0, int L1, int L2, int L3);


  void copyFermionVector(DistributedMemoryObject* input, DistributedMemoryObject* output, int L0, int L1, int L2, int L3);

  void perform_FFTWFourierTransformationOfFermionVector(DistributedMemoryObject* input, DistributedMemoryObject* output, bool forw, int L0, int L1, int L2, int L3);
  void perform_xFFTFourierTransformationOfFermionVector(DistributedMemoryObject* input, DistributedMemoryObject* output, bool forw, int threadCount, int L0, int L1, int L2, int L3);
  void tune_xFFTFourierTransformationOfFermionVector(DistributedMemoryObject* input, DistributedMemoryObject* output, int threadCount, int L0, int L1, int L2, int L3, int tuneLevel);
  void tune_xFFTFourierTransformationOfFermionVector(DistributedMemoryObject* input, DistributedMemoryObject* output, int threadCount, int L0, int L1, int L2, int L3, int tuneLevel, char* &fftPlanDescriptor);
  bool setFFTPlan_xFFTFourierTransformationOfFermionVector(int threadCount, int L0, int L1, int L2, int L3, char* fftPlanDescriptor);
  
  void perform_PiModeRemoverOperatorMultiplicationInFourierSpace(DistributedMemoryObject* input, DistributedMemoryObject* output, int L0, int L1, int L2, int L3, DistributedMemoryObject* Index_PiModes);  
  void perform_DiracTypeOperatorMultiplicationInFourierSpace(DistributedMemoryObject* input, DistributedMemoryObject* output, int L0, int L1, int L2, int L3, DistributedMemoryObject* sinP, DistributedMemoryObject* operatorData);
  void perform_yBOperatorMultiplication(DistributedMemoryObject* input, DistributedMemoryObject* output, int L0, int L1, int L2, int L3, double y, double split, DistributedMemoryObject* phiObj);
  void perform_yBDaggerOperatorMultiplication(DistributedMemoryObject* input, DistributedMemoryObject* output, int L0, int L1, int L2, int L3, double y, double split, DistributedMemoryObject* phiObj);
  void perform_RPreconditionerTypeMatrixMultiplicationInFourierSpace(DistributedMemoryObject* input, DistributedMemoryObject* output, int L0, int L1, int L2, int L3, DistributedMemoryObject* rPrecData);
  void perform_MultiplicationVectorWithDerivativesOfB(DistributedMemoryObject* leftInput, DistributedMemoryObject* rightInput, int L0, int L1, int L2, int L3, double split, DistributedMemoryObject* output);
  void perform_MultiCombinedOperator_VA1_RorCopy_VAc_RorCopy_D_InFourierSpace(DistributedMemoryObject* v1, DistributedMemoryObject* v2, DistributedMemoryObject* v3, DistributedMemoryObject* v4, int L0, int L1, int L2, int L3, bool useR, Complex alpha2, DistributedMemoryObject* sinP, DistributedMemoryObject* operatorData, DistributedMemoryObject* rPrecData1, DistributedMemoryObject* rPrecData2);
};


#endif
