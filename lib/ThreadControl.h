#ifndef ThreadControl_included
#define ThreadControl_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <numa.h>


#include "Global.h"
#include "Complex.h"
#include "DistributedMemoryObject.h"
#include "ExtremeFFT4D.h"
#include "Tools.h"
#include "xSSE.h"

#include FFTWIncludeFile

#define ThreadControl_SYNCDELAYLOOPS 100
#define ThreadControl_FFTPlanDBMax 1000
#define ThreadControl_WorkingSpaceSize 10000
#define ThreadControl_OutputSpaceSize  10000
#define ThreadControlParameterMAX 15
#define ThreadControl_ExecutionCommand_terminateThreads -1
#define ThreadControl_ExecutionCommand_allocateDistibutedMemory 1
#define ThreadControl_ExecutionCommand_scalarProductOfFermionVectors 2
#define ThreadControl_ExecutionCommand_vectorNormOfFermionVector 3
#define ThreadControl_ExecutionCommand_vectorAdditionOfFermionVectors 4
#define ThreadControl_ExecutionCommand_perform_FFTWFourierTransformationOfFermionVector 5
#define ThreadControl_ExecutionCommand_copyFermionVector 6
#define ThreadControl_ExecutionCommand_perform_PiModeRemoverOperatorMultiplicationInFourierSpace 7
#define ThreadControl_ExecutionCommand_perform_DiracTypeOperatorMultiplicationInFourierSpace 8
#define ThreadControl_ExecutionCommand_perform_yBOperatorMultiplication 9
#define ThreadControl_ExecutionCommand_perform_yBDaggerOperatorMultiplication 10
#define ThreadControl_ExecutionCommand_perform_RPreconditionerTypeMatrixMultiplicationInFourierSpace 11
#define ThreadControl_ExecutionCommand_perform_MultiplicationVectorWithDerivativesOfB 12
#define ThreadControl_ExecutionCommand_multiplyFermionVectorWithComplexNumber 13
#define ThreadControl_ExecutionCommand_zeroFermionVector 14
#define ThreadControl_ExecutionCommand_zeroUniformlyDistributedComplexVector 15 
#define ThreadControl_ExecutionCommand_vectorAdditionOfUniformlyDistributedComplexVectors 16
#define ThreadControl_ExecutionCommand_perform_xFFTFourierTransformationOfFermionVector 17
#define ThreadControl_ExecutionCommand_tune_xFFTFourierTransformationOfFermionVector 18
#define ThreadControl_ExecutionCommand_setFFTPlan_xFFTFourierTransformationOfFermionVector 19
#define ThreadControl_ExecutionCommand_perform_MultiCombinedOperator_VA1_RorCopy_VAc_RorCopy_D_InFourierSpace 20
#define ThreadControl_ExecutionCommand_multiplyFermionVectorWithTimeIndexedComplexScalars 21


struct ThreadControl_FFTPlanDB_Struct {
  Complex* input;
  Complex* output;
  int L0, L1, L2, L3;
  int blockSize;
  bool Forward;
  bool intertwined;
  fftw_plan fftwPlan;
};


class ThreadControl {
private:  
  pthread_t pThreadObject;
  pthread_mutex_t* pThreadMutex;
  int ThreadID;
  int NodeID;
  int CPUAllocationMode;
  long int ExecutionFlag;
  long int TerminationFlag;
  void** ExecutionParameters;
  int ExecutionCommand;
  Complex* workingSpace;
  Complex* outputSpace;
  ThreadControl_FFTPlanDB_Struct FFTPlanDB[ThreadControl_FFTPlanDBMax];
  int FFTPlanDBTotalEntriesCount;
  int FFTPlanDBIndexOfNextEntry;
  xSSE* xSSEObj;
  ExtremeFFT4D* xFFT;
  
  
  void allocateDistributedMemoryObject();
  void scalarProductOfFermionVectors();
  void vectorNormOfFermionVector();

  void zeroFermionVector();
  void zeroUniformlyDistributedComplexVector();
  
  
  void vectorAdditionOfFermionVectors();

  void vectorAdditionOfUniformlyDistributedComplexVectors();
  
  
  void perform_FFTWFourierTransformationOfFermionVector();
  void perform_xFFTFourierTransformationOfFermionVector();
  void tune_xFFTFourierTransformationOfFermionVector();
  void setFFTPlan_xFFTFourierTransformationOfFermionVector();
  fftw_plan getFFTPlan(Complex* input, Complex* output, int L0, int L1, int L2, int L3, int blockSize, bool forw, bool intertwined);
  
  void copyFermionVector();
  void multiplyFermionVectorWithComplexNumber();
  void multiplyFermionVectorWithTimeIndexedComplexScalars();
  void perform_PiModeRemoverOperatorMultiplicationInFourierSpace();
  void perform_DiracTypeOperatorMultiplicationInFourierSpace();
  void perform_yBOperatorMultiplication();
  void perform_yBDaggerOperatorMultiplication();
  void perform_RPreconditionerTypeMatrixMultiplicationInFourierSpace();
  void perform_MultiplicationVectorWithDerivativesOfB();
  void perform_MultiCombinedOperator_VA1_RorCopy_VAc_RorCopy_D_InFourierSpace();
  
  void resetFFTPlanEntry(int nr);
  
public:    
  ThreadControl(int threadid, int nodeid, int cpuAllocMode); 
  virtual ~ThreadControl();  

  void setThreadObject(pthread_t threadObj, pthread_mutex_t* mutObj);
  pthread_t getThreadObject();
  int getThreadID();
  int getNodeID();
  int getCPUAllocationMode();
  long int getExecutionFlag(xSSE* xSSEObjLoc);    
  long int getTerminationFlag(xSSE* xSSEObjLoc);  
  void setExecutionFlag(long int flag, xSSE* xSSEObjLoc);      
  void setTerminationFlag(long int flag, xSSE* xSSEObjLoc);      
  void setExecutionParameter(int nr, void* para);
  void* getExecutionParameter(int nr);
  void setExecutionCommand(int cmd);
  int getExecutionCommand();
  
  
  void execute();

};

#endif
