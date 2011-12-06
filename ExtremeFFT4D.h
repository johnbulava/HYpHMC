#ifndef ExtremeFFT4D_included
#define ExtremeFFT4D_included

#include <math.h>
#include <pthread.h>
#include "Global.h"
#include "Complex.h"
#include "Tools.h"
#include "ExtremeFFT4D_FFTstep.h"
#include "ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2.h"
#include "ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2.h"
#include "ExtremeFFT4D_FFTstep_Write2DSlice.h"
#include "ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2.h"
#include "ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2.h"


#define ExtremeFFT4D_TuneLevel_None 0
#define ExtremeFFT4D_TuneLevel_Low 1
#define ExtremeFFT4D_TuneLevel_Medium 2
#define ExtremeFFT4D_TuneLevel_High 3



#define ExtremeFFT4D_Forward true
#define ExtremeFFT4D_Backward false
#define ExtremeFFT4D_PrimeFactorsMAX 10
#define ExtremeFFT4D_PrefetchCombine1MAX 2
#define ExtremeFFT4D_PrefetchCombine2MAX 2
#define ExtremeFFT4D_ThreadControlLoopWaitCycles 100
#define ExtremeFFT4D_InternalEmbeddingOneMAX 8
#define ExtremeFFT4D_ThreadControlCommand_StopThread -1
#define ExtremeFFT4D_ThreadControlCommand_NoCommand 0
#define ExtremeFFT4D_ThreadControlCommand_PlanExecution 1

#define ExtremeFFT4D_InnerLoopMAX 32




class ExtremeFFT4D {
private:
  int L[4];
  int LargestL;
  int localIndexCount;
  int xtrSize1;
  int xtrSize2;
  int xtrSize3;
  int numberOfThreads;
  long int* CoreMaskForThreads;
  bool runThreadsContinously;

  int primeFactorCountL[4];
  int* primeFactorsL[4];
  int addressIncrementsInBytesL[4];
  Complex BitMask_ComplexI1;
  Complex BitMask_ComplexI2;
  
  ExtremeFFT4D_FFTstep** availFFTsteps;
  int availFFTstepCount;

  struct FFTplanType {
    bool complete;
    int LindA1;
    int LindA2;
    int LindB1;
    int LindB2;
    int InnerLoopCountA;
    int InnerLoopCountB;
    int OuterLoopCombineFactorA2;
    int OuterLoopCombineFactorB2; 
    int PrefetchCombineA1;   
    int PrefetchCombineA2;   
    int PrefetchCombineB1;   
    int PrefetchCombineB2; 
    long int inputSlice2DDistanceInBytesA;
    long int inputSlice2DDistanceInBytesB; 
    int MinorMissesA;
    int MajorMissesA;
    int MinorMissesB;
    int MajorMissesB;
    int InternalEmbeddingAOne;
    int InternalEmbeddingATwo;
    int InternalEmbeddingBOne;
    int InternalEmbeddingBTwo;          
    bool blockSimReadReadA;
    bool blockSimReadWriteA;
    bool blockSimWriteWriteA;    
    bool blockSimReadReadB;
    bool blockSimReadWriteB;
    bool blockSimWriteWriteB;        
    int FFTstepCountA;
    int FFTstepCountB;
    int FFTstepNrA[2*ExtremeFFT4D_PrimeFactorsMAX+2];
    int FFTstepNrB[2*ExtremeFFT4D_PrimeFactorsMAX+2];      
    ExtremeFFT4D_FFTstep*** FFTstepsA;
    ExtremeFFT4D_FFTstep*** FFTstepsB;    
    double performance_TotalCycles;
    double** performance_FFTstepCyclesA;
    double** performance_FFTstepCyclesB;    
    FFTplanType* next;
  };
  FFTplanType* possibleFFTplanList;
  FFTplanType* selectedFFTplan;
  

  struct ThreadControlDataStructureType {
    long int running;
    long int command;
  };
  ThreadControlDataStructureType* ThreadControlDataStructure;
  long int** ThreadControlParameterHolder;

  struct PlanExecutionControlDataStructureType {
    Complex* input;
    Complex* output;
    FFTplanType* plan;
    long int forward;
    long int doTiming;
    long int ReadFlag;
    long int WriteFlag;
    long int* SliceStatus;
    long int* WaitLoopEntered;
    long int readyTrafoA;
    long int readyTrafoB;    
  };
  PlanExecutionControlDataStructureType PlanExecutionControlDataStructure;

  
  xSSE* xSSEObj;
  Complex** Slices2D;

  

  void generateAllPossibleFFTplans(FFTplanType* &FFTplan, int state);
  void writeFFTplan(FFTplanType* fftPlan, bool writeSubsequent);
  void elaboratePlanDetails(FFTplanType* plan);
  void dismissPlanDetails(FFTplanType* plan, bool onlyFFTsteps);
  void optimizeInternalEmbeddingForPlan(FFTplanType* plan, int tuneLevel, bool enforceOpti);
  void getPlanReadyForExecution(FFTplanType* plan);
  
  void executePlan(Complex* input, Complex* output, FFTplanType* plan, int timingIterations, bool forward);
  void threadedExecutionOfPlan(int ThreadID, xSSE* xSSEObjLoc);
  void stopThreads();
  
  

  

public:
  ExtremeFFT4D(int l0, int l1, int l2, int l3, int xtrS1, int xtrS2, int xtrS3, int locIndCnt, int nrThreads, long int* threadCoreMsk, bool runContin); 
  ~ExtremeFFT4D(); 
  void ThreadControlLoopRoutine(int ThreadID);     //User never calls this routine !!!
  
  
  void startThreads();
  void stopThreadsAndWaitForTermination();
  void tune(int tuneLevel);
  void tune(int tuneLevel, Complex* input, Complex* output);
  void FastFourierTrafo(Complex* input, Complex* output, bool fft_forward);
  
  
  int getL0();
  int getL1();
  int getL2();
  int getL3();
  int getLocalIndexCount();
  int getXtraSize1();
  int getXtraSize2();
  int getXtraSize3();
  int getNumberOfThreads();
  char* getSelectedPlanDescriptor();
  bool setSelectedPlanFromDescriptor(char* descriptor);
};

#endif
