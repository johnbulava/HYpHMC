#include "MultiThreadedOperations.h"

ThreadControl* MultiThreadedOperationsThreadMainLoopThreadControlObject;
int MultiThreadedOperationsThreadMainLoopThreadID;
int MultiThreadedOperationsThreadMainLoopNodeID;
int MultiThreadedOperationsThreadMainLoopCPUAllocationMode;
void* MultiThreadedOperationsThreadMainLoopRoutine(void* para) {
  int threadID = MultiThreadedOperationsThreadMainLoopThreadID;
  int nodeID = MultiThreadedOperationsThreadMainLoopNodeID;
  int CPUallocMode = MultiThreadedOperationsThreadMainLoopCPUAllocationMode;
  
  if (LogLevel>2) printf("ThreadMainLoop called with ThreadID=%d, NodeID=%d, and CPUallocMode=%d\n",threadID,nodeID, CPUallocMode);
  if (CPUallocMode==1) {
    if (LogLevel>2) printf("Set Thread %d Affinity to 0-th Core on Node %d...",threadID,nodeID);  
    setCurrentThreadAffinityToNthCoreOnNodeID(nodeID,0);
  } else {
    if (LogLevel>2) printf("Set Thread %d Affinity to Node %d...",threadID,nodeID);  
    setCurrentThreadAffinityToNodeID(nodeID);  
  }

  double t = zeitwert();
  //wait for some seconds
  while (zeitwert()-t < MultiThreadedOperationsThreadMainLoopSecurityWaitTimeInSecs) {
    int a = 0;
    for (int I=0; I<MultiThreadedOperationsThreadMainLoopWaitCycles; I++) {
      a++;
    }
  }
  if (LogLevel>2) printf("successfully.\n");
  
  ThreadControl* threadControl = new ThreadControl(threadID, nodeID, CPUallocMode); 
  MultiThreadedOperationsThreadMainLoopThreadControlObject = threadControl;
  xSSE* xSSEObjLoc = new xSSE();
  
  if (LogLevel>2) printf("->Entering main loop of ThreadMainLoop with ThreadID=%d and NodeID=%d\n",threadID,nodeID);
  while (true) {
    long int b1 = threadControl->getTerminationFlag(xSSEObjLoc); 
    if (b1 == 1) break;
    while (true) {
      long int b2 = threadControl->getExecutionFlag(xSSEObjLoc); 
      if (b2 == 1) break;
        
      xSSEObjLoc->xSSE_PerformNopLoop_Wrapper(MultiThreadedOperationsThreadMainLoopWaitCycles);
    }

    threadControl->execute();
  }

  delete xSSEObjLoc;
  delete threadControl;
  if (LogLevel>2) printf("Leaving main loop of ThreadMainLoop with ThreadID=%d and NodeID=%d\n",threadID,nodeID);
  return NULL;
}


MultiThreadedOperations::MultiThreadedOperations(int opMode, bool enforceCoreTie) {
  if (LogLevel>2) printf("MultiThreadedOperations initializing with opMode=%d and enforceCoreTie=%d...\n",opMode, enforceCoreTie);

  initializeNUMASupport();
  xSSEObj = new xSSE();
  generalWorkSpace = createSuperAlignedComplex(MultiThreadedOperationsGeneralWorkSpaceSize);
  
  OperationMode = opMode;
  if (OperationMode<0) {
    printf("ERROR in MultiThreadedOperations: opMode=%d invalid\n",OperationMode);
    exit(0);
  }  
  if (OperationMode>2) {
    printf("ERROR in MultiThreadedOperations: opMode=%d invalid\n",OperationMode);
    exit(0);
  }  
  
  if (OperationMode==0) ThreadCount = 1;
  if (OperationMode==1) ThreadCount = 2;
  if (OperationMode==2) ThreadCount = 4;
    
  if (LogLevel>1) printf("running on %d Threads in total\n",ThreadCount);
  Threads = new ThreadControl*[ThreadCount];
  for (int I=0; I<ThreadCount; I++) {
    Threads[I] = NULL;
  }
  
  int succ = pthread_mutex_init (&pThreadMutex, NULL);
  if (succ != 0) {
    printf("ERROR in MultiThreadedOperations: Could not initialize Mutex!\n");
    exit(0);
  }
  
  
  int nodeID = 0;
  MultiThreadedOperationsThreadMainLoopCPUAllocationMode = 1;
  if (!enforceCoreTie) MultiThreadedOperationsThreadMainLoopCPUAllocationMode = 2;
  if (ThreadCount>numaNodesCount) MultiThreadedOperationsThreadMainLoopCPUAllocationMode = 2;  
  for (int I=1; I<ThreadCount; I++) {
    nodeID++;
    if (nodeID>=numaNodesCount) nodeID = 0;
    pthread_t thread;

    MultiThreadedOperationsThreadMainLoopThreadControlObject = NULL;
    if (LogLevel>2) printf("Trying to start thread nr %d on node %d...",I, nodeID);
    MultiThreadedOperationsThreadMainLoopThreadID = I;
    MultiThreadedOperationsThreadMainLoopNodeID = nodeID;
    if (pthread_create(&thread, NULL, MultiThreadedOperationsThreadMainLoopRoutine, NULL) != 0) {
      printf("ERROR in MultiThreadedOperations: Could not create thread nr %d\n",I);
      exit(0);
    }    

    ThreadControl*** dummy = new ThreadControl**[500];  //Strange construction to prevent compiler-optimization !!!
    for (int I2=0; I2<500; I2++) dummy[I2] = &MultiThreadedOperationsThreadMainLoopThreadControlObject;
    bool leaveLoop = false;
    while (!leaveLoop) {
      long int f = (long int) MultiThreadedOperationsThreadMainLoopThreadControlObject;
      int a = (int) f;
      for (int I2=0; I2<MultiThreadedOperationsThreadMainLoopWaitCycles; I2++) {
        a++;
      }
      for (int I2=0; I2<500; I2++) if ((*(dummy[I2])) != NULL) leaveLoop = true;
    }
    if (LogLevel>2) printf("successfully\n");
    Threads[I] = (*(dummy[0]));
    Threads[I]->setThreadObject(thread, &pThreadMutex);
    delete[] dummy;
  }
  
  if (MultiThreadedOperationsThreadMainLoopCPUAllocationMode == 1) {
    if (LogLevel>2) printf("Set Thread 0 Affinity to 0-th core on Node 0...");
    setCurrentThreadAffinityToNthCoreOnNodeID(0, 0);
  } else {
    if (LogLevel>2) printf("Set Thread 0 Affinity to Node 0...");
    setCurrentThreadAffinityToNodeID(0);  
  }
  double t = zeitwert();
  //wait for some seconds
  while (zeitwert()-t < MultiThreadedOperationsThreadMainLoopSecurityWaitTimeInSecs) {
    int a = 0;
    for (int I=0; I<MultiThreadedOperationsThreadMainLoopWaitCycles; I++) {
      a++;
    }
  }  
  Threads[0] = new ThreadControl(0, 0, MultiThreadedOperationsThreadMainLoopCPUAllocationMode); 
  Threads[0]->setThreadObject(pthread_self(), &pThreadMutex);
  if (LogLevel>2) printf("successfully.\n");  
  
  
  threadCompleteMask = 1;
  for (int I=0; I<ThreadCount; I++) {
    threadCompleteMask *= 2;
  }
  threadCompleteMask -= 1;
  
  xSSEObj->xSSE_PerformNopLoop_Wrapper(MultiThreadedOperationsThreadMainLoopWaitCycles);
  
  if (LogLevel>2) printf("MultiThreadedOperations initialized with opMode=%d successfully\n",opMode);
}


MultiThreadedOperations::~MultiThreadedOperations() {	
  terminateThreads();  
  delete Threads[0];
  delete[] Threads;
  pthread_mutex_destroy(&pThreadMutex);
  delete xSSEObj;
  destroySuperAlignedComplex(generalWorkSpace);
}


void MultiThreadedOperations::threadedExecute(int threadMask) {
  int dummy = threadMask;
  int thCnt = 0;
  if (LogLevel>4) printf("threadedExecute called with threadMask=%d\n",threadMask);
  while (dummy>0) {
    if ((dummy % 2) == 1) {
      if (thCnt>=ThreadCount) {
        break;
      }
      Threads[thCnt]->setExecutionFlag(1, xSSEObj); 
    }
    thCnt++;
    dummy /= 2;  
  }
  if ((threadMask % 2) == 1) {
    Threads[0]->execute();
  }
  
  if (ThreadCount>1) {
    bool stillRunning = true;
    while (stillRunning) {
      stillRunning = false;
      for (int I=0; I<ThreadCount; I++) {
        long int b = Threads[I]->getExecutionFlag(xSSEObj);
        if (b == 1) stillRunning = true;      
      }

      xSSEObj->xSSE_PerformNopLoop_Wrapper(MultiThreadedOperationsThreadMainLoopWaitCycles);
    }
  }

  if (LogLevel>4) printf("threadedExecute finished\n");
}


inline void MultiThreadedOperations::setParameterOnAllThreads(int nr, void* para) {
  for (int I=0; I<ThreadCount; I++) {
    Threads[I]->setExecutionParameter(nr, para);
  }
}


inline void MultiThreadedOperations::setExecutionCommandOnAllThreads(int cmd) {
  for (int I=0; I<ThreadCount; I++) {
    Threads[I]->setExecutionCommand(cmd);
  }
}


inline Complex MultiThreadedOperations::getComplexNumberSumFromAllThreads(int nr) {
  Complex sum(0,0);
  for (int I=0; I<ThreadCount; I++) {
    Complex* c = (Complex*) Threads[I]->getExecutionParameter(nr);
    sum = sum + (*c);
  }
  return sum;
}


void MultiThreadedOperations::terminateThreads() {
  setExecutionCommandOnAllThreads(ThreadControl_ExecutionCommand_terminateThreads);
  for (int I=0; I<ThreadCount; I++) {
    Threads[I]->setExecutionFlag(1, xSSEObj); 
  }  
  Threads[0]->execute();
  
  double t = zeitwert();
  //wait for some seconds
  while (zeitwert()-t < MultiThreadedOperationsThreadMainLoopSecurityWaitTimeInSecs) {
    xSSEObj->xSSE_PerformNopLoop_Wrapper(MultiThreadedOperationsThreadMainLoopWaitCycles);
  }  
}


DistributedMemoryObject* MultiThreadedOperations::allocateDistibutedMemory(int sizeInComplexNumbers) {
  DistributedMemoryObject* memObj = new DistributedMemoryObject(ThreadCount);
  setExecutionCommandOnAllThreads(ThreadControl_ExecutionCommand_allocateDistibutedMemory);
  setParameterOnAllThreads(0, (void*) memObj);
  int dummy = sizeInComplexNumbers;
  setParameterOnAllThreads(1, ((void*) (&dummy)));
  threadedExecute(threadCompleteMask);
  return memObj;
}


DistributedMemoryObject* MultiThreadedOperations::allocateDistibutedFermionVector(int L0, int L1, int L2, int L3) {
  int sizeInComplexNumbers = L0 * (L1+xtraSize1) * (L2+xtraSize2) * (L3+xtraSize3);
  if (OperationMode == 0) {
    sizeInComplexNumbers *= 8;
  }
  if (OperationMode == 1) {
    sizeInComplexNumbers *= 4;
  }
  if (OperationMode == 2) {
    sizeInComplexNumbers *= 2;
  }

  return allocateDistibutedMemory(sizeInComplexNumbers);
}


void MultiThreadedOperations::copyFermionVectorToDistributedFermionVector(Complex* vec, DistributedMemoryObject* memObj, int L0, int L1, int L2, int L3) {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  int i0,i1,i2,i3;
  Complex* working = new Complex[8];

  int countVec = 0;
  int countMem = 0;
  int xtrAdd1Vec = xtraSize1*8*(L3+xtraSize3)*(L2+xtraSize2);
  int xtrAdd2Vec = xtraSize2*8*(L3+xtraSize3);
  int xtrAdd3Vec = xtraSize3*8;

  int memStride = 8;
  if (OperationMode == 1) memStride = 4;
  if (OperationMode == 2) memStride = 2;
    
  int xtrAdd1Mem = xtraSize1*memStride*(L3+xtraSize3)*(L2+xtraSize2);
  int xtrAdd2Mem = xtraSize2*memStride*(L3+xtraSize3);
  int xtrAdd3Mem = xtraSize3*memStride;

  
  for (i0=0; i0<L0; i0++) {
    for (i1=0; i1<L1; i1++) {
      for (i2=0; i2<L2; i2++) {
        for (i3=0; i3<L3; i3++) {

          for (int I2=0; I2<8; I2++) {
            working[I2].x = vec[countVec+I2].x;
            working[I2].y = vec[countVec+I2].y;
          }	

          if (OperationMode == 0) {
	    Complex* v0 = memObj->dataSets[0];
            for (int I2=0; I2<8; I2++) {
	      v0[countMem+I2].x = working[I2].x;
	      v0[countMem+I2].y = working[I2].y;
            }	
          } else if (OperationMode == 1) {
	    Complex* v0 = memObj->dataSets[0];
	    Complex* v1 = memObj->dataSets[1];
            v0[countMem+0].x = working[0].x;
	    v0[countMem+0].y = working[0].y;
            v0[countMem+1].x = working[1].x;
	    v0[countMem+1].y = working[1].y;

            v1[countMem+0].x = working[2].x;
	    v1[countMem+0].y = working[2].y;
            v1[countMem+1].x = working[3].x;
	    v1[countMem+1].y = working[3].y;

            v0[countMem+2].x = working[4].x;
	    v0[countMem+2].y = working[4].y;
            v0[countMem+3].x = working[5].x;
	    v0[countMem+3].y = working[5].y;

            v1[countMem+2].x = working[6].x;
	    v1[countMem+2].y = working[6].y;
            v1[countMem+3].x = working[7].x;
	    v1[countMem+3].y = working[7].y;
	  } else if (OperationMode == 2) {
	    Complex* v0 = memObj->dataSets[0];
	    Complex* v1 = memObj->dataSets[1];
	    Complex* v2 = memObj->dataSets[2];
	    Complex* v3 = memObj->dataSets[3];
            v0[countMem+0].x = working[0].x;
	    v0[countMem+0].y = working[0].y;
            v1[countMem+0].x = working[1].x;
	    v1[countMem+0].y = working[1].y;
            v2[countMem+0].x = working[2].x;
	    v2[countMem+0].y = working[2].y;
            v3[countMem+0].x = working[3].x;
	    v3[countMem+0].y = working[3].y;

            v0[countMem+1].x = working[4].x;
	    v0[countMem+1].y = working[4].y;
            v1[countMem+1].x = working[5].x;
	    v1[countMem+1].y = working[5].y;
            v2[countMem+1].x = working[6].x;
	    v2[countMem+1].y = working[6].y;
            v3[countMem+1].x = working[7].x;
	    v3[countMem+1].y = working[7].y;
	  }

          countVec += 8;
          countMem += memStride;
	}
	countVec += xtrAdd3Vec;
        countMem += xtrAdd3Mem;
      }
      countVec += xtrAdd2Vec;
      countMem += xtrAdd2Mem;
    }
    countVec += xtrAdd1Vec;
    countMem += xtrAdd1Mem;
  }
  delete[] working;
  addPerformanceProfilingItem("MultiThreadedOperations::copyFermionVectorToDistributedFermionVector", performanceProfilerStartCycle, 0);
}


void MultiThreadedOperations::copyDistributedFermionVectorToFermionVector(DistributedMemoryObject* memObj, Complex* vec, int L0, int L1, int L2, int L3) {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  int i0,i1,i2,i3;
  Complex* working = new Complex[8];

  int countVec = 0;
  int countMem = 0;
  int xtrAdd1Vec = xtraSize1*8*(L3+xtraSize3)*(L2+xtraSize2);
  int xtrAdd2Vec = xtraSize2*8*(L3+xtraSize3);
  int xtrAdd3Vec = xtraSize3*8;

  int memStride = 8;
  if (OperationMode == 1) memStride = 4;
  if (OperationMode == 2) memStride = 2;
    
  int xtrAdd1Mem = xtraSize1*memStride*(L3+xtraSize3)*(L2+xtraSize2);
  int xtrAdd2Mem = xtraSize2*memStride*(L3+xtraSize3);
  int xtrAdd3Mem = xtraSize3*memStride;

  
  for (i0=0; i0<L0; i0++) {
    for (i1=0; i1<L1; i1++) {
      for (i2=0; i2<L2; i2++) {
        for (i3=0; i3<L3; i3++) {

          if (OperationMode == 0) {
	    Complex* v0 = memObj->dataSets[0];
            for (int I2=0; I2<8; I2++) {
	      working[I2].x = v0[countMem+I2].x;
	      working[I2].y = v0[countMem+I2].y;
            }	
          } else if (OperationMode == 1) {
	    Complex* v0 = memObj->dataSets[0];
	    Complex* v1 = memObj->dataSets[1];
            working[0].x = v0[countMem+0].x;
	    working[0].y = v0[countMem+0].y;
            working[1].x = v0[countMem+1].x;
	    working[1].y = v0[countMem+1].y;

            working[2].x = v1[countMem+0].x;
	    working[2].y = v1[countMem+0].y;
            working[3].x = v1[countMem+1].x;
	    working[3].y = v1[countMem+1].y;

            working[4].x = v0[countMem+2].x;
	    working[4].y = v0[countMem+2].y;
            working[5].x = v0[countMem+3].x;
	    working[5].y = v0[countMem+3].y;

            working[6].x = v1[countMem+2].x;
	    working[6].y = v1[countMem+2].y;
            working[7].x = v1[countMem+3].x;
	    working[7].y = v1[countMem+3].y;
	  } else if (OperationMode == 2) {
	    Complex* v0 = memObj->dataSets[0];
	    Complex* v1 = memObj->dataSets[1];
	    Complex* v2 = memObj->dataSets[2];
	    Complex* v3 = memObj->dataSets[3];
            working[0].x = v0[countMem+0].x;
	    working[0].y = v0[countMem+0].y;
            working[1].x = v1[countMem+0].x;
	    working[1].y = v1[countMem+0].y;
            working[2].x = v2[countMem+0].x;
	    working[2].y = v2[countMem+0].y;
            working[3].x = v3[countMem+0].x;
	    working[3].y = v3[countMem+0].y;

            working[4].x = v0[countMem+1].x;
	    working[4].y = v0[countMem+1].y;
            working[5].x = v1[countMem+1].x;
	    working[5].y = v1[countMem+1].y;
            working[6].x = v2[countMem+1].x;
	    working[6].y = v2[countMem+1].y;
            working[7].x = v3[countMem+1].x;
	    working[7].y = v3[countMem+1].y;
	  }

          for (int I2=0; I2<8; I2++) {
            vec[countVec+I2].x = working[I2].x;
            vec[countVec+I2].y = working[I2].y;
          }	

          countVec += 8;
          countMem += memStride;
	}
	countVec += xtrAdd3Vec;
        countMem += xtrAdd3Mem;
      }
      countVec += xtrAdd2Vec;
      countMem += xtrAdd2Mem;
    }
    countVec += xtrAdd1Vec;
    countMem += xtrAdd1Mem;
  }
  delete[] working;
  addPerformanceProfilingItem("MultiThreadedOperations::copyDistributedFermionVectorToFermionVector", performanceProfilerStartCycle, 0);
}


void MultiThreadedOperations::copyComplexVectorToUniformlyDistributedMemory(Complex* vec, DistributedMemoryObject* memObj, int size) {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  for (int I2=0; I2<memObj->dataSetCount; I2++) {
    Complex* output = memObj->dataSets[I2];

    for (int I=0; I<size; I++) {
      output[I].x = vec[I].x;
      output[I].y = vec[I].y;
    }
  }
  addPerformanceProfilingItem("MultiThreadedOperations::copyComplexVectorToUniformlyDistributedMemory", performanceProfilerStartCycle, 0);
}


void MultiThreadedOperations::calcSumOfUniformlyDistributedComplexVectors(DistributedMemoryObject* memObj, Complex* vec, int size) {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  for (int I=0; I<size; I++) {
    vec[I].x = 0;
    vec[I].y = 0;
  }
  for (int I2=0; I2<memObj->dataSetCount; I2++) {
    Complex* v = memObj->dataSets[I2];
    for (int I=0; I<size; I++) {
      vec[I].x += v[I].x;
      vec[I].y += v[I].y;
    }
  }
  addPerformanceProfilingItem("MultiThreadedOperations::calcSumOfUniformlyDistributedComplexVectors", performanceProfilerStartCycle, 0);
}


void MultiThreadedOperations::zeroFermionVector(DistributedMemoryObject* memObj, int L0, int L1, int L2, int L3) {
  setExecutionCommandOnAllThreads(ThreadControl_ExecutionCommand_zeroFermionVector);
  setParameterOnAllThreads(0, (void*) memObj);
  
  int blockSize = 8;
  if (OperationMode == 1) blockSize = 4;
  if (OperationMode == 2) blockSize = 2;  
  
  int intParas[5];
  intParas[0] = L0;
  intParas[1] = L1;
  intParas[2] = L2;
  intParas[3] = L3;  
  intParas[4] = blockSize;    
  setParameterOnAllThreads(2, ((void*) (&(intParas[0]))));

  threadedExecute(threadCompleteMask);
}


void MultiThreadedOperations::zeroUniformlyDistributedComplexVector(DistributedMemoryObject* memObj, int size) {
  setExecutionCommandOnAllThreads(ThreadControl_ExecutionCommand_zeroUniformlyDistributedComplexVector);
  setParameterOnAllThreads(0, (void*) memObj);
  
  setParameterOnAllThreads(2, ((void*) (&(size))));
  
  threadedExecute(threadCompleteMask);
}


void MultiThreadedOperations::vectorAdditionOfUniformlyDistributedComplexVectors(DistributedMemoryObject* x, DistributedMemoryObject* y, Complex& alpha, int size) {
  setExecutionCommandOnAllThreads(ThreadControl_ExecutionCommand_vectorAdditionOfUniformlyDistributedComplexVectors);
  setParameterOnAllThreads(0, (void*) x);
  setParameterOnAllThreads(1, (void*) y);
  
  setParameterOnAllThreads(2, ((void*) (&(size))));

  setParameterOnAllThreads(3, ((void*) (&(alpha))));
  
  threadedExecute(threadCompleteMask);
}


void MultiThreadedOperations::vectorAdditionOfFermionVectors(DistributedMemoryObject* x, DistributedMemoryObject* y, Complex& alpha, int L0, int L1, int L2, int L3) {
  setExecutionCommandOnAllThreads(ThreadControl_ExecutionCommand_vectorAdditionOfFermionVectors);
  setParameterOnAllThreads(0, (void*) x);
  setParameterOnAllThreads(1, (void*) y);
  
  int blockSize = 8;
  if (OperationMode == 1) blockSize = 4;
  if (OperationMode == 2) blockSize = 2;  
  
  int intParas[5];
  intParas[0] = L0;
  intParas[1] = L1;
  intParas[2] = L2;
  intParas[3] = L3;  
  intParas[4] = blockSize;    
  setParameterOnAllThreads(2, ((void*) (&(intParas[0]))));

  setParameterOnAllThreads(3, ((void*) (&(alpha))));
  
  threadedExecute(threadCompleteMask);
}


void MultiThreadedOperations::perform_FFTWFourierTransformationOfFermionVector(DistributedMemoryObject* input, DistributedMemoryObject* output, bool forw, int L0, int L1, int L2, int L3) {
  setExecutionCommandOnAllThreads(ThreadControl_ExecutionCommand_perform_FFTWFourierTransformationOfFermionVector);
  setParameterOnAllThreads(0, (void*) input);
  setParameterOnAllThreads(1, (void*) output);
  
  int blockSize = 8;
  if (OperationMode == 1) blockSize = 4;
  if (OperationMode == 2) blockSize = 2;  
  
  int intParas[7];
  intParas[0] = L0;
  intParas[1] = L1;
  intParas[2] = L2;
  intParas[3] = L3;
  intParas[4] = forw;
  intParas[5] = blockSize;    
  intParas[6] = 1;       // intertwined?  
  setParameterOnAllThreads(2, ((void*) (&(intParas[0]))));
  
  threadedExecute(threadCompleteMask);
}


void MultiThreadedOperations::perform_xFFTFourierTransformationOfFermionVector(DistributedMemoryObject* input, DistributedMemoryObject* output, bool forw, int threadCount, int L0, int L1, int L2, int L3) {
  setExecutionCommandOnAllThreads(ThreadControl_ExecutionCommand_perform_xFFTFourierTransformationOfFermionVector);
  setParameterOnAllThreads(0, (void*) input);
  setParameterOnAllThreads(1, (void*) output);
  
  int blockSize = 8;
  if (OperationMode == 1) blockSize = 4;
  if (OperationMode == 2) blockSize = 2;  
  
  int intParas[8];
  intParas[0] = L0;
  intParas[1] = L1;
  intParas[2] = L2;
  intParas[3] = L3;
  intParas[4] = forw;
  intParas[5] = blockSize;    
  intParas[6] = 1;       // intertwined?  
  intParas[7] = threadCount;       
  setParameterOnAllThreads(2, ((void*) (&(intParas[0]))));
  
  threadedExecute(threadCompleteMask);
}


void MultiThreadedOperations::scalarProductOfFermionVectors(DistributedMemoryObject* memObj1, DistributedMemoryObject* memObj2, Complex& res, int L0, int L1, int L2, int L3) {
  setExecutionCommandOnAllThreads(ThreadControl_ExecutionCommand_scalarProductOfFermionVectors);
  setParameterOnAllThreads(0, (void*) memObj1);
  setParameterOnAllThreads(1, (void*) memObj2);
  
  int blockSize = 8;
  if (OperationMode == 1) blockSize = 4;
  if (OperationMode == 2) blockSize = 2;  
  
  int intParas[5];
  intParas[0] = L0;
  intParas[1] = L1;
  intParas[2] = L2;
  intParas[3] = L3;  
  intParas[4] = blockSize;    
  setParameterOnAllThreads(2, ((void*) (&(intParas[0]))));
  
  threadedExecute(threadCompleteMask);
  
  res = getComplexNumberSumFromAllThreads(3);
}


void MultiThreadedOperations::multiplyFermionVectorWithComplexNumber(DistributedMemoryObject* input, DistributedMemoryObject* output, Complex alpha, int L0, int L1, int L2, int L3) {
  setExecutionCommandOnAllThreads(ThreadControl_ExecutionCommand_multiplyFermionVectorWithComplexNumber);
  setParameterOnAllThreads(0, (void*) input);
  setParameterOnAllThreads(1, (void*) output);
  
  int blockSize = 8;
  if (OperationMode == 1) blockSize = 4;
  if (OperationMode == 2) blockSize = 2;  
  
  int intParas[5];
  intParas[0] = L0;
  intParas[1] = L1;
  intParas[2] = L2;
  intParas[3] = L3;  
  intParas[4] = blockSize;    
  setParameterOnAllThreads(2, ((void*) (&(intParas[0]))));

  setParameterOnAllThreads(3, ((void*) (&(alpha))));
  
  threadedExecute(threadCompleteMask);
}

void MultiThreadedOperations::multiplyFermionVectorWithTimeIndexedComplexScalars(DistributedMemoryObject* input, DistributedMemoryObject* output, DistributedMemoryObject* scalars, int L0, int L1, int L2, int L3) {
  setExecutionCommandOnAllThreads(ThreadControl_ExecutionCommand_multiplyFermionVectorWithTimeIndexedComplexScalars);
  setParameterOnAllThreads(0, (void*) input);
  setParameterOnAllThreads(1, (void*) output);
  
  int blockSize = 8;
  if (OperationMode == 1) blockSize = 4;
  if (OperationMode == 2) blockSize = 2;  
  
  int intParas[5];
  intParas[0] = L0;
  intParas[1] = L1;
  intParas[2] = L2;
  intParas[3] = L3;  
  intParas[4] = blockSize;    
  setParameterOnAllThreads(2, ((void*) (&(intParas[0]))));

  setParameterOnAllThreads(3, ((void*) scalars));
  
  threadedExecute(threadCompleteMask);
}








void MultiThreadedOperations::vectorNormOfFermionVector(DistributedMemoryObject* memObj, double& res, int L0, int L1, int L2, int L3) {
  setExecutionCommandOnAllThreads(ThreadControl_ExecutionCommand_vectorNormOfFermionVector);
  setParameterOnAllThreads(0, (void*) memObj);

  int blockSize = 8;
  if (OperationMode == 1) blockSize = 4;
  if (OperationMode == 2) blockSize = 2;

  int intParas[5];
  intParas[0] = L0;
  intParas[1] = L1;
  intParas[2] = L2;
  intParas[3] = L3;  
  intParas[4] = blockSize;    
  setParameterOnAllThreads(2, ((void*) (&(intParas[0]))));
  
  threadedExecute(threadCompleteMask);
  
  res = (getComplexNumberSumFromAllThreads(3)).x;
}


void MultiThreadedOperations::copyFermionVector(DistributedMemoryObject* input, DistributedMemoryObject* output, int L0, int L1, int L2, int L3) {
  setExecutionCommandOnAllThreads(ThreadControl_ExecutionCommand_copyFermionVector);
  setParameterOnAllThreads(0, (void*) input);
  setParameterOnAllThreads(1, (void*) output);
  
  int blockSize = 8;
  if (OperationMode == 1) blockSize = 4;
  if (OperationMode == 2) blockSize = 2;  
  
  int intParas[5];
  intParas[0] = L0;
  intParas[1] = L1;
  intParas[2] = L2;
  intParas[3] = L3;  
  intParas[4] = blockSize;    
  setParameterOnAllThreads(2, ((void*) (&(intParas[0]))));

  threadedExecute(threadCompleteMask);
}


void MultiThreadedOperations::perform_PiModeRemoverOperatorMultiplicationInFourierSpace(DistributedMemoryObject* input, DistributedMemoryObject* output, int L0, int L1, int L2, int L3, DistributedMemoryObject* Index_PiModes) {
  setExecutionCommandOnAllThreads(ThreadControl_ExecutionCommand_perform_PiModeRemoverOperatorMultiplicationInFourierSpace);
  setParameterOnAllThreads(0, (void*) input);
  setParameterOnAllThreads(1, (void*) output);
  setParameterOnAllThreads(3, (void*) Index_PiModes);
  
  int blockSize = 8;
  if (OperationMode == 1) blockSize = 4;
  if (OperationMode == 2) blockSize = 2;  
  
  int intParas[5];
  intParas[0] = L0;
  intParas[1] = L1;
  intParas[2] = L2;
  intParas[3] = L3;  
  intParas[4] = blockSize;    
  setParameterOnAllThreads(2, ((void*) (&(intParas[0]))));
  
  threadedExecute(threadCompleteMask);  
}


void MultiThreadedOperations::perform_DiracTypeOperatorMultiplicationInFourierSpace(DistributedMemoryObject* input, DistributedMemoryObject* output, int L0, int L1, int L2, int L3, DistributedMemoryObject* sinP, DistributedMemoryObject* operatorData) {
  setExecutionCommandOnAllThreads(ThreadControl_ExecutionCommand_perform_DiracTypeOperatorMultiplicationInFourierSpace);
  setParameterOnAllThreads(0, (void*) input);
  setParameterOnAllThreads(1, (void*) output);
  setParameterOnAllThreads(3, (void*) sinP);
  setParameterOnAllThreads(4, (void*) operatorData);
  
  int blockSize = 8;
  if (OperationMode == 1) blockSize = 4;
  if (OperationMode == 2) blockSize = 2;  
  
  int intParas[5];
  intParas[0] = L0;
  intParas[1] = L1;
  intParas[2] = L2;
  intParas[3] = L3;  
  intParas[4] = blockSize;    
  setParameterOnAllThreads(2, ((void*) (&(intParas[0]))));
  
  int threadMask = threadCompleteMask;
  if (OperationMode == 2) threadMask = 3;
  
  threadedExecute(threadMask);    
}


void MultiThreadedOperations::perform_MultiCombinedOperator_VA1_RorCopy_VAc_RorCopy_D_InFourierSpace(DistributedMemoryObject* v1, DistributedMemoryObject* v2, DistributedMemoryObject* v3, DistributedMemoryObject* v4, int L0, int L1, int L2, int L3, bool useR, Complex alpha2, DistributedMemoryObject* sinP, DistributedMemoryObject* operatorData, DistributedMemoryObject* rPrecData1, DistributedMemoryObject* rPrecData2) {
  setExecutionCommandOnAllThreads(ThreadControl_ExecutionCommand_perform_MultiCombinedOperator_VA1_RorCopy_VAc_RorCopy_D_InFourierSpace);
  setParameterOnAllThreads(0, (void*) v1);
  setParameterOnAllThreads(1, (void*) v2);
  setParameterOnAllThreads(2, (void*) v3);
  setParameterOnAllThreads(3, (void*) v4);
  setParameterOnAllThreads(4, (void*) sinP);
  setParameterOnAllThreads(5, (void*) operatorData);
  setParameterOnAllThreads(6, (void*) rPrecData1);
  setParameterOnAllThreads(7, (void*) rPrecData2);
  
  int blockSize = 8;
  if (OperationMode == 1) blockSize = 4;
  if (OperationMode == 2) blockSize = 2;  
  
  int intParas[7];
  intParas[0] = L0;
  intParas[1] = L1;
  intParas[2] = L2;
  intParas[3] = L3;  
  intParas[4] = blockSize;    
  intParas[5] = useR;      
  intParas[6] = OperationMode;        
  setParameterOnAllThreads(8, ((void*) (&(intParas[0]))));
  setParameterOnAllThreads(9, ((void*) (&(alpha2))));
  long int* longP = (long int*) generalWorkSpace;  
  for (int I=0; I<MultiThreadedOperationsGeneralWorkSpaceSize; I++) {
    longP[I] = 0;
  }
  setParameterOnAllThreads(10, ((void*) generalWorkSpace));
  
  
  threadedExecute(threadCompleteMask);    
}


void MultiThreadedOperations::perform_yBOperatorMultiplication(DistributedMemoryObject* input, DistributedMemoryObject* output, int L0, int L1, int L2, int L3, double y, double split, DistributedMemoryObject* phiObj) {
  setExecutionCommandOnAllThreads(ThreadControl_ExecutionCommand_perform_yBOperatorMultiplication);
  setParameterOnAllThreads(0, (void*) input);
  setParameterOnAllThreads(1, (void*) output);
  setParameterOnAllThreads(3, (void*) phiObj);
  
  int blockSize = 8;
  if (OperationMode == 1) blockSize = 4;
  if (OperationMode == 2) blockSize = 2;  
  
  int intParas[5];
  intParas[0] = L0;
  intParas[1] = L1;
  intParas[2] = L2;
  intParas[3] = L3;  
  intParas[4] = blockSize;    
  setParameterOnAllThreads(2, ((void*) (&(intParas[0]))));
  
  double doubleParas[2];
  doubleParas[0] = y;
  doubleParas[1] = split;
  setParameterOnAllThreads(4, ((void*) (&(doubleParas[0]))));
  
  threadedExecute(threadCompleteMask);  
}


void MultiThreadedOperations::perform_yBDaggerOperatorMultiplication(DistributedMemoryObject* input, DistributedMemoryObject* output, int L0, int L1, int L2, int L3, double y, double split, DistributedMemoryObject* phiObj) {
  setExecutionCommandOnAllThreads(ThreadControl_ExecutionCommand_perform_yBDaggerOperatorMultiplication);
  setParameterOnAllThreads(0, (void*) input);
  setParameterOnAllThreads(1, (void*) output);
  setParameterOnAllThreads(3, (void*) phiObj);
  
  int blockSize = 8;
  if (OperationMode == 1) blockSize = 4;
  if (OperationMode == 2) blockSize = 2;  
  
  int intParas[5];
  intParas[0] = L0;
  intParas[1] = L1;
  intParas[2] = L2;
  intParas[3] = L3;  
  intParas[4] = blockSize;    
  setParameterOnAllThreads(2, ((void*) (&(intParas[0]))));

  double doubleParas[2];
  doubleParas[0] = y;
  doubleParas[1] = split;
  setParameterOnAllThreads(4, ((void*) (&(doubleParas[0]))));
  
  threadedExecute(threadCompleteMask);  
}


void MultiThreadedOperations::perform_RPreconditionerTypeMatrixMultiplicationInFourierSpace(DistributedMemoryObject* input, DistributedMemoryObject* output, int L0, int L1, int L2, int L3, DistributedMemoryObject* rPrecData) {
  setExecutionCommandOnAllThreads(ThreadControl_ExecutionCommand_perform_RPreconditionerTypeMatrixMultiplicationInFourierSpace);
  setParameterOnAllThreads(0, (void*) input);
  setParameterOnAllThreads(1, (void*) output);
  setParameterOnAllThreads(3, (void*) rPrecData);
  
  int blockSize = 8;
  if (OperationMode == 1) blockSize = 4;
  if (OperationMode == 2) blockSize = 2;  
  
  int intParas[5];
  intParas[0] = L0;
  intParas[1] = L1;
  intParas[2] = L2;
  intParas[3] = L3;  
  intParas[4] = blockSize;    
  setParameterOnAllThreads(2, ((void*) (&(intParas[0]))));
  
  threadedExecute(threadCompleteMask);  
}


void MultiThreadedOperations::perform_MultiplicationVectorWithDerivativesOfB(DistributedMemoryObject* leftInput, DistributedMemoryObject* rightInput, int L0, int L1, int L2, int L3, double split, DistributedMemoryObject* output) {
  setExecutionCommandOnAllThreads(ThreadControl_ExecutionCommand_perform_MultiplicationVectorWithDerivativesOfB);
  setParameterOnAllThreads(0, (void*) leftInput);
  setParameterOnAllThreads(1, (void*) rightInput);
  setParameterOnAllThreads(3, (void*) output);
  
  int blockSize = 8;
  if (OperationMode == 1) blockSize = 4;
  if (OperationMode == 2) blockSize = 2;  
  
  int intParas[5];
  intParas[0] = L0;
  intParas[1] = L1;
  intParas[2] = L2;
  intParas[3] = L3;  
  intParas[4] = blockSize;    
  setParameterOnAllThreads(2, ((void*) (&(intParas[0]))));

  double doubleParas[1];
  doubleParas[0] = split;
  setParameterOnAllThreads(4, ((void*) (&(doubleParas[0]))));
  
  threadedExecute(threadCompleteMask);  
}


void MultiThreadedOperations::tune_xFFTFourierTransformationOfFermionVector(DistributedMemoryObject* input, DistributedMemoryObject* output, int threadCount, int L0, int L1, int L2, int L3, int tuneLevel) {
  char* fftPlanDescriptor;
  tune_xFFTFourierTransformationOfFermionVector(input, output, threadCount, L0, L1, L2, L3, tuneLevel, fftPlanDescriptor);
  delete[] fftPlanDescriptor;
}


void MultiThreadedOperations::tune_xFFTFourierTransformationOfFermionVector(DistributedMemoryObject* input, DistributedMemoryObject* output, int threadCount, int L0, int L1, int L2, int L3, int tuneLevel, char* &fftPlanDescriptor) {
  setExecutionCommandOnAllThreads(ThreadControl_ExecutionCommand_tune_xFFTFourierTransformationOfFermionVector);
  setParameterOnAllThreads(0, (void*) input);
  setParameterOnAllThreads(1, (void*) output);
  
  int blockSize = 8;
  if (OperationMode == 1) blockSize = 4;
  if (OperationMode == 2) blockSize = 2;  
  
  int intParas[9];
  intParas[0] = L0;
  intParas[1] = L1;
  intParas[2] = L2;
  intParas[3] = L3;
  intParas[4] = true;
  intParas[5] = blockSize;    
  intParas[6] = 1;       // intertwined?  
  intParas[7] = threadCount;       
  intParas[8] = tuneLevel;       
  
  setParameterOnAllThreads(2, ((void*) (&(intParas[0]))));
  setParameterOnAllThreads(3, ((void*) fftPlanDescriptor));

  threadedExecute(1);
  fftPlanDescriptor = (char*) Threads[0]->getExecutionParameter(3);
  
  bool succ = setFFTPlan_xFFTFourierTransformationOfFermionVector(threadCount, L0, L1, L2, L3, fftPlanDescriptor);
  if (!succ) {
    printf("Error in MultiThreadedOperations::tune_xFFTFourierTransformationOfFermionVector: Could not set Plan on all threads\n");
    exit(0);
  }
}


bool MultiThreadedOperations::setFFTPlan_xFFTFourierTransformationOfFermionVector(int threadCount, int L0, int L1, int L2, int L3, char* fftPlanDescriptor) {
  setExecutionCommandOnAllThreads(ThreadControl_ExecutionCommand_setFFTPlan_xFFTFourierTransformationOfFermionVector);
  int blockSize = 8;
  if (OperationMode == 1) blockSize = 4;
  if (OperationMode == 2) blockSize = 2;  
  
  int intParas[8];
  intParas[0] = L0;
  intParas[1] = L1;
  intParas[2] = L2;
  intParas[3] = L3;
  intParas[4] = true;
  intParas[5] = blockSize;    
  intParas[6] = 1;       // intertwined?  
  intParas[7] = threadCount;       
  setParameterOnAllThreads(0, ((void*) (&(intParas[0]))));
  setParameterOnAllThreads(1, ((void*) fftPlanDescriptor));
 
  threadedExecute(threadCompleteMask);
  
  bool res = true;
  for (int I=0; I<ThreadCount; I++) {
    bool b = (bool) Threads[I]->getExecutionParameter(2);
    res = b & res;
  }
  return res; 
}
