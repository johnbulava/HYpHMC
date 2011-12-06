#include "ThreadControl.h"
#include "SSEroutines.h"

ThreadControl::ThreadControl(int threadid, int nodeid, int cpuAllocMode) {
  if (LogLevel>2) printf("Initialize ThreadControl with threadID=%d, nodeID=%d, and CPUAllocMode=%d\n",threadid, nodeid, cpuAllocMode);
  ThreadID = threadid;
  NodeID = nodeid;
  CPUAllocationMode = cpuAllocMode;
  
  xSSEObj = new xSSE();
  xFFT = NULL;
  setExecutionFlag(0, xSSEObj);    
  setTerminationFlag(0, xSSEObj);      

  ExecutionParameters = new void*[ThreadControlParameterMAX];
  for (int I=0; I<ThreadControlParameterMAX; I++) {
    ExecutionParameters[I] = NULL;
  }
  ExecutionCommand = 0;
  workingSpace = createSuperAlignedComplex(ThreadControl_WorkingSpaceSize);
  outputSpace = createSuperAlignedComplex(ThreadControl_OutputSpaceSize);  
  for (int I=0; I<ThreadControl_WorkingSpaceSize; I++) {
    workingSpace[I].x = 0;
    workingSpace[I].y = 0;    
  }
  for (int I=0; I<ThreadControl_OutputSpaceSize; I++) {
    outputSpace[I].x = 0;
    outputSpace[I].y = 0;    
  }
  
  for (int I=0; I<ThreadControl_FFTPlanDBMax; I++) {
    FFTPlanDB[I].fftwPlan = NULL;
    resetFFTPlanEntry(I);
  }
  FFTPlanDBTotalEntriesCount = 0;
  FFTPlanDBIndexOfNextEntry = 0;
  pThreadMutex = NULL;
}


ThreadControl::~ThreadControl() {	
  for (int I=0; I<ThreadControl_FFTPlanDBMax; I++) {
    resetFFTPlanEntry(I);
  }

  delete xSSEObj;
  destroySuperAlignedComplex(workingSpace);
  destroySuperAlignedComplex(outputSpace);
  delete[] ExecutionParameters;
  delete xFFT;
}


void ThreadControl::setThreadObject(pthread_t threadObj, pthread_mutex_t* mutObj) {
  pThreadObject = threadObj;
  pThreadMutex = mutObj;
}


pthread_t ThreadControl::getThreadObject() {
  return pThreadObject;
}


int ThreadControl::getThreadID() {
  return ThreadID;
}


int ThreadControl::getNodeID() {
  return NodeID;
}


int ThreadControl::getCPUAllocationMode() {
  return CPUAllocationMode;
}


long int ThreadControl::getExecutionFlag(xSSE* xSSEObjLoc) {
  long int* ExecutionFlagAddr = &(ExecutionFlag);
  long int res = 0;
  xSSEObjLoc->xSSE_ReadLongIntFromMemAddr_Wrapper(ExecutionFlagAddr, res);
  return res; 
}


long int ThreadControl::getTerminationFlag(xSSE* xSSEObjLoc) {
  long int* TerminationFlagAddr = &(TerminationFlag);
  long int res = 0;
  xSSEObjLoc->xSSE_ReadLongIntFromMemAddr_Wrapper(TerminationFlagAddr, res);
  return res;
}

void ThreadControl::setExecutionFlag(long int flag, xSSE* xSSEObjLoc) {
  long int* ExecutionFlagAddr = &(ExecutionFlag);
  xSSEObjLoc->xSSE_WriteLongIntToMemAddr_Wrapper(ExecutionFlagAddr, flag);
}


void ThreadControl::setTerminationFlag(long int flag, xSSE* xSSEObjLoc) {
  long int* TerminationFlagAddr = &(TerminationFlag);
  xSSEObjLoc->xSSE_WriteLongIntToMemAddr_Wrapper(TerminationFlagAddr, flag);
}


void ThreadControl::setExecutionParameter(int nr, void* para) {
  if ((nr<0) || (nr>=ThreadControlParameterMAX)) {
    printf("ERROR in ThreadControl::setExecutionParameter: nr = %d\n", nr);
    exit(0);
  }
  ExecutionParameters[nr] = para;
}


void* ThreadControl::getExecutionParameter(int nr) {
  if ((nr<0) || (nr>=ThreadControlParameterMAX)) {
    printf("ERROR in ThreadControl::getExecutionParameter: nr = %d\n", nr);
    exit(0);
  }
  return ExecutionParameters[nr];
}


void ThreadControl::setExecutionCommand(int cmd) {
  ExecutionCommand = cmd;
}


int ThreadControl::getExecutionCommand() {
  return ExecutionCommand;
}
  
  
void ThreadControl::allocateDistributedMemoryObject() {
  DistributedMemoryObject* memObj = (DistributedMemoryObject*) ExecutionParameters[0];
  int size = *((int*) ExecutionParameters[1]);
  if (LogLevel>3) printf("allocateDistributedMemoryObject called on Thread %d on Node %d with requested size %d\n",ThreadID,NodeID,size);

  pthread_mutex_lock (pThreadMutex);
  
  numa_set_localalloc();
  numa_set_strict(1);
  memObj->dataSets[ThreadID] = createSuperAlignedComplex(size);
  
  for (int I=0; I<size; I++) {
     memObj->dataSets[ThreadID][I].x = 0;
     memObj->dataSets[ThreadID][I].y = 0;
  }

  pthread_mutex_unlock (pThreadMutex);
}


void ThreadControl::scalarProductOfFermionVectors() {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(ThreadID);
  DistributedMemoryObject* memObj1 = (DistributedMemoryObject*) ExecutionParameters[0];
  DistributedMemoryObject* memObj2 = (DistributedMemoryObject*) ExecutionParameters[1];
  Complex* v0 = memObj1->dataSets[ThreadID];
  Complex* v1 = memObj2->dataSets[ThreadID];
  
  int* intParas = (int*) ExecutionParameters[2];;
  int L0 = intParas[0];
  int L1 = intParas[1];
  int L2 = intParas[2];
  int L3 = intParas[3];
  int blockSize = intParas[4];
    
  if (LogLevel>4) printf("scalarProductOfFermionVectors called on thread %d on node %d with L0=%d, L1=%d, L2=%d, L3=%d, and blocksize = %d\n", ThreadID, NodeID, L0, L1, L2, L3, blockSize);
  
  int i0,i1,i2,i3;

  int count = 0;
  int xtrAdd1 = xtraSize1*blockSize*(L3+xtraSize3)*(L2+xtraSize2);
  int xtrAdd2 = xtraSize2*blockSize*(L3+xtraSize3);
  int xtrAdd3 = xtraSize3*blockSize;

  Complex res(0,0);

  for (i0=0; i0<L0; i0++) {
    for (i1=0; i1<L1; i1++) {
      for (i2=0; i2<L2; i2++) {
        for (i3=0; i3<L3; i3++) {

          for (int I2=0; I2<blockSize; I2++) {
            res.x = res.x + (v0[count+I2].x*v1[count+I2].x + v0[count+I2].y*v1[count+I2].y);
            res.y = res.y + (v0[count+I2].x*v1[count+I2].y - v0[count+I2].y*v1[count+I2].x);
          } 

          count += blockSize;
	}
	count += xtrAdd3;
      }
      count += xtrAdd2;
    }
    count += xtrAdd1;
  }

  outputSpace[0] = res;
  ExecutionParameters[3] = (void*) &(outputSpace[0]);
  addPerformanceProfilingItem("ThreadControl::scalarProductOfFermionVectors", performanceProfilerStartCycle, ThreadID);
}


void ThreadControl::vectorNormOfFermionVector() {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(ThreadID);
  DistributedMemoryObject* memObj = (DistributedMemoryObject*) ExecutionParameters[0];
  Complex* v0 = memObj->dataSets[ThreadID];
  
  int* intParas = (int*) ExecutionParameters[2];;
  int L0 = intParas[0];
  int L1 = intParas[1];
  int L2 = intParas[2];
  int L3 = intParas[3];
  int blockSize = intParas[4];
  
  if (LogLevel>4) printf("vectorNormOfFermionVector called on thread %d on node %d with L0=%d, L1=%d, L2=%d, L3=%d, and blocksize = %d\n", ThreadID, NodeID, L0, L1, L2, L3, blockSize);
  
  int i0,i1,i2,i3;

  int count = 0;
  int xtrAdd1 = xtraSize1*blockSize*(L3+xtraSize3)*(L2+xtraSize2);
  int xtrAdd2 = xtraSize2*blockSize*(L3+xtraSize3);
  int xtrAdd3 = xtraSize3*blockSize;

  Complex res(0,0);

  for (i0=0; i0<L0; i0++) {
    for (i1=0; i1<L1; i1++) {
      for (i2=0; i2<L2; i2++) {
        for (i3=0; i3<L3; i3++) {

          for (int I2=0; I2<blockSize; I2++) {
            res.x = res.x + (v0[count+I2].x*v0[count+I2].x + v0[count+I2].y*v0[count+I2].y);
          } 

          count += blockSize;
	}
	count += xtrAdd3;
      }
      count += xtrAdd2;
    }
    count += xtrAdd1;
  }

  outputSpace[0] = res;
  ExecutionParameters[3] = (void*) &(outputSpace[0]);
  addPerformanceProfilingItem("ThreadControl::vectorNormOfFermionVector", performanceProfilerStartCycle, ThreadID);
}
  

void ThreadControl::multiplyFermionVectorWithComplexNumber() {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(ThreadID);
  DistributedMemoryObject* x = (DistributedMemoryObject*) ExecutionParameters[0];
  DistributedMemoryObject* y = (DistributedMemoryObject*) ExecutionParameters[1];
  Complex* v0 = x->dataSets[ThreadID];
  Complex* v1 = y->dataSets[ThreadID];
  
  int* intParas = (int*) ExecutionParameters[2];;
  int L0 = intParas[0];
  int L1 = intParas[1];
  int L2 = intParas[2];
  int L3 = intParas[3];
  int blockSize = intParas[4];
  
  Complex alpha = *((Complex*) ExecutionParameters[3]);
    
  if (LogLevel>4) printf("multiplyFermionVectorWithComplexNumber called on thread %d on node %d with alpha=(%1.3e,%1.3e), L0=%d, L1=%d, L2=%d, L3=%d, and blocksize = %d\n", ThreadID, NodeID, alpha.x, alpha.y, L0, L1, L2, L3, blockSize);
  
  int i0,i1,i2,i3;

  int count = 0;
  int xtrAdd1 = xtraSize1*blockSize*(L3+xtraSize3)*(L2+xtraSize2);
  int xtrAdd2 = xtraSize2*blockSize*(L3+xtraSize3);
  int xtrAdd3 = xtraSize3*blockSize;

  for (i0=0; i0<L0; i0++) {
    for (i1=0; i1<L1; i1++) {
      for (i2=0; i2<L2; i2++) {
        for (i3=0; i3<L3; i3++) {

          for (int I2=0; I2<blockSize; I2++) {
            double dummy =   (alpha.x*v0[count+I2].x - alpha.y*v0[count+I2].y);
            v1[count+I2].y = (alpha.y*v0[count+I2].x + alpha.x*v0[count+I2].y);
            v1[count+I2].x = dummy;
          } 

          count += blockSize;
	}
	count += xtrAdd3;
      }
      count += xtrAdd2;
    }
    count += xtrAdd1;
  }
  addPerformanceProfilingItem("ThreadControl::multiplyFermionVectorWithComplexNumber", performanceProfilerStartCycle, ThreadID);
}


void ThreadControl::multiplyFermionVectorWithTimeIndexedComplexScalars() {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(ThreadID);
  DistributedMemoryObject* x = (DistributedMemoryObject*) ExecutionParameters[0];
  DistributedMemoryObject* y = (DistributedMemoryObject*) ExecutionParameters[1];
  Complex* v0 = x->dataSets[ThreadID];
  Complex* v1 = y->dataSets[ThreadID];
  
  int* intParas = (int*) ExecutionParameters[2];;
  int L0 = intParas[0];
  int L1 = intParas[1];
  int L2 = intParas[2];
  int L3 = intParas[3];
  int blockSize = intParas[4];
  
  DistributedMemoryObject* scalars = (DistributedMemoryObject*) ExecutionParameters[3];
  Complex* vs = scalars->dataSets[ThreadID];

    
  if (LogLevel>4) printf("multiplyFermionVectorWithTimeIndexedComplexScalars called on thread %d on node %d with scalars[0]=(%1.3e,%1.3e), L0=%d, L1=%d, L2=%d, L3=%d, and blocksize = %d\n", ThreadID, NodeID, vs[0].x, vs[0].y, L0, L1, L2, L3, blockSize);
  
  int count = 0;
  int xtrAdd1 = xtraSize1*blockSize*(L3+xtraSize3)*(L2+xtraSize2);
  int xtrAdd2 = xtraSize2*blockSize*(L3+xtraSize3);
  int xtrAdd3 = xtraSize3*blockSize;
  int ind[4];

  int timeDirection = 0;
  int oneDimSizeLargest = L0;
  if (L1 > oneDimSizeLargest) { 
    timeDirection = 1;
    oneDimSizeLargest = L1;
  }
  if (L2 > oneDimSizeLargest) {
    timeDirection = 2;
    oneDimSizeLargest = L2;
  }
  if (L3 > oneDimSizeLargest) {
    timeDirection = 3;
    oneDimSizeLargest = L3;
  }

  for (ind[0]=0; ind[0]<L0; ind[0]++) {
    for (ind[1]=0; ind[1]<L1; ind[1]++) {
      for (ind[2]=0; ind[2]<L2; ind[2]++) {
        for (ind[3]=0; ind[3]<L3; ind[3]++) {

          double alphax = vs[ind[timeDirection]].x;
          double alphay = vs[ind[timeDirection]].y;
          for (int I2=0; I2<blockSize; I2++) {
            double dummy =   (alphax*v0[count+I2].x - alphay*v0[count+I2].y);
            v1[count+I2].y = (alphay*v0[count+I2].x + alphax*v0[count+I2].y);
            v1[count+I2].x = dummy;
          } 

          count += blockSize;
	}
	count += xtrAdd3;
      }
      count += xtrAdd2;
    }
    count += xtrAdd1;
  }
  addPerformanceProfilingItem("ThreadControl::multiplyFermionVectorWithTimeIndexedComplexScalars", performanceProfilerStartCycle, ThreadID);
}


void ThreadControl::zeroFermionVector() {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(ThreadID);
  DistributedMemoryObject* x = (DistributedMemoryObject*) ExecutionParameters[0];
  Complex* v0 = x->dataSets[ThreadID];
  
  int* intParas = (int*) ExecutionParameters[2];
  int L0 = intParas[0];
  int L1 = intParas[1];
  int L2 = intParas[2];
  int L3 = intParas[3];
  int blockSize = intParas[4];
     
  if (LogLevel>4) printf("zeroFermionVector called on thread %d on node %d with L0=%d, L1=%d, L2=%d, L3=%d, and blocksize = %d\n", ThreadID, NodeID, L0, L1, L2, L3, blockSize);
  
  int i0,i1,i2,i3;

  int count = 0;
  int xtrAdd1 = xtraSize1*blockSize*(L3+xtraSize3)*(L2+xtraSize2);
  int xtrAdd2 = xtraSize2*blockSize*(L3+xtraSize3);
  int xtrAdd3 = xtraSize3*blockSize;

  for (i0=0; i0<L0; i0++) {
    for (i1=0; i1<L1; i1++) {
      for (i2=0; i2<L2; i2++) {
        for (i3=0; i3<L3; i3++) {

          for (int I2=0; I2<blockSize; I2++) {
            v0[count+I2].x = 0;
            v0[count+I2].y = 0;
          } 

          count += blockSize;
	}
	count += xtrAdd3;
      }
      count += xtrAdd2;
    }
    count += xtrAdd1;
  }
  addPerformanceProfilingItem("ThreadControl::zeroFermionVector", performanceProfilerStartCycle, ThreadID);
}


void ThreadControl::zeroUniformlyDistributedComplexVector() {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(ThreadID);
  DistributedMemoryObject* x = (DistributedMemoryObject*) ExecutionParameters[0];
  Complex* v0 = x->dataSets[ThreadID];
  
  int size = *((int*) ExecutionParameters[2]);
    
  if (LogLevel>4) printf("zeroUniformlyDistributedComplexVector called on thread %d on node %d with size = %d\n", ThreadID, NodeID, size);
  
  for (int I=0; I<size; I++) {
    v0[I].x = 0;
    v0[I].y = 0;    
  }
  addPerformanceProfilingItem("ThreadControl::zeroUniformlyDistributedComplexVector", performanceProfilerStartCycle, ThreadID);
}


void ThreadControl::vectorAdditionOfUniformlyDistributedComplexVectors() {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(ThreadID);
  DistributedMemoryObject* x = (DistributedMemoryObject*) ExecutionParameters[0];
  DistributedMemoryObject* y = (DistributedMemoryObject*) ExecutionParameters[1];
  Complex* v0 = x->dataSets[ThreadID];
  Complex* v1 = y->dataSets[ThreadID];
  
  int size = *((int*) ExecutionParameters[2]);
  
  Complex alpha = *((Complex*) ExecutionParameters[3]);
    
  if (LogLevel>4) printf("vectorAdditionOfUniformlyDistributedComplexVectors called on thread %d on node %d with alpha=(%1.3e,%1.3e), size = %d\n", ThreadID, NodeID, alpha.x, alpha.y, size);
  
  for (int I=0; I<size; I++) {
    v1[I].x = v1[I].x + (alpha.x*v0[I].x - alpha.y*v0[I].y);
    v1[I].y = v1[I].y + (alpha.y*v0[I].x + alpha.x*v0[I].y);
  } 
  addPerformanceProfilingItem("ThreadControl::vectorAdditionOfUniformlyDistributedComplexVectors", performanceProfilerStartCycle, ThreadID);
}


void ThreadControl::vectorAdditionOfFermionVectors() {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(ThreadID);
  DistributedMemoryObject* x = (DistributedMemoryObject*) ExecutionParameters[0];
  DistributedMemoryObject* y = (DistributedMemoryObject*) ExecutionParameters[1];
  Complex* v0 = x->dataSets[ThreadID];
  Complex* v1 = y->dataSets[ThreadID];
  
  int* intParas = (int*) ExecutionParameters[2];;
  int L0 = intParas[0];
  int L1 = intParas[1];
  int L2 = intParas[2];
  int L3 = intParas[3];
  int blockSize = intParas[4];
  
  Complex alpha = *((Complex*) ExecutionParameters[3]);
    
  if (LogLevel>4) printf("vectorProductOfFermionVectors called on thread %d on node %d with alpha=(%1.3e,%1.3e), L0=%d, L1=%d, L2=%d, L3=%d, and blocksize = %d\n", ThreadID, NodeID, alpha.x, alpha.y, L0, L1, L2, L3, blockSize);
  
  int i0,i1,i2,i3;

  int count = 0;
  int xtrAdd1 = xtraSize1*blockSize*(L3+xtraSize3)*(L2+xtraSize2);
  int xtrAdd2 = xtraSize2*blockSize*(L3+xtraSize3);
  int xtrAdd3 = xtraSize3*blockSize;

  for (i0=0; i0<L0; i0++) {
    for (i1=0; i1<L1; i1++) {
      for (i2=0; i2<L2; i2++) {
        for (i3=0; i3<L3; i3++) {

          for (int I2=0; I2<blockSize; I2++) {
            v1[count+I2].x = v1[count+I2].x + (alpha.x*v0[count+I2].x - alpha.y*v0[count+I2].y);
            v1[count+I2].y = v1[count+I2].y + (alpha.y*v0[count+I2].x + alpha.x*v0[count+I2].y);
          } 

          count += blockSize;
	}
	count += xtrAdd3;
      }
      count += xtrAdd2;
    }
    count += xtrAdd1;
  }
  addPerformanceProfilingItem("ThreadControl::vectorAdditionOfFermionVectors", performanceProfilerStartCycle, ThreadID);
}


void ThreadControl::resetFFTPlanEntry(int nr) {
  FFTPlanDB[nr].input = NULL;
  FFTPlanDB[nr].output = NULL;
  FFTPlanDB[nr].L0 = 0;
  FFTPlanDB[nr].L1 = 0;
  FFTPlanDB[nr].L2 = 0;
  FFTPlanDB[nr].L3 = 0;
  FFTPlanDB[nr].blockSize = 0;
  FFTPlanDB[nr].Forward = false;
  FFTPlanDB[nr].intertwined = false;
  if (FFTPlanDB[nr].fftwPlan != NULL) fftw_destroy_plan(FFTPlanDB[nr].fftwPlan);
  FFTPlanDB[nr].fftwPlan = NULL;
}



fftw_plan ThreadControl::getFFTPlan(Complex* input, Complex* output, int L0, int L1, int L2, int L3, int blockSize, bool forw, bool intertwined) {
  int I;
  for (I=0; I<FFTPlanDBTotalEntriesCount; I++) {
    if ((input==FFTPlanDB[I].input) && (output==FFTPlanDB[I].output) && (forw==FFTPlanDB[I].Forward)) {
      if ((L0==FFTPlanDB[I].L0) && (L1==FFTPlanDB[I].L1) && (L2==FFTPlanDB[I].L2) && (L3==FFTPlanDB[I].L3) && (blockSize==FFTPlanDB[I].blockSize) && (intertwined==FFTPlanDB[I].intertwined)) {
        return FFTPlanDB[I].fftwPlan;
      }
    }
  }
  
  if (!intertwined) {
    printf("ERROR in ThreadControl::getFFTPlan: Only intertwined mode supported!!!\n ");
    exit(0);
  }
  
  
  int* n = new int[4];
  n[0] = L0;
  n[1] = L1;
  n[2] = L2;
  n[3] = L3;
  
  int rank = 4;
  int howmany = blockSize;
  int* inembed = new int[4];
  inembed[0] = L0;
  inembed[1] = L1+xtraSize1;
  inembed[2] = L2+xtraSize2;
  inembed[3] = L3+xtraSize3;  
  int istride = blockSize;
  int idist = 1;
  
  int* onembed = new int[4];
  onembed[0] = L0;
  onembed[1] = L1+xtraSize1;
  onembed[2] = L2+xtraSize2;
  onembed[3] = L3+xtraSize3;    
  int ostride = blockSize;
  int odist = 1;


  resetFFTPlanEntry(FFTPlanDBIndexOfNextEntry);
  
  FFTPlanDB[FFTPlanDBIndexOfNextEntry].input = input;
  FFTPlanDB[FFTPlanDBIndexOfNextEntry].output = output;
  FFTPlanDB[FFTPlanDBIndexOfNextEntry].L0 = L0;
  FFTPlanDB[FFTPlanDBIndexOfNextEntry].L1 = L1;
  FFTPlanDB[FFTPlanDBIndexOfNextEntry].L2 = L2;
  FFTPlanDB[FFTPlanDBIndexOfNextEntry].L3 = L3;
  FFTPlanDB[FFTPlanDBIndexOfNextEntry].blockSize = blockSize;
  FFTPlanDB[FFTPlanDBIndexOfNextEntry].Forward = forw;
  FFTPlanDB[FFTPlanDBIndexOfNextEntry].intertwined = intertwined;
  FFTPlanDB[FFTPlanDBIndexOfNextEntry].fftwPlan = NULL;


  int vecSize = blockSize*(L0)*(L1+xtraSize1)*(L2+xtraSize2)*(L3+xtraSize3);
  Complex* Save1 = createSuperAlignedComplex(vecSize);
  Complex* Save2 = createSuperAlignedComplex(vecSize);
  SSE_ZCopy(vecSize, input, 1, Save1, 1);
  SSE_ZCopy(vecSize, output, 1, Save2, 1);

  int measureFlag = FFTW_MEASURE | FFTW_EXHAUSTIVE;
  if (DebugMode) measureFlag = FFTW_MEASURE;
  fftw_plan fftwPlan = NULL;

  //Lock Thread during FFTW-Plan-Generation due to FFTW being non-Thread-Proof
  pthread_mutex_lock (pThreadMutex);
  if (LogLevel>2) printf("Generating FFT-plan in Thread %d on Node %d with Forward = %d\n", ThreadID, NodeID, forw);    
  if (forw == ExtremeFFT4D_Forward) {
    fftwPlan = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*)input, inembed,
                                  istride, idist,
	   			  (fftw_complex*)output, onembed, ostride, odist,
	  			  FFTW_FORWARD, measureFlag);
  } else {
    fftwPlan = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*)input, inembed,
                                  istride, idist,
	   			  (fftw_complex*)output, onembed, ostride, odist,
	  			  FFTW_BACKWARD, measureFlag);
  }
  if (LogLevel>2) printf("FFT-plan ready in Thread %d on Node %d\n", ThreadID, NodeID);      
  pthread_mutex_unlock (pThreadMutex);
  //Thread unlocked
  FFTPlanDB[FFTPlanDBIndexOfNextEntry].fftwPlan = fftwPlan;
  FFTPlanDBTotalEntriesCount++;
  FFTPlanDBIndexOfNextEntry++;
  if (FFTPlanDBTotalEntriesCount>=ThreadControl_FFTPlanDBMax) {
    FFTPlanDBTotalEntriesCount = ThreadControl_FFTPlanDBMax;  
  }
  if (FFTPlanDBIndexOfNextEntry>=ThreadControl_FFTPlanDBMax) {
    FFTPlanDBIndexOfNextEntry = 0;  
  }
    
  SSE_ZCopy(vecSize, Save1, 1, input, 1);
  SSE_ZCopy(vecSize, Save2, 1, output, 1);
  destroySuperAlignedComplex(Save1);
  destroySuperAlignedComplex(Save2);
  delete[] n;
  delete[] inembed;
  delete[] onembed;
  return fftwPlan;
}


void ThreadControl::perform_FFTWFourierTransformationOfFermionVector() {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(ThreadID);
  DistributedMemoryObject* input = (DistributedMemoryObject*) ExecutionParameters[0];
  DistributedMemoryObject* output = (DistributedMemoryObject*) ExecutionParameters[1];
  Complex* v0 = input->dataSets[ThreadID];
  Complex* v1 = output->dataSets[ThreadID];
  
  int* intParas = (int*) ExecutionParameters[2];;
  int L0 = intParas[0];
  int L1 = intParas[1];
  int L2 = intParas[2];
  int L3 = intParas[3];
  bool forw = intParas[4];
  int blockSize = intParas[5];
  bool intertwined = (bool)intParas[6];
      
  if (LogLevel>4) printf("perform_FFTWFourierTransformationOfFermionVector called on thread %d on node %d with forward=%d, L0=%d, L1=%d, L2=%d, L3=%d, blocksize=%d, and intertwined=%d\n", ThreadID, NodeID, forw, L0, L1, L2, L3, blockSize, intertwined);
  

  fftw_plan plan = getFFTPlan(v0, v1, L0, L1, L2, L3, blockSize, forw, intertwined);
  addPerformanceProfilingItem("ThreadControl::perform_FFTWFourierTransformationOfFermionVector-tune", performanceProfilerStartCycle, ThreadID);
  performanceProfilerStartCycle = getPerformanceProfilingStartCycle(ThreadID);
  
  fftw_execute(plan);   
  plan = NULL;  
  addPerformanceProfilingItem("ThreadControl::perform_FFTWFourierTransformationOfFermionVector-execute", performanceProfilerStartCycle, ThreadID);
}


void ThreadControl::perform_xFFTFourierTransformationOfFermionVector() {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(ThreadID);
  DistributedMemoryObject* input = (DistributedMemoryObject*) ExecutionParameters[0];
  DistributedMemoryObject* output = (DistributedMemoryObject*) ExecutionParameters[1];
  Complex* v0 = input->dataSets[ThreadID];
  Complex* v1 = output->dataSets[ThreadID];
  
  int* intParas = (int*) ExecutionParameters[2];;
  int L0 = intParas[0];
  int L1 = intParas[1];
  int L2 = intParas[2];
  int L3 = intParas[3];
  bool forw = intParas[4];
  int blockSize = intParas[5];
  bool intertwined = (bool)intParas[6];
  int threadCount = intParas[7];
      
  if (LogLevel>4) printf("perform_xFFTFourierTransformationOfFermionVector called on thread %d on node %d with forward=%d, L0=%d, L1=%d, L2=%d, L3=%d, blocksize=%d, intertwined=%d, and threadCount=%d\n", ThreadID, NodeID, forw, L0, L1, L2, L3, blockSize, intertwined, threadCount);
  
  if (xFFT!=NULL) {
    if ((xFFT->getL0()!=L0) || (xFFT->getL1()!=L1) || (xFFT->getL2()!=L2) || (xFFT->getL3()!=L3) || (xFFT->getLocalIndexCount()!=blockSize) || (xFFT->getNumberOfThreads()!=threadCount) || (xFFT->getXtraSize1()!=xtraSize1) || (xFFT->getXtraSize2()!=xtraSize2) || (xFFT->getXtraSize3()!=xtraSize3)) {
      delete xFFT;
      xFFT = NULL;    
    }
  }
  
  if (xFFT == NULL) {
    //Lock Thread during tuning
    pthread_mutex_lock (pThreadMutex);
    if (LogLevel>2) printf("Generating xFFT-plan in Thread %d on Node %d with Forward = %d\n", ThreadID, NodeID, forw);    
  
    long int* threadCoreMsk = new long int[threadCount];
    if (CPUAllocationMode==1) {
      for (int I=0; I<threadCount; I++) {
        threadCoreMsk[I] = getAffinityMaskFromNthCoreOnNodeID(NodeID, I);
      }
    } else {
      for (int I=0; I<threadCount; I++) {
        threadCoreMsk[I] = getAffinityMaskFromNodeID(NodeID);
      }
    } 
    xFFT = new ExtremeFFT4D(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, blockSize, threadCount, threadCoreMsk, true); 
    delete[] threadCoreMsk;

    xFFT->tune(ExtremeFFT4D_TuneLevel_Low);
    if (LogLevel>2) printf("xFFT-plan ready in Thread %d on Node %d\n", ThreadID, NodeID);

    pthread_mutex_unlock (pThreadMutex);
    //Thread unlocked
  }
  addPerformanceProfilingItem("ThreadControl::perform_xFFTFourierTransformationOfFermionVector-tune", performanceProfilerStartCycle, ThreadID);
  performanceProfilerStartCycle = getPerformanceProfilingStartCycle(ThreadID);
  
  xFFT->FastFourierTrafo(v0, v1, forw);
  addPerformanceProfilingItem("ThreadControl::perform_xFFTFourierTransformationOfFermionVector-execute", performanceProfilerStartCycle, ThreadID);
}


void ThreadControl::tune_xFFTFourierTransformationOfFermionVector() {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(ThreadID);
  DistributedMemoryObject* input = (DistributedMemoryObject*) ExecutionParameters[0];
  DistributedMemoryObject* output = (DistributedMemoryObject*) ExecutionParameters[1];
  Complex* v0 = input->dataSets[ThreadID];
  Complex* v1 = output->dataSets[ThreadID];
  
  int* intParas = (int*) ExecutionParameters[2];
  int L0 = intParas[0];
  int L1 = intParas[1];
  int L2 = intParas[2];
  int L3 = intParas[3];
  int blockSize = intParas[5];
  bool intertwined = (bool)intParas[6];
  int threadCount = intParas[7];
  int tuneLevel = intParas[8];
  
  if (LogLevel>2) printf("tune_xFFTFourierTransformationOfFermionVector called on thread %d on node %d with L0=%d, L1=%d, L2=%d, L3=%d, blocksize=%d, intertwined=%d, threadCount=%d, and tuneLevel=%d\n", ThreadID, NodeID, L0, L1, L2, L3, blockSize, intertwined, threadCount, tuneLevel);
  
  if (xFFT!=NULL) {
    if ((xFFT->getL0()!=L0) || (xFFT->getL1()!=L1) || (xFFT->getL2()!=L2) || (xFFT->getL3()!=L3) || (xFFT->getLocalIndexCount()!=blockSize) || (xFFT->getNumberOfThreads()!=threadCount) || (xFFT->getXtraSize1()!=xtraSize1) || (xFFT->getXtraSize2()!=xtraSize2) || (xFFT->getXtraSize3()!=xtraSize3)) {
      delete xFFT;
      xFFT = NULL;    
    }
  }
  
  //Lock Thread during tuning
  pthread_mutex_lock (pThreadMutex);
  if (xFFT == NULL) {  
    long int* threadCoreMsk = new long int[threadCount];
    if (CPUAllocationMode==1) {
      for (int I=0; I<threadCount; I++) {
        threadCoreMsk[I] = getAffinityMaskFromNthCoreOnNodeID(NodeID, I);
      }
    } else {
      for (int I=0; I<threadCount; I++) {
        threadCoreMsk[I] = getAffinityMaskFromNodeID(NodeID);
      }
    } 
    xFFT = new ExtremeFFT4D(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, blockSize, threadCount, threadCoreMsk, true); 
    delete[] threadCoreMsk;
  }

  xFFT->tune(tuneLevel, v0, v1);  
  
  char* fftPlanDescriptor = xFFT->getSelectedPlanDescriptor();
  ExecutionParameters[3] = (void*) fftPlanDescriptor;

  if (LogLevel>2) printf("xFFT-Tuning ready in Thread %d on Node %d\n", ThreadID, NodeID);

  pthread_mutex_unlock (pThreadMutex);
  //Thread unlocked  
  addPerformanceProfilingItem("ThreadControl::tune_xFFTFourierTransformationOfFermionVector", performanceProfilerStartCycle, ThreadID);
}


void ThreadControl::setFFTPlan_xFFTFourierTransformationOfFermionVector() {
  int* intParas = (int*) ExecutionParameters[0];
  int L0 = intParas[0];
  int L1 = intParas[1];
  int L2 = intParas[2];
  int L3 = intParas[3];
  int blockSize = intParas[5];
  bool intertwined = (bool)intParas[6];
  int threadCount = intParas[7];
  char* fftPlanDescriptor = (char*) ExecutionParameters[1];
      
  if (LogLevel>2) printf("setFFTPlan_xFFTFourierTransformationOfFermionVector called on thread %d on node %d with L0=%d, L1=%d, L2=%d, L3=%d, blocksize=%d, intertwined=%d, threadCount=%d, and fftPlanDescriptor=%s\n", ThreadID, NodeID, L0, L1, L2, L3, blockSize, intertwined, threadCount, fftPlanDescriptor);
  
  if (xFFT!=NULL) {
    if ((xFFT->getL0()!=L0) || (xFFT->getL1()!=L1) || (xFFT->getL2()!=L2) || (xFFT->getL3()!=L3) || (xFFT->getLocalIndexCount()!=blockSize) || (xFFT->getNumberOfThreads()!=threadCount) || (xFFT->getXtraSize1()!=xtraSize1) || (xFFT->getXtraSize2()!=xtraSize2) || (xFFT->getXtraSize3()!=xtraSize3)) {
      delete xFFT;
      xFFT = NULL;    
    }
  }
  
  //Lock Thread during tuning
  pthread_mutex_lock (pThreadMutex);
  if (xFFT == NULL) {
    long int* threadCoreMsk = new long int[threadCount];
    if (CPUAllocationMode==1) {
      for (int I=0; I<threadCount; I++) {
        threadCoreMsk[I] = getAffinityMaskFromNthCoreOnNodeID(NodeID, I);
      }
    } else {
      for (int I=0; I<threadCount; I++) {
        threadCoreMsk[I] = getAffinityMaskFromNodeID(NodeID);
      }
    } 
    xFFT = new ExtremeFFT4D(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, blockSize, threadCount, threadCoreMsk, true); 
    delete[] threadCoreMsk;    
  }

  bool b = xFFT->setSelectedPlanFromDescriptor(fftPlanDescriptor);
  ExecutionParameters[2] = (void*) b;
  
  pthread_mutex_unlock (pThreadMutex);
  //Thread unlocked  
}


void ThreadControl::copyFermionVector() {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(ThreadID);
  DistributedMemoryObject* x = (DistributedMemoryObject*) ExecutionParameters[0];
  DistributedMemoryObject* y = (DistributedMemoryObject*) ExecutionParameters[1];
  Complex* v0 = x->dataSets[ThreadID];
  Complex* v1 = y->dataSets[ThreadID];
  
  int* intParas = (int*) ExecutionParameters[2];;
  int L0 = intParas[0];
  int L1 = intParas[1];
  int L2 = intParas[2];
  int L3 = intParas[3];
  int blockSize = intParas[4];
  
  if (LogLevel>4) printf("copyFermionVector called on thread %d on node %d with L0=%d, L1=%d, L2=%d, L3=%d, and blocksize = %d\n", ThreadID, NodeID, L0, L1, L2, L3, blockSize);

  xSSEObj->xSSE_ComplexCopy_Wrappper(v0, v1, blockSize, L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3);
  addPerformanceProfilingItem("ThreadControl::copyFermionVector", performanceProfilerStartCycle, ThreadID);
}


void ThreadControl::perform_PiModeRemoverOperatorMultiplicationInFourierSpace() {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(ThreadID);
  DistributedMemoryObject* x = (DistributedMemoryObject*) ExecutionParameters[0];
  DistributedMemoryObject* y = (DistributedMemoryObject*) ExecutionParameters[1];
  DistributedMemoryObject* z = (DistributedMemoryObject*) ExecutionParameters[3];
  Complex* v0 = x->dataSets[ThreadID];
  Complex* v1 = y->dataSets[ThreadID];
  long int* Index_PiModes = (long int*) (z->dataSets[ThreadID]);
  
  int* intParas = (int*) ExecutionParameters[2];
  int L0 = intParas[0];
  int L1 = intParas[1];
  int L2 = intParas[2];
  int L3 = intParas[3];
  int blockSize = intParas[4];
  
  if (LogLevel>4) printf("perform_PiModeRemoverOperatorMultiplicationInFourierSpace called on thread %d on node %d with L0=%d, L1=%d, L2=%d, L3=%d, and blocksize = %d\n", ThreadID, NodeID, L0, L1, L2, L3, blockSize);
  
  if (v0 != v1) {
    copyFermionVector();
  }

  for (int I=0; I<15; I++) {
    int p = (blockSize*Index_PiModes[I])/8;
    if (p>=0) {
      for (int I2=0; I2<blockSize; I2++) {
        v1[p+I2].x = 0.0;
        v1[p+I2].y = 0.0;
      }
    }
  }
  addPerformanceProfilingItem("ThreadControl::perform_PiModeRemoverOperatorMultiplicationInFourierSpace", performanceProfilerStartCycle, ThreadID);
}


void ThreadControl::perform_DiracTypeOperatorMultiplicationInFourierSpace() {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(ThreadID);
  DistributedMemoryObject* x = (DistributedMemoryObject*) ExecutionParameters[0];
  DistributedMemoryObject* y = (DistributedMemoryObject*) ExecutionParameters[1];
  DistributedMemoryObject* z = (DistributedMemoryObject*) ExecutionParameters[3];
  DistributedMemoryObject* w = (DistributedMemoryObject*) ExecutionParameters[4];
  Complex* sinP = z->dataSets[ThreadID];
  Complex* auxData = w->dataSets[ThreadID];
  
  int* intParas = (int*) ExecutionParameters[2];
  int L0 = intParas[0];
  int L1 = intParas[1];
  int L2 = intParas[2];
  int L3 = intParas[3];
  int blockSize = intParas[4];
  
  if (LogLevel>4) printf("perform_DiracTypeOperatorMultiplicationInFourierSpace called on thread %d on node %d with L0=%d, L1=%d, L2=%d, L3=%d, and blocksize = %d\n", ThreadID, NodeID, L0, L1, L2, L3, blockSize);


  Complex *input0, *input1, *input2, *input3;
  Complex *output0, *output1, *output2, *output3; 
  input0 = input1 = input2 = input3 = x->dataSets[ThreadID];
  output0 = output1 = output2 = output3 = y->dataSets[ThreadID];  
  
  int innerLoop1 = 1;
  int blsHalf = blockSize/2;
  int countAdd = blockSize;
  int countStart = 0;
  int ind0, ind1, ind2, ind3;
  ind0 = ind1 = ind2 = ind3 = 0;
  
  if (blockSize==8) {
    input0 = input1 = input2 = input3 = x->dataSets[ThreadID];
    output0 = output1 = output2 = output3 = y->dataSets[ThreadID];  
    innerLoop1 = 2;
    countAdd = blsHalf;
    ind0=0;
    ind1=1;
    ind2=2;
    ind3=3;    
  }
  if (blockSize==4) {
    input0 = input1 = x->dataSets[0];
    input2 = input3 = x->dataSets[1];
    output0 = output1 = y->dataSets[0];
    output2 = output3 = y->dataSets[1];  
    ind0 = ind2 = 0;
    ind1 = ind3 = 1;
    countStart = 2*ThreadID;
  } 
  if (blockSize==2) {  
    input0 = x->dataSets[0];
    input1 = x->dataSets[1];
    input2 = x->dataSets[2];
    input3 = x->dataSets[3];
    output0 = y->dataSets[0];
    output1 = y->dataSets[1];
    output2 = y->dataSets[2];  
    output3 = y->dataSets[3];  
    ind0 = ind1 = ind2 = ind3 = 0;
    if (ThreadID!=0) countStart = 1;
  }
  

  int count = countStart;
  int countAux = 0;
  int xtrAdd1 = xtraSize1*blockSize*(L3+xtraSize3)*(L2+xtraSize2);
  int xtrAdd2 = xtraSize2*blockSize*(L3+xtraSize3);
  int xtrAdd3 = xtraSize3*blockSize;
  double p0,p1,p2,p3;

  for (int I0=2*(L0-1); I0>=0; I0-=2) {
    p0 = sinP[0*256+I0].x;
    for (int I1=2*(L1-1); I1>=0; I1-=2) {
      p1 = sinP[1*256+I1].x;
      for (int I2=2*(L2-1); I2>=0; I2-=2) {
        p2 = sinP[2*256+I2].x;
        for (int I3=2*(L3-1); I3>=0; I3-=2) {
          p3 = sinP[3*256+I3].x;
	  
	  for (int I=0; I<innerLoop1; I++) {
            workingSpace[0].x = auxData[countAux].x * input0[count+ind0].x + auxData[countAux].y*(p3*input2[count+ind2].x - p0*input2[count+ind2].y 
	                                                                   + p1*input3[count+ind3].x + p2*input3[count+ind3].y);
            workingSpace[0].y = auxData[countAux].x * input0[count+ind0].y + auxData[countAux].y*(p3*input2[count+ind2].y + p0*input2[count+ind2].x 
	                                                                   + p1*input3[count+ind3].y - p2*input3[count+ind3].x);
							       
            workingSpace[1].x = auxData[countAux].x * input1[count+ind1].x + auxData[countAux].y*(p1*input2[count+ind2].x - p2*input2[count+ind2].y 
	                                                                   - p3*input3[count+ind3].x - p0*input3[count+ind3].y);
            workingSpace[1].y = auxData[countAux].x * input1[count+ind1].y + auxData[countAux].y*(p1*input2[count+ind2].y + p2*input2[count+ind2].x 
	                                                                   - p3*input3[count+ind3].y + p0*input3[count+ind3].x);

            workingSpace[2].x = auxData[countAux].x * input2[count+ind2].x + auxData[countAux].y*(-p3*input0[count+ind0].x - p0*input0[count+ind0].y 
	                                                                   - p1*input1[count+ind1].x - p2*input1[count+ind1].y);
            workingSpace[2].y = auxData[countAux].x * input2[count+ind2].y + auxData[countAux].y*(-p3*input0[count+ind0].y + p0*input0[count+ind0].x 
	                                                                   - p1*input1[count+ind1].y + p2*input1[count+ind1].x);

            workingSpace[3].x = auxData[countAux].x * input3[count+ind3].x + auxData[countAux].y*(-p1*input0[count+ind0].x + p2*input0[count+ind0].y 
	                                                                   + p3*input1[count+ind1].x - p0*input1[count+ind1].y);
            workingSpace[3].y = auxData[countAux].x * input3[count+ind3].y + auxData[countAux].y*(-p1*input0[count+ind0].y - p2*input0[count+ind0].x 
	                                                                   + p3*input1[count+ind1].y + p0*input1[count+ind1].x);
				  
	    output0[count+ind0].x = workingSpace[0].x;
  	    output0[count+ind0].y = workingSpace[0].y;
	    output1[count+ind1].x = workingSpace[1].x;
  	    output1[count+ind1].y = workingSpace[1].y;
	    output2[count+ind2].x = workingSpace[2].x;
  	    output2[count+ind2].y = workingSpace[2].y;
	    output3[count+ind3].x = workingSpace[3].x;
  	    output3[count+ind3].y = workingSpace[3].y;
	  
	    count += countAdd;
	  }
          countAux++;
        }
	count += xtrAdd3;
      }
      count += xtrAdd2;
    }
    count += xtrAdd1;
  }
  addPerformanceProfilingItem("ThreadControl::perform_DiracTypeOperatorMultiplicationInFourierSpace", performanceProfilerStartCycle, ThreadID);
}


void ThreadControl::perform_yBOperatorMultiplication() {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(ThreadID);
  DistributedMemoryObject* m1 = (DistributedMemoryObject*) ExecutionParameters[0];
  DistributedMemoryObject* m2 = (DistributedMemoryObject*) ExecutionParameters[1];
  DistributedMemoryObject* z = (DistributedMemoryObject*) ExecutionParameters[3];
  Complex* in = m1->dataSets[ThreadID];
  Complex* out = m2->dataSets[ThreadID];
  double* phi = (double*) z->dataSets[ThreadID];
  
  int* intParas = (int*) ExecutionParameters[2];
  int L0 = intParas[0];
  int L1 = intParas[1];
  int L2 = intParas[2];
  int L3 = intParas[3];
  int blockSize = intParas[4];


  double* doubleParas = (double*) ExecutionParameters[4];
  double y = doubleParas[0];
  double split = doubleParas[1];
  
  if (LogLevel>4) printf("perform_yBOperatorMultiplication called on thread %d on node %d with L0=%d, L1=%d, L2=%d, L3=%d, y=%f, split=%f, and blocksize = %d\n", ThreadID, NodeID, L0, L1, L2, L3, y, split, blockSize);
  
  int i0,i1,i2,i3;

  int innerLoopEnd2 = 1;
  if (blockSize==8) innerLoopEnd2 = 2; 
  int innerLoopEnd1 = blockSize/(2*innerLoopEnd2);
  int blShalf = blockSize/2;
  double fp = 0.5*(1+split);
  double fm = 0.5*(1-split);
  

  double globalRelFac = 1;
  if ((blockSize==4) && (ThreadID==1)) globalRelFac = -1;
  if ((blockSize==2) && (ThreadID==2)) globalRelFac = -1;
  if ((blockSize==2) && (ThreadID==3)) globalRelFac = -1;
  

  int count = 0;
  int phiCount = 0;
  int xtrAdd1 = xtraSize1*blockSize*(L3+xtraSize3)*(L2+xtraSize2);
  int xtrAdd2 = xtraSize2*blockSize*(L3+xtraSize3);
  int xtrAdd3 = xtraSize3*blockSize;

  for (i0=0; i0<L0; i0++) {
    for (i1=0; i1<L1; i1++) {
      for (i2=0; i2<L2; i2++) {
        for (i3=0; i3<L3; i3++) {
	
          double relFac = globalRelFac;
	  int c = 0;
          for (int I=0; I<innerLoopEnd2; I++) {
            for (int i=0; i<innerLoopEnd1; i++) {
              	  
              workingSpace[c].x = y*(phi[phiCount+0]*in[count+c].x - relFac*phi[phiCount+3]*in[count+c].y
                                + fm*(phi[phiCount+2]*in[count+c+blShalf].x - phi[phiCount+1]*in[count+c+blShalf].y)
                                + fp*relFac*(phi[phiCount+2]*in[count+c+blShalf].x - phi[phiCount+1]*in[count+c+blShalf].y));
	   
              workingSpace[c].y = y*(phi[phiCount+0]*in[count+c].y + relFac*phi[phiCount+3]*in[count+c].x
                                + fm*(phi[phiCount+2]*in[count+c+blShalf].y + phi[phiCount+1]*in[count+c+blShalf].x)
                                + fp*relFac*(phi[phiCount+2]*in[count+c+blShalf].y + phi[phiCount+1]*in[count+c+blShalf].x));

              workingSpace[c+blShalf].x = y*(fm*(phi[phiCount+2]*in[count+c].x + phi[phiCount+1]*in[count+c].y)
                                        + fp*relFac*(-phi[phiCount+2]*in[count+c].x - phi[phiCount+1]*in[count+c].y)
                                        + split*(phi[phiCount+0]*in[count+c+blShalf].x + relFac*phi[phiCount+3]*in[count+c+blShalf].y));
		       
              workingSpace[c+blShalf].y = y*(fm*(phi[phiCount+2]*in[count+c].y - phi[phiCount+1]*in[count+c].x)
                                        + fp*relFac*(-phi[phiCount+2]*in[count+c].y + phi[phiCount+1]*in[count+c].x)
                                        + split*(phi[phiCount+0]*in[count+c+blShalf].y - relFac*phi[phiCount+3]*in[count+c+blShalf].x));

	      c++;  
	    }
	    relFac *= -1;	  
	  }

          for (int I=0; I<blockSize; I++) {
  	    out[count+I].x = workingSpace[I].x;
  	    out[count+I].y = workingSpace[I].y;
	  }

          count += blockSize;
          phiCount +=4;
	}
	count += xtrAdd3;
      }
      count += xtrAdd2;
    }
    count += xtrAdd1;
  }  
  addPerformanceProfilingItem("ThreadControl::perform_yBOperatorMultiplication", performanceProfilerStartCycle, ThreadID);
}


void ThreadControl::perform_yBDaggerOperatorMultiplication() {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(ThreadID);
  DistributedMemoryObject* m1 = (DistributedMemoryObject*) ExecutionParameters[0];
  DistributedMemoryObject* m2 = (DistributedMemoryObject*) ExecutionParameters[1];
  DistributedMemoryObject* z = (DistributedMemoryObject*) ExecutionParameters[3];
  Complex* in = m1->dataSets[ThreadID];
  Complex* out = m2->dataSets[ThreadID];
  double* phi = (double*) z->dataSets[ThreadID];
  
  int* intParas = (int*) ExecutionParameters[2];
  int L0 = intParas[0];
  int L1 = intParas[1];
  int L2 = intParas[2];
  int L3 = intParas[3];
  int blockSize = intParas[4];


  double* doubleParas = (double*) ExecutionParameters[4];
  double y = doubleParas[0];
  double split = doubleParas[1];
  
  if (LogLevel>4) printf("perform_yBDaggerOperatorMultiplication called on thread %d on node %d with L0=%d, L1=%d, L2=%d, L3=%d, y=%f, split=%f, and blocksize = %d\n", ThreadID, NodeID, L0, L1, L2, L3, y, split, blockSize);
  
  int i0,i1,i2,i3;

  int innerLoopEnd2 = 1;
  if (blockSize==8) innerLoopEnd2 = 2; 
  int innerLoopEnd1 = blockSize/(2*innerLoopEnd2);
  int blShalf = blockSize/2;
  double fp = 0.5*(1+split);
  double fm = 0.5*(1-split);
  

  double globalRelFac = 1;
  if ((blockSize==4) && (ThreadID==1)) globalRelFac = -1;
  if ((blockSize==2) && (ThreadID==2)) globalRelFac = -1;
  if ((blockSize==2) && (ThreadID==3)) globalRelFac = -1;
  

  int count = 0;
  int phiCount = 0;
  int xtrAdd1 = xtraSize1*blockSize*(L3+xtraSize3)*(L2+xtraSize2);
  int xtrAdd2 = xtraSize2*blockSize*(L3+xtraSize3);
  int xtrAdd3 = xtraSize3*blockSize;

  for (i0=0; i0<L0; i0++) {
    for (i1=0; i1<L1; i1++) {
      for (i2=0; i2<L2; i2++) {
        for (i3=0; i3<L3; i3++) {
	
          double relFac = globalRelFac;
	  int c = 0;
          for (int I=0; I<innerLoopEnd2; I++) {
            for (int i=0; i<innerLoopEnd1; i++) {
              	  
              workingSpace[c].x = y*(phi[phiCount+0]*in[count+c].x + relFac*phi[phiCount+3]*in[count+c].y
                                + fm*(phi[phiCount+2]*in[count+c+blShalf].x - phi[phiCount+1]*in[count+c+blShalf].y)
                                + fp*relFac*(-phi[phiCount+2]*in[count+c+blShalf].x + phi[phiCount+1]*in[count+c+blShalf].y));
	   
              workingSpace[c].y = y*(phi[phiCount+0]*in[count+c].y - relFac*phi[phiCount+3]*in[count+c].x
                                + fm*(phi[phiCount+2]*in[count+c+blShalf].y + phi[phiCount+1]*in[count+c+blShalf].x)
                                + fp*relFac*(-phi[phiCount+2]*in[count+c+blShalf].y - phi[phiCount+1]*in[count+c+blShalf].x));

              workingSpace[c+blShalf].x = y*(fm*(phi[phiCount+2]*in[count+c].x + phi[phiCount+1]*in[count+c].y)
                                        + fp*relFac*(phi[phiCount+2]*in[count+c].x + phi[phiCount+1]*in[count+c].y)
                                        + split*(phi[phiCount+0]*in[count+c+blShalf].x - relFac*phi[phiCount+3]*in[count+c+blShalf].y));
		       
              workingSpace[c+blShalf].y = y*(fm*(phi[phiCount+2]*in[count+c].y - phi[phiCount+1]*in[count+c].x)
                                        + fp*relFac*(phi[phiCount+2]*in[count+c].y - phi[phiCount+1]*in[count+c].x)
                                        + split*(phi[phiCount+0]*in[count+c+blShalf].y + relFac*phi[phiCount+3]*in[count+c+blShalf].x));

	      c++;  
	    }
	    relFac *= -1;	  
	  }

          for (int I=0; I<blockSize; I++) {
  	    out[count+I].x = workingSpace[I].x;
  	    out[count+I].y = workingSpace[I].y;
	  }

          count += blockSize;
          phiCount +=4;
	}
	count += xtrAdd3;
      }
      count += xtrAdd2;
    }
    count += xtrAdd1;    
  }  
  addPerformanceProfilingItem("ThreadControl::perform_yBDaggerOperatorMultiplication", performanceProfilerStartCycle, ThreadID);
}


void ThreadControl::perform_RPreconditionerTypeMatrixMultiplicationInFourierSpace() {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(ThreadID);
  DistributedMemoryObject* x = (DistributedMemoryObject*) ExecutionParameters[0];
  DistributedMemoryObject* y = (DistributedMemoryObject*) ExecutionParameters[1];
  DistributedMemoryObject* z = (DistributedMemoryObject*) ExecutionParameters[3];
  Complex* v0 = x->dataSets[ThreadID];
  Complex* v1 = y->dataSets[ThreadID];
  Complex* mulData = z->dataSets[ThreadID];
  
  int* intParas = (int*) ExecutionParameters[2];
  int L0 = intParas[0];
  int L1 = intParas[1];
  int L2 = intParas[2];
  int L3 = intParas[3];
  int blockSize = intParas[4];
  
  if (LogLevel>4) printf("perform_RPreconditionerTypeMatrixMultiplicationInFourierSpace called on thread %d on node %d with L0=%d, L1=%d, L2=%d, L3=%d, and blocksize = %d\n", ThreadID, NodeID, L0, L1, L2, L3, blockSize);
  
  int i0,i1,i2,i3;

  int count = 0;
  int xtrAdd1 = xtraSize1*blockSize*(L3+xtraSize3)*(L2+xtraSize2);
  int xtrAdd2 = xtraSize2*blockSize*(L3+xtraSize3);
  int xtrAdd3 = xtraSize3*blockSize;
  int internalLoopCount = blockSize/2;
  int mulCount = 0;

  for (i0=0; i0<L0; i0++) {
    for (i1=0; i1<L1; i1++) {
      for (i2=0; i2<L2; i2++) {
        for (i3=0; i3<L3; i3++) {

          for (int i=0; i<internalLoopCount; i++) {
	    //Imagin�rteil von mulData enth�t Faktor fr bottom-quark.
            v1[count].x = v0[count].x * mulData[mulCount].x;
            v1[count].y = v0[count].y * mulData[mulCount].x;
            count++;
          }
          for (int i=0; i<internalLoopCount; i++) {
	    //Imagin�rteil von mulData enth�t Faktor fr bottom-quark.
            v1[count].x = v0[count].x * mulData[mulCount].y;
            v1[count].y = v0[count].y * mulData[mulCount].y;
            count++;
          }

          mulCount++;
	}
	count += xtrAdd3;
      }
      count += xtrAdd2;
    }
    count += xtrAdd1;
  }  
  addPerformanceProfilingItem("ThreadControl::perform_RPreconditionerTypeMatrixMultiplicationInFourierSpace", performanceProfilerStartCycle, ThreadID);
}


void ThreadControl::perform_MultiplicationVectorWithDerivativesOfB() {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(ThreadID);
  DistributedMemoryObject* m1 = (DistributedMemoryObject*) ExecutionParameters[0];
  DistributedMemoryObject* m2 = (DistributedMemoryObject*) ExecutionParameters[1];
  DistributedMemoryObject* z = (DistributedMemoryObject*) ExecutionParameters[3];
  Complex* left = m1->dataSets[ThreadID];
  Complex* right = m2->dataSets[ThreadID];
  Complex* output = z->dataSets[ThreadID];
  
  int* intParas = (int*) ExecutionParameters[2];
  int L0 = intParas[0];
  int L1 = intParas[1];
  int L2 = intParas[2];
  int L3 = intParas[3];
  int blockSize = intParas[4];


  double* doubleParas = (double*) ExecutionParameters[4];
  double split = doubleParas[0];
  
  if (LogLevel>4) printf("perform_MultiplicationVectorWithDerivativesOfB called on thread %d on node %d with L0=%d, L1=%d, L2=%d, L3=%d, split=%f, and blocksize = %d\n", ThreadID, NodeID, L0, L1, L2, L3, split, blockSize);
  
  int i0,i1,i2,i3;

  int innerLoopEnd2 = 1;
  if (blockSize==8) innerLoopEnd2 = 2; 
  int innerLoopEnd1 = blockSize/(2*innerLoopEnd2);
  int blShalf = blockSize/2;
  double fp = 0.5*(1+split);
  double fm = 0.5*(1-split);
  

  double globalRelFac = 1;
  if ((blockSize==4) && (ThreadID==1)) globalRelFac = -1;
  if ((blockSize==2) && (ThreadID==2)) globalRelFac = -1;
  if ((blockSize==2) && (ThreadID==3)) globalRelFac = -1;
  

  int count = 0;
  int outputCount = 0;
  int xtrAdd1 = xtraSize1*blockSize*(L3+xtraSize3)*(L2+xtraSize2);
  int xtrAdd2 = xtraSize2*blockSize*(L3+xtraSize3);
  int xtrAdd3 = xtraSize3*blockSize;

  Complex s0,s1,s2,s3;

  for (i0=0; i0<L0; i0++) {
    for (i1=0; i1<L1; i1++) {
      for (i2=0; i2<L2; i2++) {
        for (i3=0; i3<L3; i3++) {
          output[outputCount+0].x = 0;
	  output[outputCount+0].y = 0;
	  output[outputCount+1].x = 0;
	  output[outputCount+1].y = 0;
	  output[outputCount+2].x = 0;
	  output[outputCount+2].y = 0;
	  output[outputCount+3].x = 0;
	  output[outputCount+3].y = 0;
	
          double relFac = globalRelFac;
	  int c = 0;
          for (int I=0; I<innerLoopEnd2; I++) {
            for (int i=0; i<innerLoopEnd1; i++) {
              	  
              s0.x = left[count+c].x*right[count+c].x + left[count+c].y*right[count+c].y;
              s0.y = left[count+c].x*right[count+c].y - left[count+c].y*right[count+c].x;
              s1.x = left[count+c].x*right[count+c+blShalf].x + left[count+c].y*right[count+c+blShalf].y;
              s1.y = left[count+c].x*right[count+c+blShalf].y - left[count+c].y*right[count+c+blShalf].x;
              s2.x = left[count+c+blShalf].x*right[count+c].x + left[count+c+blShalf].y*right[count+c].y;
              s2.y = left[count+c+blShalf].x*right[count+c].y - left[count+c+blShalf].y*right[count+c].x;
              s3.x = left[count+c+blShalf].x*right[count+c+blShalf].x + left[count+c+blShalf].y*right[count+c+blShalf].y;
              s3.y = left[count+c+blShalf].x*right[count+c+blShalf].y - left[count+c+blShalf].y*right[count+c+blShalf].x;	      
	      
	      output[outputCount+0].x +=  s0.x + split*s3.x;
	      output[outputCount+0].y +=  s0.y + split*s3.y;
	      output[outputCount+1].x +=  fm*(-s1.y+s2.y) + fp*relFac*(-s1.y-s2.y);
	      output[outputCount+1].y +=  fm*(s1.x-s2.x) + fp*relFac*(s1.x+s2.x);
	      output[outputCount+2].x += fm*(s1.x+s2.x) + fp*relFac*(s1.x-s2.x);
	      output[outputCount+2].y += fm*(s1.y+s2.y) + fp*relFac*(s1.y-s2.y); 
	      output[outputCount+3].x += relFac*(-s0.y + split*s3.y);
	      output[outputCount+3].y += relFac*(s0.x - split*s3.x); 
		  
	      c++;  
	    }
	    relFac *= -1;	  
	  }

          count += blockSize;
          outputCount +=4;
	}
	count += xtrAdd3;
      }
      count += xtrAdd2;
    }
    count += xtrAdd1;
  }  
  addPerformanceProfilingItem("ThreadControl::perform_MultiplicationVectorWithDerivativesOfB", performanceProfilerStartCycle, ThreadID);
}
   

void ThreadControl::perform_MultiCombinedOperator_VA1_RorCopy_VAc_RorCopy_D_InFourierSpace() {
  printf("ThreadControl::perform_MultiCombinedOperator_VA1_RorCopy_VAc_RorCopy_D_InFourierSpace: not implemented!!!\n");
  exit(0);
/*  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(ThreadID);
  DistributedMemoryObject* m1 = (DistributedMemoryObject*) ExecutionParameters[0];
  DistributedMemoryObject* m2 = (DistributedMemoryObject*) ExecutionParameters[1];
  DistributedMemoryObject* m3 = (DistributedMemoryObject*) ExecutionParameters[2];
  DistributedMemoryObject* m4 = (DistributedMemoryObject*) ExecutionParameters[3];
  DistributedMemoryObject* m5 = (DistributedMemoryObject*) ExecutionParameters[4];
  DistributedMemoryObject* m6 = (DistributedMemoryObject*) ExecutionParameters[5];
  DistributedMemoryObject* m7 = (DistributedMemoryObject*) ExecutionParameters[6];
  DistributedMemoryObject* m8 = (DistributedMemoryObject*) ExecutionParameters[7];

  Complex* v1 = m1->dataSets[ThreadID];
  Complex* v2 = m2->dataSets[ThreadID];
  Complex* v3 = m3->dataSets[ThreadID];
  Complex* v4 = m4->dataSets[ThreadID];
  Complex* sinP = m5->dataSets[ThreadID];
  Complex* auxData = m6->dataSets[ThreadID];
  Complex* RmulData1 = m7->dataSets[ThreadID];
  Complex* RmulData2 = m8->dataSets[ThreadID];
  
  int* intParas = (int*) ExecutionParameters[8];
  int L0 = intParas[0];
  int L1 = intParas[1];
  int L2 = intParas[2];
  int L3 = intParas[3];
  int blockSize = intParas[4];
  bool useR = (bool) intParas[5];
  int opMode = (bool) intParas[6];
  Complex alpha2 = *((Complex*) ExecutionParameters[9]);
  long int* syncData = (long int*) (ExecutionParameters[10]);
  
  if (LogLevel>4) printf("perform_MultiCombinedOperator_VA1_RorCopy_VAc_RorCopy_D_InFourierSpace called on thread %d on node %d with L0=%d, L1=%d, L2=%d, L3=%d, blocksize=%d, opMode=%d, useR=%d, and al.x=%f, al.y=%f\n", ThreadID, NodeID, L0, L1, L2, L3, blockSize, useR, opMode, alpha2.x, alpha2.y);
  
  int syncCount = 1;
  if (opMode==1) syncCount = 2; 
  if (opMode==2) syncCount = 4;   

  int countLocali2 = 0;
  int countLocali1 = 0;
  int countGlobal = 0;
  int RmulCountLocal = 0;
  int RmulCountGlobal = 0;
  int xtrAdd1 = xtraSize1*blockSize*(L3+xtraSize3)*(L2+xtraSize2);
  int xtrAdd2 = xtraSize2*blockSize*(L3+xtraSize3);
  int xtrAdd3 = xtraSize3*blockSize;
  int blockSizeHalf = blockSize/2;
  
  int mulCount = 0;
  int i2CombineCount = (L1CacheSizeInBytes / (16*blockSize*(L3+xtraSize3))) / 4;
  int i1CombineCount = (L2CacheSizeInBytes / (16*blockSize*(L3+xtraSize3)*(L2+xtraSize2))) / 4;
 

  for (int i0=0; i0<L0; i0++) {
    for (int i1Combined=0; i1Combined<L1; i1Combined+=i1CombineCount) {

      countLocali1 = countGlobal;
      for (int i1=i1Combined; ((i1<i1Combined+i1CombineCount) && (i1<L1)); i1++) {
        for (int i2Combined=0; i2Combined<L2; i2Combined+=i2CombineCount) {

          //Vector-Addition without Complex-Factor, ie alpha = 1
          countLocali2 = countLocali1;
          for (int i2=i2Combined; ((i2<i2Combined+i2CombineCount) && (i2<L2)); i2++) {
            for (int i3=0; i3<L3; i3++) {
              for (int i=0; i<blockSize; i++) {
                v3[countLocali2].x = v3[countLocali2].x + v1[countLocali2].x;
                v3[countLocali2].y = v3[countLocali2].y + v1[countLocali2].y;
		countLocali2++;
              } 
	    }
  	    countLocali2 += xtrAdd3;
	  }
	  
	  //Multiplication with R
          countLocali2 = countLocali1;
	  RmulCountLocal = RmulCountGlobal;
          for (int i2=i2Combined; ((i2<i2Combined+i2CombineCount) && (i2<L2)); i2++) {
            for (int i3=0; i3<L3; i3++) {
              for (int i=0; i<blockSizeHalf; i++) {
  	        //Imagin�rteil von mulData enth�t Faktor fr bottom-quark.
                v3[countLocali2].x = v3[countLocali2].x * RmulData1[RmulCountLocal].x;
                v3[countLocali2].y = v3[countLocali2].y * RmulData1[RmulCountLocal].x;
		countLocali2++;
	      }
              for (int i=0; i<blockSizeHalf; i++) {
  	        //Imagin�rteil von mulData enth�t Faktor fr bottom-quark.
                v3[countLocali2].x = v3[countLocali2].x * RmulData1[RmulCountLocal].y;
                v3[countLocali2].y = v3[countLocali2].y * RmulData1[RmulCountLocal].y;
		countLocali2++;
	      }
	      RmulCountLocal++;
	    }
  	    countLocali2 += xtrAdd3;
          }

          //Vector-Addition without Complex-Factor alpha2
          countLocali2 = countLocali1;
          for (int i2=i2Combined; ((i2<i2Combined+i2CombineCount) && (i2<L2)); i2++) {
            for (int i3=0; i3<L3; i3++) {
              for (int i=0; i<blockSize; i++) {
                v3[countLocali2].x = v3[countLocali2].x + (alpha2.x*v4[countLocali2].x - alpha2.y*v4[countLocali2].y);
                v3[countLocali2].y = v3[countLocali2].y + (alpha2.y*v4[countLocali2].x + alpha2.x*v4[countLocali2].y);
		countLocali2++;
              } 
	    }
  	    countLocali2 += xtrAdd3;
	  }


          if (useR) {
  	    //Multiplication with R
            countLocali2 = countLocali1;
	    RmulCountLocal = RmulCountGlobal;
            for (int i2=i2Combined; ((i2<i2Combined+i2CombineCount) && (i2<L2)); i2++) {
              for (int i3=0; i3<L3; i3++) {
                for (int i=0; i<blockSizeHalf; i++) {
  	          //Imagin�rteil von mulData enth�t Faktor fr bottom-quark.
                  v1[countLocali2].x = v3[countLocali2].x * RmulData2[RmulCountLocal].x;
                  v1[countLocali2].y = v3[countLocali2].y * RmulData2[RmulCountLocal].x;
	  	  countLocali2++;
	        }
                for (int i=0; i<blockSizeHalf; i++) {
  	          //Imagin�rteil von mulData enth�t Faktor fr bottom-quark.
                  v1[countLocali2].x = v3[countLocali2].x * RmulData2[RmulCountLocal].y;
                  v1[countLocali2].y = v3[countLocali2].y * RmulData2[RmulCountLocal].y;
 		  countLocali2++;
	        }
	        RmulCountLocal++;
	      }
  	      countLocali2 += xtrAdd3;
            }
	  } else {
  	    //Pure Copying
            countLocali2 = countLocali1;
            for (int i2=i2Combined; ((i2<i2Combined+i2CombineCount) && (i2<L2)); i2++) {
              for (int i3=0; i3<L3; i3++) {
                for (int i=0; i<blockSize; i++) {
                  v1[countLocali2].x = v3[countLocali2].x;
                  v1[countLocali2].y = v3[countLocali2].y;
	  	  countLocali2++;
	        }
	      }
  	      countLocali2 += xtrAdd3;
            }
	  }
	  
	  RmulCountGlobal = RmulCountLocal;
	  countLocali1 = countLocali2;
	}
	
	countLocali1 += xtrAdd2;
      }

      //Wait for synchronization
      xSSEObj->xSSE_WriteLongIntToMemAddr_Wrapper(&(syncData[ThreadID]), 1+i1Combined);
      while (true) {
        bool ready = true;
	
        for (int I=0; I<syncCount; I++) {
	  long int val = 0;
          xSSEObj->xSSE_ReadLongIntFromMemAddr_Wrapper(&(syncData[I]), val);
	  if (val != 1+i1Combined) ready = false;
        }
	
	if (ready) break;
	xSSEObj->xSSE_PerformNopLoop_Wrapper(ThreadControl_SYNCDELAYLOOPS);
      }
      
      //perform Dirac-operator
      





      countGlobal = countLocali1;
    }
    
    countGlobal += xtrAdd1;
  }  

  addPerformanceProfilingItem("ThreadControl::perform_MultiCombinedOperator_VA1_RorCopy_VAc_RorCopy_D_InFourierSpace", performanceProfilerStartCycle, ThreadID);
*/
}



void ThreadControl::execute() {
  if (ExecutionCommand == ThreadControl_ExecutionCommand_terminateThreads) {
    setTerminationFlag(1, xSSEObj);      
  }
  if (ExecutionCommand == ThreadControl_ExecutionCommand_allocateDistibutedMemory) {
    allocateDistributedMemoryObject();
  }
  if (ExecutionCommand == ThreadControl_ExecutionCommand_scalarProductOfFermionVectors) {
    scalarProductOfFermionVectors();
  }
  if (ExecutionCommand == ThreadControl_ExecutionCommand_vectorNormOfFermionVector) {
    vectorNormOfFermionVector();
  }
  if (ExecutionCommand == ThreadControl_ExecutionCommand_vectorAdditionOfFermionVectors) {
    vectorAdditionOfFermionVectors();
  }
  if (ExecutionCommand == ThreadControl_ExecutionCommand_perform_FFTWFourierTransformationOfFermionVector) {
    perform_FFTWFourierTransformationOfFermionVector();
  }
  if (ExecutionCommand == ThreadControl_ExecutionCommand_copyFermionVector) {
    copyFermionVector();
  }
  if (ExecutionCommand == ThreadControl_ExecutionCommand_perform_PiModeRemoverOperatorMultiplicationInFourierSpace) {
    perform_PiModeRemoverOperatorMultiplicationInFourierSpace();
  }
  if (ExecutionCommand == ThreadControl_ExecutionCommand_perform_DiracTypeOperatorMultiplicationInFourierSpace) {
    perform_DiracTypeOperatorMultiplicationInFourierSpace();
  }
  if (ExecutionCommand == ThreadControl_ExecutionCommand_perform_yBOperatorMultiplication) {
    perform_yBOperatorMultiplication();
  }
  if (ExecutionCommand == ThreadControl_ExecutionCommand_perform_yBDaggerOperatorMultiplication) {
    perform_yBDaggerOperatorMultiplication();
  }
  if (ExecutionCommand == ThreadControl_ExecutionCommand_perform_RPreconditionerTypeMatrixMultiplicationInFourierSpace) {
    perform_RPreconditionerTypeMatrixMultiplicationInFourierSpace();
  }
  if (ExecutionCommand == ThreadControl_ExecutionCommand_perform_MultiplicationVectorWithDerivativesOfB) {
    perform_MultiplicationVectorWithDerivativesOfB();
  } 
  if (ExecutionCommand == ThreadControl_ExecutionCommand_multiplyFermionVectorWithComplexNumber) {
    multiplyFermionVectorWithComplexNumber();
  } 
  
  if (ExecutionCommand == ThreadControl_ExecutionCommand_zeroFermionVector) {
    zeroFermionVector();
  } 
  if (ExecutionCommand == ThreadControl_ExecutionCommand_zeroUniformlyDistributedComplexVector) {
    zeroUniformlyDistributedComplexVector();
  } 
  if (ExecutionCommand == ThreadControl_ExecutionCommand_vectorAdditionOfUniformlyDistributedComplexVectors) {
    vectorAdditionOfUniformlyDistributedComplexVectors();
  } 
  if (ExecutionCommand == ThreadControl_ExecutionCommand_perform_xFFTFourierTransformationOfFermionVector) {
    perform_xFFTFourierTransformationOfFermionVector();
  } 
  if (ExecutionCommand == ThreadControl_ExecutionCommand_tune_xFFTFourierTransformationOfFermionVector) {
    tune_xFFTFourierTransformationOfFermionVector();
  } 
  if (ExecutionCommand == ThreadControl_ExecutionCommand_setFFTPlan_xFFTFourierTransformationOfFermionVector) {
    setFFTPlan_xFFTFourierTransformationOfFermionVector();
  } 
  if (ExecutionCommand == ThreadControl_ExecutionCommand_perform_MultiCombinedOperator_VA1_RorCopy_VAc_RorCopy_D_InFourierSpace) {
    perform_MultiCombinedOperator_VA1_RorCopy_VAc_RorCopy_D_InFourierSpace();
  } 
  if (ExecutionCommand == ThreadControl_ExecutionCommand_multiplyFermionVectorWithTimeIndexedComplexScalars) {
    multiplyFermionVectorWithTimeIndexedComplexScalars();
  } 
  
  setExecutionFlag(0, xSSEObj);    
}
