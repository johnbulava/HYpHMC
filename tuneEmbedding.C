#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>


#include "Complex.h"
#include "Quat.h"
#include "Global.h"
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include "Tools.h"
#include "FermionMatrixOperations.h"
#include "TuningDataBase.h"

#define scanTimeInSecs 100



//Variables
int Parameter_L0 = 0;
int Parameter_L1 = 0;
int Parameter_L2 = 0;
int Parameter_L3 = 0;
double Parameter_RHO = 0;
double Parameter_R = 0;
double Parameter_Y = 0;
int Parameter_ParaOpMode = 0;
int Parameter_ThreadsPerNodeMax = 0;
int Parameter_FLAG_xFFT = 0;
int Parameter_MAXXTRSIZE = 0;
int Parameter_FLAG_USE_P = 0;
int Parameter_FLAG_USE_Q = 0;
int Parameter_FLAG_USE_R = 0;
int Parameter_FLAG_USE_QHM = 0;
char* Parameter_HOSTNAME = NULL;
char* tuningDBFileName = NULL;
TuningDataBase* tuningDB = NULL;
FermionMatrixOperations* fermiOps = NULL;
double startTime;


void startTimer() {
  startTime = zeitwert();
}


double timePassed() {
  double time = (zeitwert()-startTime);
  return time;
}


int main(int argc,char **argv) {
  LogLevel = 0;
  fftw_init_threads();
  fftw_plan_with_nthreads(1);  

  if (LogLevel>1) printf("Number of arguments = %d\n",argc);
  for (int I=0; I<argc; I++) {
    if (LogLevel>1) printf("Argument %d: %s\n",I+1,argv[I]);  
  }
  
  if (argc<13) {
    printf("Lattice Size Parameters, Parallel-OpMode, ThreadsPerNode, Use_xFFT, Use_P, Use_Q, Use_R, UseQuasiHermitMode and MaxXtrSize required!!!\n");
    exit(0);
  }

  bool error = false;
  if (sscanf(argv[1],"%d",&Parameter_L0)!=1) error = true;
  if (sscanf(argv[2],"%d",&Parameter_L1)!=1) error = true;
  if (sscanf(argv[3],"%d",&Parameter_L2)!=1) error = true;
  if (sscanf(argv[4],"%d",&Parameter_L3)!=1) error = true;  
  if (sscanf(argv[5],"%d",&Parameter_ParaOpMode)!=1) error = true;
  if (sscanf(argv[6],"%d",&Parameter_ThreadsPerNodeMax)!=1) error = true;
  if (sscanf(argv[7],"%d",&Parameter_FLAG_xFFT)!=1) error = true;
  if (sscanf(argv[8],"%d",&Parameter_FLAG_USE_P)!=1) error = true;
  if (sscanf(argv[9],"%d",&Parameter_FLAG_USE_Q)!=1) error = true;  
  if (sscanf(argv[10],"%d",&Parameter_FLAG_USE_R)!=1) error = true;  
  if (sscanf(argv[11],"%d",&Parameter_FLAG_USE_QHM)!=1) error = true;  
  if (sscanf(argv[12],"%d",&Parameter_MAXXTRSIZE)!=1) error = true;
  if (error) {
    printf("Parameters could not be read!!!\n");
    exit(0);  
  }
  
  if (Parameter_L0<=0) error = true;
  if (Parameter_L1<=0) error = true;
  if (Parameter_L2<=0) error = true;
  if (Parameter_L3<=0) error = true;
  if ((Parameter_ParaOpMode<0) || (Parameter_ParaOpMode>2))  error = true;
  if (Parameter_ThreadsPerNodeMax<=0) error = true;
  if ((Parameter_FLAG_xFFT<0) || (Parameter_FLAG_xFFT>1)) error = true;
  if ((Parameter_FLAG_USE_P<0) || (Parameter_FLAG_USE_P>1)) error = true;
  if ((Parameter_FLAG_USE_Q<0) || (Parameter_FLAG_USE_Q>1)) error = true;
  if ((Parameter_FLAG_USE_R<0) || (Parameter_FLAG_USE_R>1)) error = true;
  if ((Parameter_FLAG_USE_QHM<0) || (Parameter_FLAG_USE_QHM>1)) error = true;
  if (Parameter_MAXXTRSIZE<0) error = true;
  if (error) {
    printf("Invalid Parameter setting !!!\n");
    exit(0);  
  }
  
  if (DebugMode) {
    printf("ERROR: System is in debug-mode !!!\n");
    exit(0);    
  }
  
  printf("Parameters: \n");
  printf("  --> L0                   = %d\n", Parameter_L0);
  printf("  --> L1                   = %d\n", Parameter_L1);
  printf("  --> L2                   = %d\n", Parameter_L2);
  printf("  --> L3                   = %d\n", Parameter_L3);
  printf("  --> ParaOpMode           = %d\n", Parameter_ParaOpMode);
  printf("  --> ThreadsPerNodeMax    = %d\n", Parameter_ThreadsPerNodeMax);
  printf("  --> xFFT                 = %d\n", Parameter_FLAG_xFFT);
  printf("  --> Use_P                = %d\n", Parameter_FLAG_USE_P);
  printf("  --> Use_Q                = %d\n",Parameter_FLAG_USE_Q);
  printf("  --> Use_R                = %d\n",Parameter_FLAG_USE_R);
  printf("  --> Use_QHM              = %d\n",Parameter_FLAG_USE_QHM);
  printf("  --> MaxXtrSize           = %d\n", Parameter_MAXXTRSIZE);
  printf("\n\n");

  iniTools(5517);
  
  Parameter_RHO = 1.0;
  Parameter_R = 0.5;
  Parameter_Y = 1.0;
  Parameter_HOSTNAME = getHostName();
  printf("Benchmark is running on: %s\n",Parameter_HOSTNAME);

  tuningDBFileName = getTuningDBFileName(Parameter_L0,Parameter_L1,Parameter_L2,Parameter_L3, (bool)Parameter_FLAG_xFFT);
  printf("Tuning-DB is %s\n", tuningDBFileName);
  tuningDB = new TuningDataBase(tuningDBFileName);
  delete[] tuningDBFileName;
  
  char* performanceProfileFileName = getTuningPerformanceProfileFileName(Parameter_L0,Parameter_L1,Parameter_L2,Parameter_L3, (bool)Parameter_FLAG_xFFT);
  printf("Tuning-Performance-Profile is %s\n", performanceProfileFileName);
  initializePerformanceProfiler(performanceProfileFileName);
  delete[] performanceProfileFileName;
  
  for (xtraSize1=0; xtraSize1<1+Parameter_MAXXTRSIZE; xtraSize1++) {
    for (xtraSize2=0; xtraSize2<1+Parameter_MAXXTRSIZE; xtraSize2++) {
      for (xtraSize3=0; xtraSize3<1+Parameter_MAXXTRSIZE; xtraSize3++) {
        for (int threadsCount=1; threadsCount<=Parameter_ThreadsPerNodeMax; threadsCount++) {
	  if (!tuningDB->isDBEntryAvail(Parameter_HOSTNAME, numaCoresCount, Parameter_L0,Parameter_L1,Parameter_L2,Parameter_L3, Parameter_ParaOpMode, threadsCount, Parameter_FLAG_xFFT, Parameter_FLAG_USE_P, Parameter_FLAG_USE_Q, Parameter_FLAG_USE_R, Parameter_FLAG_USE_QHM, xtraSize1, xtraSize2, xtraSize3)) {
            fftw_plan_with_nthreads(threadsCount);  
            fermiOps = new FermionMatrixOperations(Parameter_L0,Parameter_L1,Parameter_L2,Parameter_L3, Parameter_RHO, Parameter_R, Parameter_Y);
            fermiOps->setxFFTusage(Parameter_FLAG_xFFT);
	    fermiOps->setxFFT_DistributedFFT_ThreadCount(threadsCount);
	    fermiOps->setPreconditioner((bool) Parameter_FLAG_USE_P, 1.1, 0);
            fermiOps->setQPreconditioner((bool) Parameter_FLAG_USE_Q, 0.25, 0.25);
  	    fermiOps->setRPreconditioner((bool) Parameter_FLAG_USE_R, 1.0, 1.0);
            fermiOps->activateMultiThreadedOps(Parameter_ParaOpMode, false);  
	  

            Complex* input = fermiOps->createFermionVector();
            Complex* output = fermiOps->createFermionVector();
            Complex* phiFieldComplex = fermiOps->createFermionVector();
	    double* phiField = (double*)phiFieldComplex;
	
  	    for (int I=0; I<fermiOps->getVectorLengthXtrSize(); I++) {
    	      input[I].x = 1.0;
    	      input[I].y = 1.0;	  
  	    }
	    for (int I=0; I<4*Parameter_L0*Parameter_L1*Parameter_L2*Parameter_L3; I++) {
	      phiField[I] = 1.0;
  	    }
            DistributedMemoryObject* memObj1 = fermiOps->createDistributedFermionVector();
            DistributedMemoryObject* memObj2 = fermiOps->createDistributedFermionVector();
            DistributedMemoryObject* memObj3 = fermiOps->createDistributedUniformLatticeComplexVector(2); 
            fermiOps->copyFermionVectorToDistributedFermionVector(input, memObj1);
            fermiOps->copyFermionVectorToDistributedFermionVector(output, memObj2);
            fermiOps->copyPhiFieldUniformlyToDistributedPhiField(phiField, memObj3);

            //Give FFTW time to optimize...
	    char* fftPlanDescriptor = NULL;
	    if (Parameter_FLAG_USE_QHM==1) {
              fermiOps->tuneDistributedFFT(memObj1, memObj2, ExtremeFFT4D_TuneLevel_Low, fftPlanDescriptor);
	      fermiOps->executeDistributedFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(memObj1, memObj2, memObj3, true, true);
	    } else {
	      printf("executeFermionMatrixFermionDaggerMatrixMultiplication not implemented in QHM-Mode!!!\n");
	      exit(0);
	    }
	    if (fftPlanDescriptor == NULL) {
	      fftPlanDescriptor = new char[100];
	      snprintf(fftPlanDescriptor, 100, " ");
	    }

            printf("Testing %dx(%d+%d)x(%d+%d)x(%d+%d) threadsPerNode=%d POM=%d xFFT=%d P=%d Q=%d R=%d QHM=%d", Parameter_L0,Parameter_L1,xtraSize1,Parameter_L2,xtraSize2,Parameter_L3,xtraSize3,threadsCount, Parameter_ParaOpMode, Parameter_FLAG_xFFT,Parameter_FLAG_USE_P, Parameter_FLAG_USE_Q,Parameter_FLAG_USE_R,Parameter_FLAG_USE_QHM);
    	    startTimer();
	    int runCount = 0;
  	    while (timePassed() < scanTimeInSecs) {
	      for (int I=0; I<10; I++) {
  	        if (Parameter_FLAG_USE_QHM==1) {
  	          fermiOps->executeDistributedFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(memObj1, memObj2, memObj3, true, true);
 	        } else {
	          printf("executeFermionMatrixFermionDaggerMatrixMultiplication not implemented in QHM-Mode!!!\n");
	          exit(0);
	        }
	        runCount++;
	      }
	    }
	    double deltaT = timePassed();
	    printf(" ==> Runs=%d in %f secs. ==> speed=%f/sec.\n",runCount, deltaT,runCount/deltaT);

            tuningDB->writeDBEntry(Parameter_HOSTNAME, numaCoresCount, Parameter_L0,Parameter_L1,Parameter_L2,Parameter_L3, Parameter_ParaOpMode, threadsCount, Parameter_FLAG_xFFT, Parameter_FLAG_USE_P, Parameter_FLAG_USE_Q, Parameter_FLAG_USE_R, Parameter_FLAG_USE_QHM, xtraSize1, xtraSize2, xtraSize3, runCount,deltaT,fftPlanDescriptor);
            writePerformanceProfilingDataToDisk();
	    
            delete memObj1;
            delete memObj2;
            delete memObj3;
	    delete[] fftPlanDescriptor;
            fermiOps->deactivateMultiThreadedOps();

	    fermiOps->destroyFermionVector(input);
  	    fermiOps->destroyFermionVector(output);
	    fermiOps->destroyFermionVector(phiFieldComplex);
            delete fermiOps;
	    fermiOps = NULL;
	  }
	}
      }
    }
  }
  delete tuningDB;

  LogLevel = 1;
  printf("\n\n");
  
  int threadCountPerNode = 1;
  char* fftPlanDescriptor = NULL;
  readOptimalFermionVectorEmbeddingAndFFTPlanFromTuningDB(Parameter_L0, Parameter_L1, Parameter_L2, Parameter_L3, threadCountPerNode, Parameter_ParaOpMode, Parameter_FLAG_xFFT, Parameter_FLAG_USE_P, Parameter_FLAG_USE_Q, Parameter_FLAG_USE_R, Parameter_FLAG_USE_QHM, fftPlanDescriptor);
  delete[] fftPlanDescriptor;

  delete[] Parameter_HOSTNAME;
  desiniTools();
}
