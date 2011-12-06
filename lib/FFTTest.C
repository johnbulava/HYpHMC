#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <iostream>

#include "Global.h"
#include FFTWIncludeFile
#include "Complex.h"
#include "ComplexVector.h"
#include "ExtremeFFT4D.h"
#include "Tools.C"
#include "xSSE.h"
#include "SSEroutines.C"
#include "MultiThreadedOperations.h"



long int CycleCounter;
int L0,L1,L2,L3;
Complex* input;
Complex* output;
int usexFFT = true;
int optimizeEmbedding = 0;
int threadCountPerNode = 1;
int ParaOpMode = 0;
int innerIndex = 8;
int tuneLevel = 1;
fftw_plan fftwPlan = NULL;
char* fileName = NULL;
ExtremeFFT4D* xFFT;


void startCycleCounter() {
  CycleCounter = getCPUCycleCounter();
}


long int cyclesPassed() {
  long int cycles = getCPUCycleCounter()-CycleCounter;
  return cycles;
}


void randomInput() {
  for (int I=0; I<8*40*40*40*40; I++) {
    input[I].x = NaN;
    input[I].y = NaN;    
  }
  for (int I=0; I<8*40*40*40*40; I++) {
    input[I] = Complex(2*(zufall()-0.5),2*(zufall()-0.5));
  }
}


void zeroOutput() {
  for (int I=0; I<8*40*40*40*40; I++) {
    output[I].x = 0;
    output[I].y = 0;    
  }
}


void generateFileName(char* tag) {
  delete[] fileName;
  fileName = new char[1000];
  snprintf(fileName,1000,"FFTTest_%s_XFFT%dOptiEmbed%dThreads%dParaOp%dInnerInd%dTuneLevel%d_%s.dat", getHostName(), usexFFT, optimizeEmbedding, threadCountPerNode, ParaOpMode, innerIndex, tuneLevel, tag);
}


void performFFT() {
  bool forw = true;
  if (usexFFT) {
    xFFT->FastFourierTrafo(input, output, forw);
  
  } else {
    fftw_execute(fftwPlan);   
  }
}


void setHyperCubicVolume(int L) {
  L0 = L;
  L1 = L;
  L2 = L;
  L3 = L;
}


void set1DVolume(int L) {
  L0 = L;
  L1 = 1;
  L2 = 1;
  L3 = 1;
}


long int calcTheoreticalFlops() {
  long int cycles = (long int) ((log(L3)/log(2)) * 8*L0*L1*L2*L3* 4);
  return cycles;
}


void tuneFFT() {
  bool forw = true;
  if (LogLevel>2) printf("Generating FFT-plan. Forward = %d\n", forw);  

  if (optimizeEmbedding==0) {
    xtraSize1 = 0;
    xtraSize2 = 0;
    xtraSize3 = 0;
  }
  if (optimizeEmbedding==1) {
    char* fftPlanDescriptor = NULL;  
    readOptimalFermionVectorEmbeddingAndFFTPlanFromTuningDB(L0, L1, L2, L3, threadCountPerNode, ParaOpMode, usexFFT, 1, 1, 1, 1, fftPlanDescriptor);
    delete[] fftPlanDescriptor;
  } 

  if (usexFFT) {
    delete xFFT;  
    
    long int* threadCoreMsk = new long int[threadCountPerNode];
    for (int I=0; I<threadCountPerNode; I++) threadCoreMsk[I] = getAffinityMaskFromNodeID(0);      
    
    xFFT = new ExtremeFFT4D(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, innerIndex, threadCountPerNode, threadCoreMsk, true);  
    xFFT->tune(tuneLevel);
    delete[] threadCoreMsk;
  } else {
    int* n = new int[4];
    n[0] = L0;
    n[1] = L1;
    n[2] = L2;
    n[3] = L3;
  
    int rank = 4;
    int howmany = innerIndex;
    int* inembed = new int[4];
    inembed[0] = L0;
    inembed[1] = L1+xtraSize1;
    inembed[2] = L2+xtraSize2;
    inembed[3] = L3+xtraSize3;  
    int istride = innerIndex;
    int idist = 1;
  
    int* onembed = new int[4];
    onembed[0] = L0;
    onembed[1] = L1+xtraSize1;
    onembed[2] = L2+xtraSize2;
    onembed[3] = L3+xtraSize3;    
    int ostride = innerIndex;
    int odist = 1;

    int measureFlag = FFTW_MEASURE | FFTW_EXHAUSTIVE;
    if (tuneLevel==0) measureFlag = FFTW_MEASURE;
    fftw_destroy_plan(fftwPlan);
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
    delete[] n;
    delete[] inembed;
    delete[] onembed;
  }
}



void tune1DFFT() {
  bool forw = true;
  if (LogLevel>2) printf("Generating 1D-FFT-plan. Forward = %d\n", forw);  

  xtraSize1 = 0;
  xtraSize2 = 0;
  xtraSize3 = 0;

  if (usexFFT) {
    delete xFFT;  
    
    printf("1D xFFT not implemented!!!\n");
    exit(0);
  } else {
    int* n = new int[1];
    n[0] = L0;
  
    int rank = 1;
    int howmany = innerIndex;
    int* inembed = new int[1];
    inembed[0] = L0;
    int istride = innerIndex;
    int idist = 1;
  
    int* onembed = new int[1];
    onembed[0] = L0;
    int ostride = innerIndex;
    int odist = 1;

    int measureFlag = FFTW_MEASURE | FFTW_EXHAUSTIVE;
    if (tuneLevel==0) measureFlag = FFTW_MEASURE;
    fftw_destroy_plan(fftwPlan);
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
    delete[] n;
    delete[] inembed;
    delete[] onembed;
  }
}



void FourierTrafo4DVolumeScan(int iter) {
  printf("\nVolume Scan of 4D-Fast-Fourier-Transformation ...\n");
  setCurrentThreadAffinityToNodeID(0);
  randomInput();
  zeroOutput();
  generateFileName("VolumeScan");
  FILE* file = fopen(fileName, "w");  
  fprintf(file, "#Volume Scan\n");
  fclose(file);
  
  
  for (int L=4; L<=32; L*=2) {
    setHyperCubicVolume(L);
    int xsmax1 = 1;
    int xsmax2 = 1;
    int xsmax3 = 1;
    if (optimizeEmbedding>1) {
      xsmax1 = optimizeEmbedding;
      xsmax2 = optimizeEmbedding;
      xsmax3 = optimizeEmbedding;
    }
    long int optimal_Cycles = -1;
    long int optimal_xs1 = 0;
    long int optimal_xs2 = 0;
    long int optimal_xs3 = 0;
    
    for (int xs1=0; xs1<xsmax1; xs1++) {
      for (int xs2=0; xs2<xsmax2; xs2++) {
        for (int xs3=0; xs3<xsmax3; xs3++) {
          xtraSize1 = xs1;
          xtraSize2 = xs2;
          xtraSize3 = xs3;
    
          printf("Tuning for Volume: %d with xtrSize (%d, %d ,%d)\n", L, xtraSize1, xtraSize2, xtraSize3);
          tuneFFT();
          performFFT();
    
          printf("Start measuring...");
          startCycleCounter();
          for (int I=0; I<iter; I++) {
            performFFT();    
          }
          long int cycles = cyclesPassed() / iter;
          long int theorCycles = calcTheoreticalFlops();
          printf("   %ld cycles on average, and %ld theoretically: %1.2f percent performance.\n", cycles, theorCycles, 100.0*theorCycles/cycles);
          file = fopen(fileName, "a");  
	  if (optimizeEmbedding>1) fprintf(file, "#");
          fprintf(file, "%d %ld %d %d %d\n", L, cycles, xtraSize1, xtraSize2, xtraSize3);
          fclose(file);
	  
	  if ((optimal_Cycles<0) || (optimal_Cycles>cycles)) {
  	    optimal_Cycles = cycles;
	    optimal_xs1 = xtraSize1;
	    optimal_xs2 = xtraSize2;
	    optimal_xs3 = xtraSize3;	    
	  }
	}
      }
    }
    if (optimizeEmbedding>1) {
      file = fopen(fileName, "a");  
      fprintf(file, "%d %ld %d %d %d\n", L, optimal_Cycles, optimal_xs1, optimal_xs2, optimal_xs3);
      fclose(file);
    }
    iter /= 16;
  }
}



void FourierTrafo4DParaScan(int iter, int L) {
  printf("\nPara Scan of 4D-Fast-Fourier-Transformation ...\n");
  randomInput();
  zeroOutput();
  generateFileName("ParaScan");
  FILE* file = fopen(fileName, "w");  
  fprintf(file, "#Para Scan\n");
  fclose(file);
  
  
  setHyperCubicVolume(L);
  xtraSize1 = 0;
  xtraSize2 = 0;
  xtraSize3 = 0;

  printf("Tuning for Volume: %d with xtrSize (%d, %d ,%d)\n", L, xtraSize1, xtraSize2, xtraSize3);
  tuneFFT();
  performFFT();
    
  printf("Start measuring...");
  startCycleCounter();
  for (int I=0; I<iter; I++) {
    performFFT();    
  }
  long int cycles = cyclesPassed() / iter;
  long int theorCycles = calcTheoreticalFlops();
  printf("   %ld cycles on average, and %ld theoretically: %1.2f percent performance.\n", cycles, theorCycles, 100.0*theorCycles/cycles);
  file = fopen(fileName, "a");  
  fprintf(file, "1 1 %d %ld %d %d %d\n", L, cycles, xtraSize1, xtraSize2, xtraSize3);
  fclose(file);
  
  MultiThreadedOperations* threadedOps = NULL;
  DistributedMemoryObject* memObj1 = NULL;
  DistributedMemoryObject* memObj2 = NULL;

  for (int I=1; I<5; I++) {
    int paraOp = 1;
    int threadCNT = 1;
    int nodeCNT = 2;
    if (I==1) {
      paraOp = 1;
      nodeCNT = 2;
      threadCNT = 1;
    }
    if (I==2) {
      paraOp = 1;
      nodeCNT = 2;
      threadCNT = 2;
    }
    if (I==3) {
      paraOp = 2;
      nodeCNT = 4;
      threadCNT = 1;
    }
    if (I==4) {
      paraOp = 2;
      nodeCNT = 4;
      threadCNT = 2;
    }
    
    threadedOps = new MultiThreadedOperations(paraOp, false);
    memObj1 = threadedOps->allocateDistibutedFermionVector(L, L, L, L);
    memObj2 = threadedOps->allocateDistibutedFermionVector(L, L, L, L);
    threadedOps->copyFermionVectorToDistributedFermionVector(input, memObj1, L, L, L, L);
    threadedOps->copyFermionVectorToDistributedFermionVector(output, memObj2, L, L, L, L);

    printf("Tuning for Volume: %d with xtrSize (%d, %d ,%d) and modus %d\n", L, xtraSize1, xtraSize2, xtraSize3, I);
    if (usexFFT) {
      threadedOps->perform_xFFTFourierTransformationOfFermionVector(memObj1, memObj2, ExtremeFFT4D_Forward, threadCNT, L, L, L, L);
    } else {
      fftw_plan_with_nthreads(threadCNT);    
      threadedOps->perform_FFTWFourierTransformationOfFermionVector(memObj1, memObj2, ExtremeFFT4D_Forward, L, L, L, L);
    }

    printf("Start measuring...");
    startCycleCounter();
    if (usexFFT) {
      for (int I=0; I<iter; I++) {
        threadedOps->perform_xFFTFourierTransformationOfFermionVector(memObj1, memObj2, ExtremeFFT4D_Forward, threadCNT, L, L, L, L);
      }
    } else {
      for (int I=0; I<iter; I++) {
        threadedOps->perform_FFTWFourierTransformationOfFermionVector(memObj1, memObj2, ExtremeFFT4D_Forward, L, L, L, L);
      }
    }
    long int cycles = cyclesPassed() / iter;
    long int theorCycles = calcTheoreticalFlops();
    printf("   %ld cycles on average, and %ld theoretically: %1.2f percent performance.\n", cycles, theorCycles, 100.0*theorCycles/cycles);
    file = fopen(fileName, "a");  
    fprintf(file, "%d %d %d %ld %d %d %d\n", nodeCNT, threadCNT, L, cycles, xtraSize1, xtraSize2, xtraSize3);
    fclose(file);

    delete threadedOps;
    delete memObj1;
    delete memObj2;
  }
}


void FourierTrafo1DVolumeScan(int iter) {
  printf("\nVolume Scan of 1D-Fast-Fourier-Transformation ...\n");
  setCurrentThreadAffinityToNodeID(0);
  randomInput();
  zeroOutput();
  generateFileName("1DVolumeScan");
  FILE* file = fopen(fileName, "w");  
  fprintf(file, "#1DVolume Scan\n");
  fclose(file);
  
  
  for (int L=2; L<=262144; L*=2) {
    set1DVolume(L);

    xtraSize1 = 0;
    xtraSize2 = 0;
    xtraSize3 = 0;
    
    printf("Tuning for Volume: %d with xtrSize (%d, %d ,%d)\n", L, xtraSize1, xtraSize2, xtraSize3);
    tune1DFFT();
    performFFT();
    
    printf("Start measuring...");
    startCycleCounter();
    for (int I=0; I<iter; I++) {
      performFFT();    
    }
    long int cycles = cyclesPassed() / iter;
    long int theorCycles = calcTheoreticalFlops();
    printf("   %ld cycles on average, and %ld theoretically: %1.2f percent performance.\n", cycles, theorCycles, 100.0*theorCycles/cycles);
    file = fopen(fileName, "a");  
    fprintf(file, "%d %ld %d %d %d\n", L, cycles, xtraSize1, xtraSize2, xtraSize3);
    fclose(file);
	  
    iter /= 2;
  }
}






int main(int argc,char **argv) {
  LogLevel = 3;

  bool error = false;
  if (argc<7) { 
    printf("Parameters could not be read!!!\n");
    error = true;
    exit(0);  
  }
  if (sscanf(argv[1],"%d",&usexFFT)!=1) error = true;
  if (sscanf(argv[2],"%d",&optimizeEmbedding)!=1) error = true;
  if (sscanf(argv[3],"%d",&threadCountPerNode)!=1) error = true;
  if (sscanf(argv[4],"%d",&ParaOpMode)!=1) error = true;  
  if (sscanf(argv[5],"%d",&innerIndex)!=1) error = true;
  if (sscanf(argv[6],"%d",&tuneLevel)!=1) error = true;
  if (error) {
    printf("Parameters could not be read!!!\n");
    exit(0);  
  }

  printf("Read Parameters:\n");
  printf("  -> usexFFT:  %d\n",usexFFT);
  printf("  -> optimizeEmbedding:  %d\n",optimizeEmbedding);
  printf("  -> threadCountPerNode:  %d\n",threadCountPerNode);
  printf("  -> ParaOpMode:  %d\n",ParaOpMode);
  printf("  -> innerIndex:  %d\n",innerIndex);
  printf("  -> tuneLevel:  %d\n",tuneLevel);
  

  if (LogLevel>1) printf("Initialisiere FFTW with %d threads.\n",threadCountPerNode);
  fftw_init_threads();
  fftw_plan_with_nthreads(threadCountPerNode);  
  
  iniTools(1);
  setHyperCubicVolume(1);
  
  input = createSuperAlignedComplex(8*40*40*40*40);
  output = createSuperAlignedComplex(8*40*40*40*40);
  
  
  FourierTrafo1DVolumeScan(10*262144);
//  FourierTrafo4DVolumeScan(10*16*16*16);
//  FourierTrafo4DParaScan(10, 32);
  
  
  destroySuperAlignedComplex(input);
  destroySuperAlignedComplex(output);
  
  delete[] fileName;
  desiniTools();
  fftw_destroy_plan(fftwPlan);
  fftw_cleanup_threads();
  if (LogLevel>1) printf("FFTTest terminated correctly.\n");
}
