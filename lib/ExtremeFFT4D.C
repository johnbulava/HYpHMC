#include "ExtremeFFT4D.h"

ExtremeFFT4D::ExtremeFFT4D(int l0, int l1, int l2, int l3, int xtrS1, int xtrS2, int xtrS3, int locIndCnt, int nrThreads, long int* threadCoreMsk, bool runContin) {
  L[0] = l0;
  L[1] = l1;
  L[2] = l2;
  L[3] = l3;
  xtrSize1 = xtrS1;
  xtrSize2 = xtrS2;
  xtrSize3 = xtrS3;
  localIndexCount = locIndCnt;
  numberOfThreads = nrThreads;  
  CoreMaskForThreads = new long int[numberOfThreads];
  if (threadCoreMsk != NULL) {
    for (int I=0; I<numberOfThreads; I++) CoreMaskForThreads[I] = threadCoreMsk[I];  
  } else {
    for (int I=0; I<numberOfThreads; I++) CoreMaskForThreads[I] = 255;    
  }
  runThreadsContinously = runContin;
  xSSEObj = new xSSE();
  
  long int* BitMask_ComplexI1_Pointer = (long int*)&BitMask_ComplexI1;
  long int* BitMask_ComplexI2_Pointer = (long int*)&BitMask_ComplexI2;
  BitMask_ComplexI1_Pointer[0] = (long int) 0x8000000000000000;
  BitMask_ComplexI1_Pointer[1] = (long int) 0;
  BitMask_ComplexI2_Pointer[0] = (long int) 0;
  BitMask_ComplexI2_Pointer[1] = (long int) 0x8000000000000000;

  if (LogLevel>1) printf("Initializing ExtremeFFT4D with parameters L0: %d, L1:%d, L2:%d, L3:%d, xtrSize1: %d, xtrSize2: %d, xtrSize3: %d, locIndCnt: %d, nrThreads: %d\n", L[0], L[1], L[2], L[3], xtrSize1, xtrSize2, xtrSize3, localIndexCount, numberOfThreads);

  for (int I=0; I<4; I++) {
    if (L[I] <= 0) {
      printf("ExtremeFFT4D: L%d not positive!\n", I);
      exit(0);  
    }
  }

  addressIncrementsInBytesL[3] = 16*localIndexCount;
  addressIncrementsInBytesL[2] = 16*localIndexCount*(L[3]+xtrSize3);
  addressIncrementsInBytesL[1] = 16*localIndexCount*(L[3]+xtrSize3)*(L[2]+xtrSize2);
  addressIncrementsInBytesL[0] = 16*localIndexCount*(L[3]+xtrSize3)*(L[2]+xtrSize2)*(L[1]+xtrSize1);
  
  for (int I=0; I<4; I++) {
    findPrimeFactors(L[I], primeFactorsL[I], primeFactorCountL[I]);
    if (LogLevel>0*4) printf("Prime factors of L%d: ", I);  
    for (int I2=0; I2<primeFactorCountL[I]; I2++) {      
      if (LogLevel>0*4) printf("%d ", primeFactorsL[I][I2]);
    }
    if (LogLevel>0*4) printf("\n");
  }

  LargestL = L[0];
  for (int I=0; I<4; I++) {
    if (L[I]>LargestL) LargestL = L[I];
  }
  
  Slices2D = new Complex*[numberOfThreads];
  ThreadControlDataStructure = new ThreadControlDataStructureType[numberOfThreads];
  ThreadControlParameterHolder = new long int*[numberOfThreads];
  for (int I=0; I<numberOfThreads; I++) {
    ThreadControlParameterHolder[I] = new long int[2];
    long int slice2Dsize = L1CacheSizePerWayInBytes/16 + ExtremeFFT4D_InnerLoopMAX*LargestL*LargestL + 4*L1CacheSizePerWayInBytes*LargestL/(4*16) + 4*LargestL*ExtremeFFT4D_InternalEmbeddingOneMAX;
    Slices2D[I] = createSuperAlignedComplex(slice2Dsize);

    ThreadControlDataStructure[I].running = 0;
    ThreadControlDataStructure[I].command = ExtremeFFT4D_ThreadControlCommand_NoCommand;
  }

  PlanExecutionControlDataStructure.SliceStatus = new long int[LargestL*LargestL];
  PlanExecutionControlDataStructure.WaitLoopEntered = new long int[numberOfThreads];
  
  availFFTstepCount = 5;
  availFFTsteps = new ExtremeFFT4D_FFTstep*[availFFTstepCount];  
  availFFTsteps[0] = new ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2();
  availFFTsteps[1] = new ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2();
  availFFTsteps[2] = new ExtremeFFT4D_FFTstep_Write2DSlice();
  availFFTsteps[3] = new ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2();
  availFFTsteps[4] = new ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2();


  possibleFFTplanList = new FFTplanType;
  FFTplanType* dummy = possibleFFTplanList;
  dummy->complete = false;
  dummy->LindA1 = 2;
  dummy->LindA2 = 3;
  dummy->LindB1 = 0;
  dummy->LindB2 = 1;
  dummy->inputSlice2DDistanceInBytesA = 0;
  dummy->inputSlice2DDistanceInBytesB = 0;  
  dummy->MinorMissesA = -1;
  dummy->MajorMissesA = -1;
  dummy->MinorMissesB = -1;
  dummy->MajorMissesB = -1;
  dummy->InternalEmbeddingAOne = 0;
  dummy->InternalEmbeddingATwo = 0;
  dummy->InternalEmbeddingBOne = 0;
  dummy->InternalEmbeddingBTwo = 0;  
  dummy->FFTstepCountA = 0;
  dummy->FFTstepCountB = 0;
  dummy->FFTstepsA = NULL;
  dummy->FFTstepsB = NULL;
  dummy->performance_TotalCycles = NaN;
  dummy->performance_FFTstepCyclesA = NULL;
  dummy->performance_FFTstepCyclesB = NULL;
  dummy->next = NULL;  
  generateAllPossibleFFTplans(dummy, 0);
  
  dummy->complete = false;
  dummy->LindA1 = 1;
  dummy->LindA2 = 3;
  dummy->LindB1 = 0;
  dummy->LindB2 = 2;
  dummy->inputSlice2DDistanceInBytesA = 0;
  dummy->inputSlice2DDistanceInBytesB = 0;  
  dummy->MinorMissesA = -1;
  dummy->MajorMissesA = -1;
  dummy->MinorMissesB = -1;
  dummy->MajorMissesB = -1;
  dummy->InternalEmbeddingAOne = 0;
  dummy->InternalEmbeddingATwo = 0;
  dummy->InternalEmbeddingBOne = 0;
  dummy->InternalEmbeddingBTwo = 0;  
  dummy->FFTstepCountA = 0;
  dummy->FFTstepCountB = 0;
  dummy->FFTstepsA = NULL;
  dummy->FFTstepsB = NULL;  
  dummy->performance_TotalCycles = NaN;
  dummy->performance_FFTstepCyclesA = NULL;
  dummy->performance_FFTstepCyclesB = NULL;
  dummy->next = NULL;  
  generateAllPossibleFFTplans(dummy, 0);
  
  dummy->complete = false;
  dummy->LindA1 = 0;
  dummy->LindA2 = 3;
  dummy->LindB1 = 1;
  dummy->LindB2 = 2;
  dummy->inputSlice2DDistanceInBytesA = 0;
  dummy->inputSlice2DDistanceInBytesB = 0;  
  dummy->MinorMissesA = -1;
  dummy->MajorMissesA = -1;
  dummy->MinorMissesB = -1;
  dummy->MajorMissesB = -1;
  dummy->InternalEmbeddingAOne = 0;
  dummy->InternalEmbeddingATwo = 0;
  dummy->InternalEmbeddingBOne = 0;
  dummy->InternalEmbeddingBTwo = 0;  
  dummy->FFTstepCountA = 0;
  dummy->FFTstepCountB = 0;
  dummy->FFTstepsA = NULL;
  dummy->FFTstepsB = NULL;  
  dummy->performance_TotalCycles = NaN;
  dummy->performance_FFTstepCyclesA = NULL;
  dummy->performance_FFTstepCyclesB = NULL;
  dummy->next = NULL;  
  generateAllPossibleFFTplans(dummy, 0);
  
  if (!possibleFFTplanList->complete) {
    printf("ExtremeFFT4D: No FFT-plan available for (%dx%dx%dx%d)!\n",L[0],L[1],L[2],L[3]);
    exit(0);
  }
  selectedFFTplan = possibleFFTplanList;    
  elaboratePlanDetails(selectedFFTplan);
  optimizeInternalEmbeddingForPlan(selectedFFTplan, ExtremeFFT4D_TuneLevel_None, true);
  
  getPlanReadyForExecution(selectedFFTplan);
}


ExtremeFFT4D::~ExtremeFFT4D() {
  if (LogLevel>2) printf("Desctructor called of ExtremeFFT4D. Terminating threads...");
  stopThreadsAndWaitForTermination();
  if (LogLevel>2) printf("successfully.\n");
  for (int I=0; I<numberOfThreads; I++) {
    destroySuperAlignedComplex(Slices2D[I]);
    delete[] ThreadControlParameterHolder[I];
  }
  delete xSSEObj;
  delete[] Slices2D;
  delete[] ThreadControlParameterHolder;
  
  for (int I=0; I<availFFTstepCount; I++) {
    delete availFFTsteps[I];
  }
  delete[] availFFTsteps;
  delete[] ThreadControlDataStructure;
  
  while (possibleFFTplanList->next!=NULL) {
    FFTplanType* dummy = possibleFFTplanList;
    while ((dummy->next)->next!=NULL) {
      dummy = dummy->next;
    }
    dismissPlanDetails(dummy->next, false);
    delete dummy->next;
    dummy->next = NULL;
  }
  dismissPlanDetails(possibleFFTplanList, false);
  delete possibleFFTplanList;
  possibleFFTplanList=NULL;

  delete[] PlanExecutionControlDataStructure.SliceStatus;
  delete[] PlanExecutionControlDataStructure.WaitLoopEntered;
  
  for (int I=0; I<4; I++) {
    delete[] primeFactorsL[I];
  }
  delete[] CoreMaskForThreads;
}


void ExtremeFFT4D::dismissPlanDetails(FFTplanType* plan, bool onlyFFTsteps) {
  if (plan == NULL) return;
  if (plan->FFTstepsA != NULL) {
    for (int I2=0; I2<numberOfThreads; I2++) {
      for (int I=0; I<plan->FFTstepCountA; I++) {
        delete plan->FFTstepsA[I2][I];
      }
      delete[] plan->FFTstepsA[I2];
    }
    delete[] plan->FFTstepsA;
    plan->FFTstepsA = NULL;
  }
  if (plan->FFTstepsB != NULL) {
    for (int I2=0; I2<numberOfThreads; I2++) {  
      for (int I=0; I<plan->FFTstepCountB; I++) {
        delete plan->FFTstepsB[I2][I];
      }
      delete[] plan->FFTstepsB[I2];
    }
    delete[] plan->FFTstepsB;
    plan->FFTstepsB = NULL;
  }
  if (!onlyFFTsteps) {
    plan->performance_TotalCycles = NaN;
    if (plan->performance_FFTstepCyclesA != NULL) {
      for (int I2=0; I2<numberOfThreads; I2++) {  
        delete[] plan->performance_FFTstepCyclesA[I2];
      }
      delete[] plan->performance_FFTstepCyclesA;
      plan->performance_FFTstepCyclesA = NULL;
    }
    if (plan->performance_FFTstepCyclesB != NULL) {
      for (int I2=0; I2<numberOfThreads; I2++) {  
        delete[] plan->performance_FFTstepCyclesB[I2];
      }
      delete[] plan->performance_FFTstepCyclesB;
      plan->performance_FFTstepCyclesB = NULL;
    }  
  }
}


void ExtremeFFT4D::elaboratePlanDetails(FFTplanType* plan) {
  dismissPlanDetails(plan, false);
  
  long int* InvertedAddressBaseL[4];
  InvertedAddressBaseL[0] = new long int[LargestL];
  InvertedAddressBaseL[1] = new long int[LargestL];
  InvertedAddressBaseL[2] = new long int[LargestL];
  InvertedAddressBaseL[3] = new long int[LargestL];
  int* orderedPrimeFactorsL[4];
  int orderedPrimeFactorCountL[4];
  for (int I=0; I<4; I++) orderedPrimeFactorCountL[I] = 0;
  orderedPrimeFactorsL[plan->LindA1] = new int[primeFactorCountL[plan->LindA1]];
  orderedPrimeFactorsL[plan->LindA2] = new int[primeFactorCountL[plan->LindA2]];
  orderedPrimeFactorsL[plan->LindB1] = new int[primeFactorCountL[plan->LindB1]];
  orderedPrimeFactorsL[plan->LindB2] = new int[primeFactorCountL[plan->LindB2]];
  
  for (int I=0; I<plan->FFTstepCountA; I++) {
    for (int I2=0; I2<availFFTsteps[plan->FFTstepNrA[I]]->getPrimeFacCount1(); I2++) {
      orderedPrimeFactorsL[plan->LindA1][orderedPrimeFactorCountL[plan->LindA1]] = availFFTsteps[plan->FFTstepNrA[I]]->getPrimeFacs1()[I2];
      orderedPrimeFactorCountL[plan->LindA1]++;
    }
    for (int I2=0; I2<availFFTsteps[plan->FFTstepNrA[I]]->getPrimeFacCount2(); I2++) {
      orderedPrimeFactorsL[plan->LindA2][orderedPrimeFactorCountL[plan->LindA2]] = availFFTsteps[plan->FFTstepNrA[I]]->getPrimeFacs2()[I2];
      orderedPrimeFactorCountL[plan->LindA2]++;
    }
  }
  for (int I=0; I<plan->FFTstepCountB; I++) {
    for (int I2=0; I2<availFFTsteps[plan->FFTstepNrB[I]]->getPrimeFacCount1(); I2++) {
      orderedPrimeFactorsL[plan->LindB1][orderedPrimeFactorCountL[plan->LindB1]] = availFFTsteps[plan->FFTstepNrB[I]]->getPrimeFacs1()[I2];
      orderedPrimeFactorCountL[plan->LindB1]++;
    }
    for (int I2=0; I2<availFFTsteps[plan->FFTstepNrB[I]]->getPrimeFacCount2(); I2++) {
      orderedPrimeFactorsL[plan->LindB2][orderedPrimeFactorCountL[plan->LindB2]] = availFFTsteps[plan->FFTstepNrB[I]]->getPrimeFacs2()[I2];
      orderedPrimeFactorCountL[plan->LindB2]++;
    }
  }
  
  for (int I=0; I<L[plan->LindA1]; I++) {
    InvertedAddressBaseL[plan->LindA1][I] = addressIncrementsInBytesL[plan->LindA1] * primeNumberBasedInverter(I, orderedPrimeFactorsL[plan->LindA1], primeFactorCountL[plan->LindA1]);
  }
  for (int I=0; I<L[plan->LindA2]; I++) {
    InvertedAddressBaseL[plan->LindA2][I] = addressIncrementsInBytesL[plan->LindA2] * primeNumberBasedInverter(I, orderedPrimeFactorsL[plan->LindA2], primeFactorCountL[plan->LindA2]);
  }
  for (int I=0; I<L[plan->LindB1]; I++) {
    InvertedAddressBaseL[plan->LindB1][I] = addressIncrementsInBytesL[plan->LindB1] * primeNumberBasedInverter(I, orderedPrimeFactorsL[plan->LindB1], primeFactorCountL[plan->LindB1]);
  }
   for (int I=0; I<L[plan->LindB2]; I++) {
    InvertedAddressBaseL[plan->LindB2][I] = addressIncrementsInBytesL[plan->LindB2] * primeNumberBasedInverter(I, orderedPrimeFactorsL[plan->LindB2], primeFactorCountL[plan->LindB2]);
  }

  plan->FFTstepsA = new ExtremeFFT4D_FFTstep**[numberOfThreads];
  for (int I3=0; I3<numberOfThreads; I3++) {
    plan->FFTstepsA[I3] = new ExtremeFFT4D_FFTstep*[plan->FFTstepCountA];
    int alreadyPerfPrimeFacProd1 = 1;
    int alreadyPerfPrimeFacProd2 = 1;
    for (int I=0; I<plan->FFTstepCountA; I++) {
      plan->FFTstepsA[I3][I] = availFFTsteps[plan->FFTstepNrA[I]]->deriveSpecificFFTstepInstance(L[plan->LindA1],L[plan->LindA2], plan->InnerLoopCountA, addressIncrementsInBytesL[plan->LindA1],addressIncrementsInBytesL[plan->LindA2], alreadyPerfPrimeFacProd1, alreadyPerfPrimeFacProd2, InvertedAddressBaseL[plan->LindA1], InvertedAddressBaseL[plan->LindA2]);
      plan->FFTstepsA[I3][I]->setInternalEmbedding(0, 0);
    
      for (int I2=0; I2<plan->FFTstepsA[I3][I]->getPrimeFacCount1(); I2++) {
        alreadyPerfPrimeFacProd1 *= plan->FFTstepsA[I3][I]->getPrimeFacs1()[I2];
      }
      for (int I2=0; I2<plan->FFTstepsA[I3][I]->getPrimeFacCount2(); I2++) {
        alreadyPerfPrimeFacProd2 *= plan->FFTstepsA[I3][I]->getPrimeFacs2()[I2];
      }
    }
  }

  plan->FFTstepsB = new ExtremeFFT4D_FFTstep**[numberOfThreads];
  for (int I3=0; I3<numberOfThreads; I3++) {
    plan->FFTstepsB[I3] = new ExtremeFFT4D_FFTstep*[plan->FFTstepCountB];
    int alreadyPerfPrimeFacProd1 = 1;
    int alreadyPerfPrimeFacProd2 = 1;
    for (int I=0; I<plan->FFTstepCountB; I++) {
      plan->FFTstepsB[I3][I] = availFFTsteps[plan->FFTstepNrB[I]]->deriveSpecificFFTstepInstance(L[plan->LindB1],L[plan->LindB2], plan->InnerLoopCountB, addressIncrementsInBytesL[plan->LindB1],addressIncrementsInBytesL[plan->LindB2], alreadyPerfPrimeFacProd1, alreadyPerfPrimeFacProd2, InvertedAddressBaseL[plan->LindB1], InvertedAddressBaseL[plan->LindB2]);
      plan->FFTstepsB[I3][I]->setInternalEmbedding(0, 0);
    
      for (int I2=0; I2<plan->FFTstepsB[I3][I]->getPrimeFacCount1(); I2++) {
        alreadyPerfPrimeFacProd1 *= plan->FFTstepsB[I3][I]->getPrimeFacs1()[I2];
      }
      for (int I2=0; I2<plan->FFTstepsB[I3][I]->getPrimeFacCount2(); I2++) {
        alreadyPerfPrimeFacProd2 *= plan->FFTstepsB[I3][I]->getPrimeFacs2()[I2];
      }
    }
  }

  for (int I=0; I<4; I++) {
    delete[] orderedPrimeFactorsL[I];
  }
  delete[] InvertedAddressBaseL[0];
  delete[] InvertedAddressBaseL[1];
  delete[] InvertedAddressBaseL[2];
  delete[] InvertedAddressBaseL[3];
}


void ExtremeFFT4D::optimizeInternalEmbeddingForPlan(FFTplanType* plan, int tuneLevel, bool enforceOpti) {
  long int sizeOfWhole4DFieldInBytes = 16*localIndexCount*L[0]*L[1]*L[2]*L[3];
  
  if (enforceOpti || (plan->InternalEmbeddingAOne<0) || (plan->InternalEmbeddingATwo<0) || (plan->inputSlice2DDistanceInBytesA<0)) {
    for (int I=0; I<plan->FFTstepCountA; I++) {
      plan->FFTstepsA[0][I]->setPrefetchCombineFacs(plan->PrefetchCombineA1, plan->PrefetchCombineA2);
    }
  
    if (tuneLevel != ExtremeFFT4D_TuneLevel_None) {
      int bestMinorMisses = -1;
      int bestMajorMisses = -1;
      int bestEmbedOne = -1;
      int bestEmbedTwo = -1;
      long int bestInputSliceDis = -1;
      int embedOneMax = ExtremeFFT4D_InternalEmbeddingOneMAX;
      if (tuneLevel < ExtremeFFT4D_TuneLevel_High) embedOneMax = 1;
      int embedTwoMax = L1CacheSizePerWayInBytes/L1CacheLineSizeInBytes;
      if (tuneLevel < ExtremeFFT4D_TuneLevel_Medium) embedTwoMax = 50;

      for (int embedOne=0; embedOne<embedOneMax; embedOne++) {
        for (int embedTwo=0; embedTwo<embedTwoMax; embedTwo++) {
          for (int I=0; I<plan->FFTstepCountA; I++) {
            plan->FFTstepsA[0][I]->setInternalEmbedding(embedOne, embedTwo);
	  }
          for (long int inputSliceDis=0; inputSliceDis<=L1CacheSizePerWayInBytes/2; inputSliceDis+=L1CacheSizePerWayInBytes/16) {
            int minorMisses = 0;
            int majorMisses = 0;
    
            for (int I=0; I<plan->FFTstepCountA; I++) {
              int majMis = bestMajorMisses;
	      int minMis = bestMinorMisses;
              plan->FFTstepsA[0][I]->calcL1ExcessMisses(sizeOfWhole4DFieldInBytes, inputSliceDis, minMis, majMis, true);
              if (minMis > 0) minorMisses += minMis;
              if (majMis > 0) majorMisses += majMis;
            }
	    if ((bestMinorMisses<0) || (bestMajorMisses<0) || ((minorMisses<bestMinorMisses) && (majorMisses<=bestMajorMisses)) || ((minorMisses<=bestMinorMisses) && (majorMisses<bestMajorMisses))) {
              bestMinorMisses = minorMisses;
              bestMajorMisses = majorMisses;
              bestEmbedOne = embedOne;
              bestEmbedTwo = embedTwo;
              bestInputSliceDis = inputSliceDis;
	    }
  	  }
        }
      }    
      plan->InternalEmbeddingAOne = bestEmbedOne;
      plan->InternalEmbeddingATwo = bestEmbedTwo;
      plan->inputSlice2DDistanceInBytesA = bestInputSliceDis;
      plan->MinorMissesA = bestMinorMisses;
      plan->MajorMissesA = bestMajorMisses;
    } else {
      plan->InternalEmbeddingAOne = 0;
      plan->InternalEmbeddingATwo = 0;
      plan->inputSlice2DDistanceInBytesA = 0;
      plan->MinorMissesA = -1;
      plan->MajorMissesA = -1;
    }
  }  
  
  if (enforceOpti || (plan->InternalEmbeddingBOne<0) || (plan->InternalEmbeddingBTwo<0) || (plan->inputSlice2DDistanceInBytesB<0)) {
    for (int I=0; I<plan->FFTstepCountB; I++) {
      plan->FFTstepsB[0][I]->setPrefetchCombineFacs(plan->PrefetchCombineB1, plan->PrefetchCombineB2);
    }

    if (tuneLevel != ExtremeFFT4D_TuneLevel_None) {
      int bestMinorMisses = -1;
      int bestMajorMisses = -1;
      int bestEmbedOne = -1;
      int bestEmbedTwo = -1;
      long int bestInputSliceDis = -1;
      int embedOneMax = ExtremeFFT4D_InternalEmbeddingOneMAX;
      if (tuneLevel < ExtremeFFT4D_TuneLevel_High) embedOneMax = 1;
      int embedTwoMax = L1CacheSizePerWayInBytes/L1CacheLineSizeInBytes;
      if (tuneLevel < ExtremeFFT4D_TuneLevel_Medium) embedTwoMax = 50;

      for (int embedOne=0; embedOne<embedOneMax; embedOne++) {
        for (int embedTwo=0; embedTwo<embedTwoMax; embedTwo++) {
          for (int I=0; I<plan->FFTstepCountB; I++) {
            plan->FFTstepsB[0][I]->setInternalEmbedding(embedOne, embedTwo);
  	  }
          for (long int inputSliceDis=0; inputSliceDis<=L1CacheSizePerWayInBytes/2; inputSliceDis+=L1CacheSizePerWayInBytes/16) {
            int minorMisses = 0;
            int majorMisses = 0;
    
            for (int I=0; I<plan->FFTstepCountB; I++) {
              int majMis = bestMajorMisses;
	      int minMis = bestMinorMisses;
              plan->FFTstepsB[0][I]->calcL1ExcessMisses(sizeOfWhole4DFieldInBytes, inputSliceDis, minMis, majMis, true);
              if (minMis > 0) minorMisses += minMis;
              if (majMis > 0) majorMisses += majMis;
            }
	    if ((bestMinorMisses<0) || (bestMajorMisses<0) || ((minorMisses<bestMinorMisses) && (majorMisses<=bestMajorMisses)) || ((minorMisses<=bestMinorMisses) && (majorMisses<bestMajorMisses))) {
              bestMinorMisses = minorMisses;
              bestMajorMisses = majorMisses;
              bestEmbedOne = embedOne;
              bestEmbedTwo = embedTwo;
              bestInputSliceDis = inputSliceDis;
	    }
  	  }
        }
      }
      plan->InternalEmbeddingBOne = bestEmbedOne;
      plan->InternalEmbeddingBTwo = bestEmbedTwo;
      plan->inputSlice2DDistanceInBytesB = bestInputSliceDis;
      plan->MinorMissesB = bestMinorMisses;
      plan->MajorMissesB = bestMajorMisses;
    } else {
      plan->InternalEmbeddingBOne = 0;
      plan->InternalEmbeddingBTwo = 0;
      plan->inputSlice2DDistanceInBytesB = 0;
      plan->MinorMissesB = -1;
      plan->MajorMissesB = -1;
    }
  }
  
  getPlanReadyForExecution(plan);
}


void ExtremeFFT4D::getPlanReadyForExecution(FFTplanType* plan) {
  for (int I2=0; I2<numberOfThreads; I2++) {
    for (int I=0; I<plan->FFTstepCountA; I++) {
      plan->FFTstepsA[I2][I]->setPrefetchCombineFacs(plan->PrefetchCombineA1, plan->PrefetchCombineA2);
    }
    for (int I=0; I<plan->FFTstepCountB; I++) {
      plan->FFTstepsB[I2][I]->setPrefetchCombineFacs(plan->PrefetchCombineB1, plan->PrefetchCombineB2);
    }
  }
  
  long int sizeOfWhole4DFieldInBytes = 16*localIndexCount*L[0]*L[1]*L[2]*L[3];
    
  for (int I2=0; I2<numberOfThreads; I2++) {
    plan->MinorMissesA = 0;
    plan->MajorMissesA = 0;
    plan->MinorMissesB = 0;
    plan->MajorMissesB = 0;
    
    for (int I=0; I<plan->FFTstepCountA; I++) {
      plan->FFTstepsA[I2][I]->setInternalEmbedding(plan->InternalEmbeddingAOne, plan->InternalEmbeddingATwo);
      if (I2==0) {
        int majMis = 0;
	int minMis = 0; 
	plan->FFTstepsA[0][I]->calcL1ExcessMisses(sizeOfWhole4DFieldInBytes, plan->inputSlice2DDistanceInBytesA, minMis, majMis, false);      
	if (minMis>0) plan->MinorMissesA += minMis;
	if (majMis>0) plan->MajorMissesA += majMis;      
      }      
    }
    for (int I=0; I<plan->FFTstepCountB; I++) {
      plan->FFTstepsB[I2][I]->setInternalEmbedding(plan->InternalEmbeddingBOne, plan->InternalEmbeddingBTwo);
      if (I2==0) {
        int majMis = 0;
	int minMis = 0; 
	plan->FFTstepsB[0][I]->calcL1ExcessMisses(sizeOfWhole4DFieldInBytes, plan->inputSlice2DDistanceInBytesB, minMis, majMis, false);      
	if (minMis>0) plan->MinorMissesB += minMis;
	if (majMis>0) plan->MajorMissesB += majMis;      
      }      
    }
  }
}


void ExtremeFFT4D::writeFFTplan(FFTplanType* fftPlan, bool writeSubsequent) {
  int nr = 0;
  printf("\n");
  while (fftPlan!=NULL) {
    if (!fftPlan->complete) break;
    if (writeSubsequent) {
      printf("----------- Plan Nr. %d ---------\n", nr);
    } else {
      printf("----------- Plan ---------\n");
    }
    nr++;
    printf("  Dimension order: (%d,%d) (%d,%d)\n",fftPlan->LindA1,fftPlan->LindA2,fftPlan->LindB1,fftPlan->LindB2);
    printf("  with Inner-Loop-Count/LoopCombineFactor (%d,%d) (%d,%d)\n", fftPlan->InnerLoopCountA, fftPlan->OuterLoopCombineFactorA2, fftPlan->InnerLoopCountB, fftPlan->OuterLoopCombineFactorB2);
    printf("  Prefetch-Combining: (%d,%d) (%d,%d)\n",fftPlan->PrefetchCombineA1,fftPlan->PrefetchCombineA2,fftPlan->PrefetchCombineB1,fftPlan->PrefetchCombineB2);
    printf("  Internal Embedding: (%ld,%d,%d) (%ld,%d,%d)\n",fftPlan->inputSlice2DDistanceInBytesA,fftPlan->InternalEmbeddingAOne,fftPlan->InternalEmbeddingATwo,fftPlan->inputSlice2DDistanceInBytesB,fftPlan->InternalEmbeddingBOne,fftPlan->InternalEmbeddingBTwo);    
    if ((fftPlan->MinorMissesA>=0) && (fftPlan->MajorMissesA>=0) && (fftPlan->MinorMissesB>=0) && (fftPlan->MajorMissesB>=0)) {
      printf("  L1-Misses: (%d,%d) (%d,%d)\n", fftPlan->MinorMissesA,fftPlan->MajorMissesA, fftPlan->MinorMissesB,fftPlan->MajorMissesB);    
    }
    printf("  Bus Blocking (RR, RW, WW): (%d,%d,%d) (%d,%d,%d)\n", fftPlan->blockSimReadReadA,fftPlan->blockSimReadWriteA,fftPlan->blockSimWriteWriteA,fftPlan->blockSimReadReadB,fftPlan->blockSimReadWriteB,fftPlan->blockSimWriteWriteB);
    printf("  FFTstep order:\n");
    for (int I=0; I<fftPlan->FFTstepCountA; I++) {
      printf("    ");
      if (availFFTsteps[fftPlan->FFTstepNrA[I]]->getReadFromInput()) printf("READ + ");
      printf("(");
      for (int I2=0; I2<availFFTsteps[fftPlan->FFTstepNrA[I]]->getPrimeFacCount1(); I2++) {
        printf("%d", availFFTsteps[fftPlan->FFTstepNrA[I]]->getPrimeFacs1()[I2]);
	if (I2<availFFTsteps[fftPlan->FFTstepNrA[I]]->getPrimeFacCount1()-1) printf(", ");
      }
      printf(") / (");
      for (int I2=0; I2<availFFTsteps[fftPlan->FFTstepNrA[I]]->getPrimeFacCount2(); I2++) {
        printf("%d", availFFTsteps[fftPlan->FFTstepNrA[I]]->getPrimeFacs2()[I2]);
	if (I2<availFFTsteps[fftPlan->FFTstepNrA[I]]->getPrimeFacCount2()-1) printf(", ");
      }
      printf(")");
      if (availFFTsteps[fftPlan->FFTstepNrA[I]]->getWriteToOutput()) printf(" + WRITE");
      printf("\n");
    }

    for (int I=0; I<fftPlan->FFTstepCountB; I++) {
      printf("    ");
      if (availFFTsteps[fftPlan->FFTstepNrB[I]]->getReadFromInput()) printf("READ + ");
      printf("(");
      for (int I2=0; I2<availFFTsteps[fftPlan->FFTstepNrB[I]]->getPrimeFacCount1(); I2++) {
        printf("%d", availFFTsteps[fftPlan->FFTstepNrB[I]]->getPrimeFacs1()[I2]);
	if (I2<availFFTsteps[fftPlan->FFTstepNrB[I]]->getPrimeFacCount1()-1) printf(", ");
      }
      printf(") / (");
      for (int I2=0; I2<availFFTsteps[fftPlan->FFTstepNrB[I]]->getPrimeFacCount2(); I2++) {
        printf("%d", availFFTsteps[fftPlan->FFTstepNrB[I]]->getPrimeFacs2()[I2]);
	if (I2<availFFTsteps[fftPlan->FFTstepNrB[I]]->getPrimeFacCount2()-1) printf(", ");
      }
      printf(")");
      if (availFFTsteps[fftPlan->FFTstepNrB[I]]->getWriteToOutput()) printf(" + WRITE");
      printf("\n");
    }

    if (!isNaN(fftPlan->performance_TotalCycles)) {
      printf("  Performance: Cycles per FFT: %1.2f Mega Cycles\n", fftPlan->performance_TotalCycles/1E6);
      for (int I2=0; I2<numberOfThreads; I2++) {
        printf("    Thread %d [(Perf, Admin)]: ", I2);
        for (int I=0; I<fftPlan->FFTstepCountA-1; I++) printf("(%1.2f, %1.2f), ", fftPlan->performance_FFTstepCyclesA[I2][2*I+1]/1E6, fftPlan->performance_FFTstepCyclesA[I2][2*I+0]/1E6);
        printf("(%1.2f, %1.2f) ", fftPlan->performance_FFTstepCyclesA[I2][2*fftPlan->FFTstepCountA-1]/1E6, fftPlan->performance_FFTstepCyclesA[I2][2*fftPlan->FFTstepCountA-2]/1E6);
	printf(" + %1.2f / ", fftPlan->performance_FFTstepCyclesA[I2][2*fftPlan->FFTstepCountA]/1E6);
        for (int I=0; I<fftPlan->FFTstepCountB-1; I++) printf("(%1.2f, %1.2f), ", fftPlan->performance_FFTstepCyclesB[I2][2*I+1]/1E6, fftPlan->performance_FFTstepCyclesB[I2][2*I+0]/1E6); 
	printf("(%1.2f, %1.2f) ", fftPlan->performance_FFTstepCyclesB[I2][2*fftPlan->FFTstepCountB-1]/1E6, fftPlan->performance_FFTstepCyclesB[I2][2*fftPlan->FFTstepCountB-2]/1E6);
	printf(" + %1.2f\n", fftPlan->performance_FFTstepCyclesB[I2][2*fftPlan->FFTstepCountB]/1E6);
      }
    }
   
    printf("\n");
    if (!writeSubsequent) break;
    fftPlan = fftPlan->next;
  }
}


void ExtremeFFT4D::generateAllPossibleFFTplans(FFTplanType* &FFTplan, int state) {
  //Check whether plan is reasonable / complete
  int primeFacsA1[ExtremeFFT4D_PrimeFactorsMAX+10];
  int primeFacCountA1 = 0;
  int primeFacsA2[ExtremeFFT4D_PrimeFactorsMAX+10];
  int primeFacCountA2 = 0;
  int primeFacsB1[ExtremeFFT4D_PrimeFactorsMAX+10];
  int primeFacCountB1 = 0;
  int primeFacsB2[ExtremeFFT4D_PrimeFactorsMAX+10];
  int primeFacCountB2 = 0;
  for (int I=0; I<FFTplan->FFTstepCountA; I++) {
    for (int I2=0; I2<availFFTsteps[FFTplan->FFTstepNrA[I]]->getPrimeFacCount1(); I2++) {
      primeFacsA1[primeFacCountA1] = availFFTsteps[FFTplan->FFTstepNrA[I]]->getPrimeFacs1()[I2];
      primeFacCountA1++;
    }
    for (int I2=0; I2<availFFTsteps[FFTplan->FFTstepNrA[I]]->getPrimeFacCount2(); I2++) {
      primeFacsA2[primeFacCountA2] = availFFTsteps[FFTplan->FFTstepNrA[I]]->getPrimeFacs2()[I2];
      primeFacCountA2++;
    }
  }
  for (int I=0; I<FFTplan->FFTstepCountB; I++) {
    for (int I2=0; I2<availFFTsteps[FFTplan->FFTstepNrB[I]]->getPrimeFacCount1(); I2++) {
      primeFacsB1[primeFacCountB1] = availFFTsteps[FFTplan->FFTstepNrB[I]]->getPrimeFacs1()[I2];
      primeFacCountB1++;
    }
    for (int I2=0; I2<availFFTsteps[FFTplan->FFTstepNrB[I]]->getPrimeFacCount2(); I2++) {
      primeFacsB2[primeFacCountB2] = availFFTsteps[FFTplan->FFTstepNrB[I]]->getPrimeFacs2()[I2];
      primeFacCountB2++;
    }
  }
  //1. Check whether there are too many transforms in one direction
  if (primeFacCountA1>primeFactorCountL[FFTplan->LindA1]) return;
  if (primeFacCountA2>primeFactorCountL[FFTplan->LindA2]) return;
  if (primeFacCountB1>primeFactorCountL[FFTplan->LindB1]) return;
  if (primeFacCountB2>primeFactorCountL[FFTplan->LindB2]) return;
  
  //Check completed if state indicates completion
  if (state>=2) {
    if (primeFacCountA1!=primeFactorCountL[FFTplan->LindA1]) return;
    if (primeFacCountA2!=primeFactorCountL[FFTplan->LindA2]) return;
    
    for (int I=0; I<primeFacCountA1; I++) {
      bool found = false;
      for (int I2=0; I2<primeFactorCountL[FFTplan->LindA1]; I2++) {
        if (primeFacsA1[I]==primeFactorsL[FFTplan->LindA1][I2]) {
	  primeFacsA1[I] = -1;
	  found = true;
	}
      }
      if (!found) return;
    }
    for (int I=0; I<primeFacCountA2; I++) {
      bool found = false;
      for (int I2=0; I2<primeFactorCountL[FFTplan->LindA2]; I2++) {
        if (primeFacsA2[I]==primeFactorsL[FFTplan->LindA2][I2]) {
	  primeFacsA2[I] = -1;
	  found = true;
	}
      }
      if (!found) return;
    }
  }
  if (state>=4) {
    if (primeFacCountB1!=primeFactorCountL[FFTplan->LindB1]) return;
    if (primeFacCountB2!=primeFactorCountL[FFTplan->LindB2]) return;
    
    for (int I=0; I<primeFacCountB1; I++) {
      bool found = false;
      for (int I2=0; I2<primeFactorCountL[FFTplan->LindB1]; I2++) {
        if (primeFacsB1[I]==primeFactorsL[FFTplan->LindB1][I2]) {
	  primeFacsB1[I] = -1;
	  found = true;
	}
      }
      if (!found) return;
    }
    for (int I=0; I<primeFacCountB2; I++) {
      bool found = false;
      for (int I2=0; I2<primeFactorCountL[FFTplan->LindB2]; I2++) {
        if (primeFacsB2[I]==primeFactorsL[FFTplan->LindB2][I2]) {
	  primeFacsB2[I] = -1;
	  found = true;
	}
      }
      if (!found) return;
    }
  }
  
  //Save Plan if state idicates completion
  if (state==4) {
    int locIndFacCombine=1;
    while ((locIndFacCombine<=L[3]) && (locIndFacCombine*localIndexCount<=ExtremeFFT4D_InnerLoopMAX)) {
      if ((L[3]%locIndFacCombine)==0) {
	FFTplan->next = new FFTplanType;		
        *(FFTplan->next) = *FFTplan;
	(FFTplan->next)->next = NULL;

        FFTplan->complete = true;
        FFTplan->OuterLoopCombineFactorA2 = 1;
	FFTplan->OuterLoopCombineFactorB2 = 1;
	if (FFTplan->LindA2==3) FFTplan->OuterLoopCombineFactorB2 = locIndFacCombine;
	if (FFTplan->LindB2==3) FFTplan->OuterLoopCombineFactorB2 = locIndFacCombine;
	FFTplan->InnerLoopCountA = localIndexCount * FFTplan->OuterLoopCombineFactorA2;
	FFTplan->InnerLoopCountB = localIndexCount * FFTplan->OuterLoopCombineFactorB2;
        FFTplan->inputSlice2DDistanceInBytesA = 0;
        FFTplan->inputSlice2DDistanceInBytesB = 0;  
	FFTplan->MinorMissesA = -1;
	FFTplan->MajorMissesA = -1;
	FFTplan->MinorMissesB = -1;
	FFTplan->MajorMissesB = -1;
	FFTplan->InternalEmbeddingAOne = 0;
	FFTplan->InternalEmbeddingATwo = 0;
	FFTplan->InternalEmbeddingBOne = 0;
	FFTplan->InternalEmbeddingBTwo = 0;
	FFTplan->PrefetchCombineA1 = 1;
	FFTplan->PrefetchCombineA2 = 1;
	FFTplan->PrefetchCombineB1 = 1;
	FFTplan->PrefetchCombineB2 = 1;
	FFTplan->blockSimReadReadA = false;
	FFTplan->blockSimReadWriteA = false;
	FFTplan->blockSimWriteWriteA = false;
	FFTplan->blockSimReadReadB = false;
	FFTplan->blockSimReadWriteB = false;
	FFTplan->blockSimWriteWriteB = false;

	FFTplan = FFTplan->next;
      }
    
      locIndFacCombine++;
      if (locIndFacCombine*LargestL*LargestL*localIndexCount*16 > L2CacheSizeInBytes/2) break;
    }
    return;
  }
  
  
  for (int I=0; I<availFFTstepCount; I++) {
    bool eligible = true;
    
    //Check MustBeFirstTrafo
    if (state<2) {
      for (int I2=0; I2<FFTplan->FFTstepCountA; I2++) {
        if ((availFFTsteps[FFTplan->FFTstepNrA[I2]]->getPrimeFacCount1()>0) && (availFFTsteps[I]->getMustBeFirstTrafo1())) {
	  eligible=false;
	}
        if ((availFFTsteps[FFTplan->FFTstepNrA[I2]]->getPrimeFacCount2()>0) && (availFFTsteps[I]->getMustBeFirstTrafo2())) {
	  eligible=false;
	}
      }
    } else {
      for (int I2=0; I2<FFTplan->FFTstepCountB; I2++) {
        if ((availFFTsteps[FFTplan->FFTstepNrB[I2]]->getPrimeFacCount1()>0) && (availFFTsteps[I]->getMustBeFirstTrafo1())) {
	  eligible=false;
	}
        if ((availFFTsteps[FFTplan->FFTstepNrB[I2]]->getPrimeFacCount2()>0) && (availFFTsteps[I]->getMustBeFirstTrafo2())) {
	  eligible=false;
	}
      }
    }
    
    if (((state%2)==0) && (!availFFTsteps[I]->getReadFromInput())) eligible = false;
    if (((state%2)==1) && (availFFTsteps[I]->getReadFromInput())) eligible = false;
    
    if (eligible) {
      int stateNEW = state;
      if (((state%2)==0) && (availFFTsteps[I]->getReadFromInput())) {  //Read from input
        stateNEW++;
      }
      if (((state%2)==1) && (!availFFTsteps[I]->getReadFromInput())) {  //Read from slice
	if (availFFTsteps[I]->getWriteToOutput()) stateNEW++;
      }
  
      //Recursive call
      if (state<2) {
	FFTplan->FFTstepNrA[FFTplan->FFTstepCountA] = I;
        FFTplan->FFTstepCountA++;
        generateAllPossibleFFTplans(FFTplan, stateNEW);
        FFTplan->FFTstepCountA--;
      } else {
	FFTplan->FFTstepNrB[FFTplan->FFTstepCountB] = I;
        FFTplan->FFTstepCountB++;
        generateAllPossibleFFTplans(FFTplan, stateNEW);
        FFTplan->FFTstepCountB--;
      }
    }
  }
}


void ExtremeFFT4D::tune(int tuneLevel) {
  int VLxtrSize = localIndexCount*L[0]*(L[1]+xtrSize1)*(L[2]+xtrSize2)*(L[3]+xtrSize3);
  Complex* input = createSuperAlignedComplex(VLxtrSize);
  Complex* output = createSuperAlignedComplex(VLxtrSize);

  tune(tuneLevel, input, output);

  destroySuperAlignedComplex(input);
  destroySuperAlignedComplex(output);
}


void ExtremeFFT4D::tune(int tuneLevel, Complex* input, Complex* output) {
  if (LogLevel>2) printf("Tuning FFT...\n");
  int VLxtrSize = localIndexCount*L[0]*(L[1]+xtrSize1)*(L[2]+xtrSize2)*(L[3]+xtrSize3);  
  for (int I=0; I<VLxtrSize; I++) {
    input[I].x = 0;
    input[I].y = 0;
    output[I].x = NaN;
    output[I].y = NaN;    
  }
  int n[4];
  n[0] = 1 % L[0];
  n[1] = 2 % L[1];
  n[2] = 3 % L[2];
  n[3] = 4 % L[3];  
  int deltaPos = localIndexCount*(n[3] + n[2]*(L[3]+xtrSize3) + n[1]*(L[3]+xtrSize3)*(L[2]+xtrSize2) + n[0]*(L[3]+xtrSize3)*(L[2]+xtrSize2)*(L[1]+xtrSize1));
  int volume = L[0]*L[1]*L[2]*L[3];
  int p = 0;
  for (int I0=0; I0<L[0]; I0++) {
    double angle0 = (n[0]*I0*2*pi) / L[0];
    for (int I1=0; I1<L[1]; I1++) {
      double angle1 = (n[1]*I1*2*pi) / L[1];
      for (int I2=0; I2<L[2]; I2++) {
        double angle2 = (n[2]*I2*2*pi) / L[2];
        for (int I3=0; I3<L[3]; I3++) {
          double angle3 = (n[3]*I3*2*pi) / L[3];
	  double c = cos(angle0+angle1+angle2+angle3);
	  double s = sin(angle0+angle1+angle2+angle3);

	  for (int i=0; i<localIndexCount; i++) {
	    input[p].x = c;
	    input[p].y = s;
	    p++;
	  }
	}
	p += localIndexCount*xtrSize3;
      }
      p += localIndexCount*xtrSize2*(L[3]+xtrSize3);
    }  
    p += localIndexCount*xtrSize1*(L[2]+xtrSize2)*(L[3]+xtrSize3);
  }
  
  int CopyCount = 0;  
  for (int rCA1=1; rCA1<=ExtremeFFT4D_PrefetchCombine1MAX; rCA1*=2) {
    for (int rCA2=1; rCA2<=ExtremeFFT4D_PrefetchCombine2MAX; rCA2*=2) {
      CopyCount++;
    }
  }
  for (int rCB1=1; rCB1<=ExtremeFFT4D_PrefetchCombine1MAX; rCB1*=2) {
    for (int rCB2=1; rCB2<=ExtremeFFT4D_PrefetchCombine2MAX; rCB2*=2) {
      CopyCount++;
    }
  }
  CopyCount*=8;
  
  //Set all plans to uninitialized status
  FFTplanType* plan = possibleFFTplanList;
  while (plan!=NULL) {
    dismissPlanDetails(plan, false);
    plan->PrefetchCombineA1 = 1;
    plan->PrefetchCombineA2 = 1;
    plan->PrefetchCombineB1 = 1;
    plan->PrefetchCombineB2 = 1;
    plan->inputSlice2DDistanceInBytesA = -1;
    plan->inputSlice2DDistanceInBytesB = -1; 
    plan->MinorMissesA = -1;
    plan->MajorMissesA = -1;
    plan->MinorMissesB = -1;
    plan->MajorMissesB = -1;
    plan->InternalEmbeddingAOne = -1;
    plan->InternalEmbeddingATwo = -1;
    plan->InternalEmbeddingBOne = -1;
    plan->InternalEmbeddingBTwo = -1;
    plan->blockSimReadReadA = false;
    plan->blockSimReadWriteA = false;
    plan->blockSimWriteWriteA = false;
    plan->blockSimReadReadB = false;
    plan->blockSimReadWriteB = false;
    plan->blockSimWriteWriteB = false;    
    
    plan = plan->next;
  }
  
  plan = possibleFFTplanList;
  FFTplanType** planCopies = new FFTplanType*[CopyCount];
  for (int I=0; I<CopyCount; I++) planCopies[I] = new FFTplanType;
  while (plan!=NULL) {
    if (!plan->complete) break;
    elaboratePlanDetails(plan);
    
    int count = 0;
    if (plan->performance_FFTstepCyclesA == NULL) {
      for (int rCA1=1; rCA1<=ExtremeFFT4D_PrefetchCombine1MAX; rCA1*=2) {
        for (int rCA2=1; rCA2<=ExtremeFFT4D_PrefetchCombine2MAX; rCA2*=2) {
	  plan->PrefetchCombineA1 = rCA1;
	  plan->PrefetchCombineA2 = rCA2;
	  plan->inputSlice2DDistanceInBytesA = -1;
	  plan->MinorMissesA = -1;
	  plan->MajorMissesA = -1;
	  plan->InternalEmbeddingAOne = -1;
	  plan->InternalEmbeddingATwo = -1;
          optimizeInternalEmbeddingForPlan(plan, tuneLevel, false);
	  
	  for (int blockRR=0; blockRR<2; blockRR++) {
  	    for (int blockRW=0; blockRW<2; blockRW++) {
	      for (int blockWW=0; blockWW<2; blockWW++) {
                *(planCopies[count]) = *plan;
                planCopies[count]->next = NULL;

		planCopies[count]->blockSimReadReadA = blockRR;
		planCopies[count]->blockSimReadWriteA = blockRW;
		planCopies[count]->blockSimWriteWriteA = blockWW;
		
		count++;
	      }
	    }
	  }
	}
      }
    }

    if (plan->performance_FFTstepCyclesB == NULL) {
      for (int rCB1=1; rCB1<=ExtremeFFT4D_PrefetchCombine1MAX; rCB1*=2) {
        for (int rCB2=1; rCB2<=ExtremeFFT4D_PrefetchCombine2MAX; rCB2*=2) {
	  plan->PrefetchCombineB1 = rCB1;
	  plan->PrefetchCombineB2 = rCB2;
	  plan->inputSlice2DDistanceInBytesB = -1;  
	  plan->MinorMissesB = -1;
	  plan->MajorMissesB = -1;
	  plan->InternalEmbeddingBOne = -1;
	  plan->InternalEmbeddingBTwo = -1;
          optimizeInternalEmbeddingForPlan(plan, tuneLevel, false);
	
	  for (int blockRR=0; blockRR<2; blockRR++) {
  	    for (int blockRW=0; blockRW<2; blockRW++) {
	      for (int blockWW=0; blockWW<2; blockWW++) {
                *(planCopies[count]) = *plan;
                planCopies[count]->next = NULL;

		planCopies[count]->blockSimReadReadB = blockRR;
		planCopies[count]->blockSimReadWriteB = blockRW;
		planCopies[count]->blockSimWriteWriteB = blockWW;

		count++;
	      }
	    }
	  }
	}
      }
    }

    double bestPerformanceA = NaN;
    double bestPerformanceB = NaN;    
    int bestPlanA = -1;      
    int bestPlanB = -1;      
    for (int I=0; I<count; I++) {
      planCopies[I]->performance_TotalCycles = NaN;     
      planCopies[I]->performance_FFTstepCyclesA = NULL;
      planCopies[I]->performance_FFTstepCyclesB = NULL;

      if (LogLevel>2) printf("Testing FFT...");
      getPlanReadyForExecution(planCopies[I]);
      for (int i=0; i<VLxtrSize; i++) {
        output[i].x = NaN;
        output[i].y = NaN;
      }      
      executePlan(input, output, planCopies[I], 2, ExtremeFFT4D_Forward);    
      if (LogLevel>2) printf(" takes %1.2f Mega Cycles\n", planCopies[I]->performance_TotalCycles/1E6);
      
      double perfCycA = 0;
      double perfCycB = 0;      
      for (int i=0; i<1+2*planCopies[I]->FFTstepCountA; i++) {
        perfCycA += planCopies[I]->performance_FFTstepCyclesA[0][i];
      }
      for (int i=0; i<1+2*planCopies[I]->FFTstepCountB; i++) {
        perfCycB += planCopies[I]->performance_FFTstepCyclesB[0][i];
      }

      if (isNaN(bestPerformanceA) || (bestPerformanceA>perfCycA)) {
        bestPerformanceA = perfCycA;
	bestPlanA = I;
      }
      if (isNaN(bestPerformanceB) || (bestPerformanceB>perfCycB)) {
        bestPerformanceB = perfCycB;
	bestPlanB = I;
      }
	
      //Check correctness of output
      for (int i=0; i<localIndexCount; i++) {
        if (fabs(output[deltaPos+i].x-volume)>1E-15*volume) {
	  printf("ERROR in ExtremeFFT4D: FFT-malfunction (1: val=%e, local-index=%d)!!!\n", output[deltaPos+i].x,i);
	  writeFFTplan(planCopies[I], false);
	  exit(0);
	}
	output[deltaPos+i].x=0;
	output[deltaPos+i].y=0;
      }

      int pos = 0;
      for (int I0=0; I0<L[0]; I0++) {
        for (int I1=0; I1<L[1]; I1++) {
          for (int I2=0; I2<L[2]; I2++) {
            for (int I3=0; I3<L[3]; I3++) {
  	      for (int i=0; i<localIndexCount; i++) {
                if ((isNaN(output[pos].x)) || (isNaN(output[pos].y)) ||(fabs(output[pos].x)>1E-15*volume) || (fabs(output[pos].y)>1E-15*volume)) {	
	          printf("ERROR in ExtremeFFT4D: FFT-malfunction (2: val=%e, %e, index=%d)!!!\n", output[pos].x, output[pos].y,pos);
	          writeFFTplan(planCopies[I], false);
	          exit(0);
 	        }
  	        pos++;
	      }
  	    }
   	    pos += localIndexCount*xtrSize3;
          }
          pos += localIndexCount*xtrSize2*(L[3]+xtrSize3);
        }  
        pos += localIndexCount*xtrSize1*(L[2]+xtrSize2)*(L[3]+xtrSize3);
      }
    }


    //Take over results for best Trafo A to equivalent transformations
    if (!isNaN(bestPerformanceA)) {
      FFTplanType* plan2 = plan;
      while (plan2 != NULL) {
        if ((plan2->complete) && (plan2->performance_FFTstepCyclesA==NULL)) {
	  if ((plan->LindA1==plan2->LindA1) && (plan->LindA2==plan2->LindA2) && (plan->InnerLoopCountA==plan2->InnerLoopCountA)
  	   && (plan->OuterLoopCombineFactorA2==plan2->OuterLoopCombineFactorA2) && (plan->FFTstepCountA==plan2->FFTstepCountA)) {
            bool fftStepsIdentical = true;
	    for (int i=0; i<plan->FFTstepCountA; i++) {
	      if (plan->FFTstepNrA[i]!=plan2->FFTstepNrA[i]) fftStepsIdentical = false;
            }	
	    if (fftStepsIdentical) {
	      plan2->PrefetchCombineA1 = planCopies[bestPlanA]->PrefetchCombineA1;
	      plan2->PrefetchCombineA2 = planCopies[bestPlanA]->PrefetchCombineA2;
	      plan2->inputSlice2DDistanceInBytesA = planCopies[bestPlanA]->inputSlice2DDistanceInBytesA;
  	      plan2->MinorMissesA = planCopies[bestPlanA]->MinorMissesA;
	      plan2->MajorMissesA = planCopies[bestPlanA]->MajorMissesA;
	      plan2->InternalEmbeddingAOne = planCopies[bestPlanA]->InternalEmbeddingAOne;
	      plan2->InternalEmbeddingATwo = planCopies[bestPlanA]->InternalEmbeddingATwo;	      
	      plan2->blockSimReadReadA = planCopies[bestPlanA]->blockSimReadReadA;
	      plan2->blockSimReadWriteA = planCopies[bestPlanA]->blockSimReadWriteA;
	      plan2->blockSimWriteWriteA = planCopies[bestPlanA]->blockSimWriteWriteA;
	      
	      plan2->performance_FFTstepCyclesA = new double*[numberOfThreads];
	      for (int i=0; i<numberOfThreads; i++) {
	        plan2->performance_FFTstepCyclesA[i] = new double[1+2*plan->FFTstepCountA];
		for (int i2=0; i2<1+2*plan->FFTstepCountA; i2++) {
		  plan2->performance_FFTstepCyclesA[i][i2] = planCopies[bestPlanA]->performance_FFTstepCyclesA[i][i2];
		}
	      }
	    }
	  }
	}
	plan2=plan2->next;
      }
    }
    
    //Take over results for best Trafo B to equivalent transformations
    if (!isNaN(bestPerformanceB)) {
      FFTplanType* plan2 = plan;
      while (plan2 != NULL) {
        if ((plan2->complete) && (plan2->performance_FFTstepCyclesB==NULL)) {
	  if ((plan->LindB1==plan2->LindB1) && (plan->LindB2==plan2->LindB2) && (plan->InnerLoopCountB==plan2->InnerLoopCountB)
  	   && (plan->OuterLoopCombineFactorB2==plan2->OuterLoopCombineFactorB2) && (plan->FFTstepCountB==plan2->FFTstepCountB)) {
            bool fftStepsIdentical = true;
	    for (int i=0; i<plan->FFTstepCountB; i++) {
	      if (plan->FFTstepNrB[i]!=plan2->FFTstepNrB[i]) fftStepsIdentical = false;
            }	
	    if (fftStepsIdentical) {
	      plan2->PrefetchCombineB1 = planCopies[bestPlanB]->PrefetchCombineB1;
	      plan2->PrefetchCombineB2 = planCopies[bestPlanB]->PrefetchCombineB2;
	      plan2->inputSlice2DDistanceInBytesB = planCopies[bestPlanB]->inputSlice2DDistanceInBytesB;
  	      plan2->MinorMissesB = planCopies[bestPlanB]->MinorMissesB;
	      plan2->MajorMissesB = planCopies[bestPlanB]->MajorMissesB;
	      plan2->InternalEmbeddingBOne = planCopies[bestPlanB]->InternalEmbeddingBOne;
	      plan2->InternalEmbeddingBTwo = planCopies[bestPlanB]->InternalEmbeddingBTwo;	      
	      plan2->blockSimReadReadB = planCopies[bestPlanB]->blockSimReadReadB;
	      plan2->blockSimReadWriteB = planCopies[bestPlanB]->blockSimReadWriteB;
	      plan2->blockSimWriteWriteB = planCopies[bestPlanB]->blockSimWriteWriteB;
	      
	      plan2->performance_FFTstepCyclesB = new double*[numberOfThreads];
	      for (int i=0; i<numberOfThreads; i++) {
	        plan2->performance_FFTstepCyclesB[i] = new double[1+2*plan->FFTstepCountB];
		for (int i2=0; i2<1+2*plan->FFTstepCountB; i2++) {
		  plan2->performance_FFTstepCyclesB[i][i2] = planCopies[bestPlanB]->performance_FFTstepCyclesB[i][i2];
		}
	      }
	    }
	  }
	}
	plan2=plan2->next;
      }
    }

    for (int I=0; I<count; I++) {
      for (int i=0; i<numberOfThreads; i++) {
        delete[] planCopies[I]->performance_FFTstepCyclesA[i];
        delete[] planCopies[I]->performance_FFTstepCyclesB[i];	
      }
      delete[] planCopies[I]->performance_FFTstepCyclesA;      
      delete[] planCopies[I]->performance_FFTstepCyclesB;      
    }
    
    dismissPlanDetails(plan, true);
    plan = plan->next;
  }
  for (int I=0; I<CopyCount; I++) delete planCopies[I];
  delete[] planCopies;


  //Calc Total-Cycles
  plan = possibleFFTplanList;
  double bestPerformance = NaN;
  selectedFFTplan = NULL;
  while (plan!=NULL) {
    if (!plan->complete) break;
    double perfCyc = 0;
    for (int i=0; i<1+2*plan->FFTstepCountA; i++) {
      perfCyc += plan->performance_FFTstepCyclesA[0][i];
    }
    for (int i=0; i<1+2*plan->FFTstepCountB; i++) {
      perfCyc += plan->performance_FFTstepCyclesB[0][i];
    }
    plan->performance_TotalCycles = perfCyc;
    if (isNaN(bestPerformance) || (plan->performance_TotalCycles<bestPerformance)) {
      bestPerformance = plan->performance_TotalCycles;
      selectedFFTplan = plan;
    }
    
    plan = plan->next;
  }
  
  elaboratePlanDetails(selectedFFTplan);
  optimizeInternalEmbeddingForPlan(selectedFFTplan, tuneLevel, false);
  getPlanReadyForExecution(selectedFFTplan); 
  executePlan(input, output, selectedFFTplan, 2, ExtremeFFT4D_Forward);    

  if (LogLevel>2) {
    printf("\nSelected Plan after tuning: \n");
    writeFFTplan(selectedFFTplan, false);
  }
}


void ExtremeFFT4D::FastFourierTrafo(Complex* input, Complex* output, bool fft_forward) {
  executePlan(input, output, selectedFFTplan, 1, fft_forward);
//writeFFTplan(selectedFFTplan, false);
}


void ExtremeFFT4D::threadedExecutionOfPlan(int ThreadID, xSSE* xSSEObjLoc) {
  Complex* input = PlanExecutionControlDataStructure.input;
  Complex* output = PlanExecutionControlDataStructure.output;
  FFTplanType* plan = PlanExecutionControlDataStructure.plan;
  bool forward = PlanExecutionControlDataStructure.forward;
  bool doTiming = PlanExecutionControlDataStructure.doTiming;
  int counter;
  long int* WaitLoopEnteredAddr = &(PlanExecutionControlDataStructure.WaitLoopEntered[ThreadID]);
  long int WaitLoopEntered;
  long int* ReadFlagAddr = &(PlanExecutionControlDataStructure.ReadFlag);
  long int* WriteFlagAddr = &(PlanExecutionControlDataStructure.WriteFlag);
  long int ReadFlag = 0;
  long int ReadFlag2 = 0;
  long int WriteFlag = 0;
  long int WriteFlag2 = 0;
  long int* readyTrafoAddr = &(PlanExecutionControlDataStructure.readyTrafoA);
  long int readyTrafo = 0;
  int OuterLoopInd1 = plan->LindB1;
  int OuterLoopInd2 = plan->LindB2;
  int OuterLoopCombineFactor = plan->OuterLoopCombineFactorA2;
  int localFFTstepCount = plan->FFTstepCountA;
  ExtremeFFT4D_FFTstep** localFFTsteps = plan->FFTstepsA[ThreadID];
  double* localPerformance_FFTstepCycles = plan->performance_FFTstepCyclesA[ThreadID];
  double lastCPUCycleCount = NaN;
  if (doTiming) {
    lastCPUCycleCount = getCPUCycleCounter();
  }
  bool blockSimReadRead = plan->blockSimReadReadA;
  bool blockSimReadWrite = plan->blockSimReadWriteA;
  bool blockSimWriteWrite = plan->blockSimWriteWriteA;
  long int inputSlice2DDistanceInBytes = plan->inputSlice2DDistanceInBytesA;
    
 
  for (int TrafoType=1; TrafoType<=2; TrafoType++) {
    for (int ModuloSelector=1; ModuloSelector>=1; ModuloSelector--) {
      counter = 0;  
      
      for (int I1=0; I1<L[OuterLoopInd1]; I1++) {
        for (int I2=0; I2<L[OuterLoopInd2]; I2+=OuterLoopCombineFactor) {
          if ((counter % ((ModuloSelector*(numberOfThreads-1))+1)) == (ModuloSelector*ThreadID)) {
            long int* sliceStatusAddr = &(PlanExecutionControlDataStructure.SliceStatus[counter]);
	    long int sliceStatus = 0;
    	    xSSEObjLoc->xSSE_ReadLongIntFromMemAddr_Wrapper(sliceStatusAddr, sliceStatus);
    	    if (sliceStatus == 0) {
              sliceStatus = ThreadID+1;
	      xSSEObjLoc->xSSE_WriteLongIntToMemAddr_Wrapper(sliceStatusAddr, sliceStatus);
  	      xSSEObjLoc->xSSE_ReadLongIntFromMemAddr_Wrapper(sliceStatusAddr, sliceStatus);
       	      xSSEObjLoc->xSSE_ReadLongIntFromMemAddr_Wrapper(sliceStatusAddr, sliceStatus);
	      if (sliceStatus == ThreadID+1) {
                int index = (I1*addressIncrementsInBytesL[OuterLoopInd1] + I2*addressIncrementsInBytesL[OuterLoopInd2]) / 16;
                Complex* pIn = &(input[index]);
                Complex* pOut = &(output[index]);
		int slice2DIndex = (inputSlice2DDistanceInBytes - (((long int) Slices2D[ThreadID]) - ((long int) pIn))) % L1CacheSizePerWayInBytes;
		if (slice2DIndex<0) slice2DIndex = (L1CacheSizePerWayInBytes + slice2DIndex) % L1CacheSizePerWayInBytes;
		slice2DIndex /= 16;
 	        for (int I=0; I<localFFTstepCount; I++) {
	          bool readFromInput = localFFTsteps[I]->getReadFromInput();
	          bool writeToOutput = localFFTsteps[I]->getWriteToOutput();		  
		  
		  while (true) {		  
		    xSSEObjLoc->xSSE_ReadLongIntFromMemAddr_Wrapper(ReadFlagAddr, ReadFlag);
		    xSSEObjLoc->xSSE_ReadLongIntFromMemAddr_Wrapper(WriteFlagAddr, WriteFlag);
 		    if ((!blockSimReadRead) || (!readFromInput) || (ReadFlag == 0)) {
   		      if ((!blockSimWriteWrite) || (!writeToOutput) || (WriteFlag == 0)) {
 		        if ((!blockSimReadWrite) || (((!readFromInput) || (WriteFlag == 0)) && ((!writeToOutput) || (ReadFlag == 0)))) {			
		          if (readFromInput) {
		  	    ReadFlag = ThreadID+1;
  		            xSSEObjLoc->xSSE_WriteLongIntToMemAddr_Wrapper(ReadFlagAddr, ReadFlag);			  
			  }
		          if (writeToOutput) {
			    WriteFlag = ThreadID+1;
  		            xSSEObjLoc->xSSE_WriteLongIntToMemAddr_Wrapper(WriteFlagAddr, WriteFlag);			  
		  	  }

      		          xSSEObjLoc->xSSE_ReadLongIntFromMemAddr_Wrapper(ReadFlagAddr, ReadFlag2);
		          xSSEObjLoc->xSSE_ReadLongIntFromMemAddr_Wrapper(WriteFlagAddr, WriteFlag2);
      		          xSSEObjLoc->xSSE_ReadLongIntFromMemAddr_Wrapper(ReadFlagAddr, ReadFlag2);
		          xSSEObjLoc->xSSE_ReadLongIntFromMemAddr_Wrapper(WriteFlagAddr, WriteFlag2);
		  
		          bool breakLoop = false;
		          if ((ReadFlag == ReadFlag2) && (WriteFlag == WriteFlag2)) {
	  		    if (doTiming) {
			      double currentCPUCycleCount = getCPUCycleCounter();
			      localPerformance_FFTstepCycles[2*I+0] += currentCPUCycleCount - lastCPUCycleCount;
                              lastCPUCycleCount = currentCPUCycleCount;
			    }			    			   			    
                            localFFTsteps[I]->executeFFTstep(pIn, pOut, &(Slices2D[ThreadID][slice2DIndex]), forward);			    
                            if (doTiming) {
			      double currentCPUCycleCount = getCPUCycleCounter();
			      localPerformance_FFTstepCycles[2*I+1] += currentCPUCycleCount - lastCPUCycleCount;
                              lastCPUCycleCount = currentCPUCycleCount;
			    }
			    breakLoop = true;
		  	  } 
			
  	 	          if (readFromInput) {
      		            xSSEObjLoc->xSSE_ReadLongIntFromMemAddr_Wrapper(ReadFlagAddr, ReadFlag2);
		  	    if (ReadFlag2 == 1+ThreadID) {
			      ReadFlag2 = 0;
   		              xSSEObjLoc->xSSE_WriteLongIntToMemAddr_Wrapper(ReadFlagAddr, ReadFlag2);			  
			    }
		 	  }
			  if (writeToOutput) {
      		            xSSEObjLoc->xSSE_ReadLongIntFromMemAddr_Wrapper(WriteFlagAddr, WriteFlag2);
			    if (WriteFlag2 == 1+ThreadID) {
			      WriteFlag2 = 0;
   		              xSSEObjLoc->xSSE_WriteLongIntToMemAddr_Wrapper(WriteFlagAddr, WriteFlag2);			  
			    }
		  	  }
			  if (breakLoop) break;
		        }
		      }
		    } 

		    xSSEObjLoc->xSSE_PerformNopLoop_Wrapper(ThreadID+ExtremeFFT4D_ThreadControlLoopWaitCycles);
		  } 
	        }
		
                sliceStatus = -(ThreadID+1);
	        xSSEObjLoc->xSSE_WriteLongIntToMemAddr_Wrapper(sliceStatusAddr, sliceStatus);
 	      }
            }  
          }
          counter++;
        }
      }
    }
    WaitLoopEntered = TrafoType;
    xSSEObjLoc->xSSE_WriteLongIntToMemAddr_Wrapper(WaitLoopEnteredAddr, WaitLoopEntered);
    
    while (true) {
      xSSEObjLoc->xSSE_ReadLongIntFromMemAddr_Wrapper(readyTrafoAddr, readyTrafo);
      if (readyTrafo == 1) break;
    
      if (ThreadID == 0) {
        bool allThreadsWaiting = true;
        for (int I=1; I<numberOfThreads; I++) {
	  long int* threadWaitAddr = &(PlanExecutionControlDataStructure.WaitLoopEntered[I]);
	  long int threadWait = 0;
	  xSSEObjLoc->xSSE_ReadLongIntFromMemAddr_Wrapper(threadWaitAddr, threadWait);
	  if (threadWait != TrafoType) allThreadsWaiting = false;
	}
      
        if (allThreadsWaiting) {
          bool notReady = false;
          counter = 0;  
          for (int I1=0; I1<L[OuterLoopInd1]; I1++) {
            for (int I2=0; I2<L[OuterLoopInd2]; I2+=OuterLoopCombineFactor) {
              long int* sliceStatusAddr = &(PlanExecutionControlDataStructure.SliceStatus[counter]);
	      long int sliceStatus = 0;
      	      xSSEObjLoc->xSSE_ReadLongIntFromMemAddr_Wrapper(sliceStatusAddr, sliceStatus);
	
    	      if (sliceStatus >= 0) {
	        notReady = true;
	      }
	      counter++;
  	    }
          }
          if (!notReady) {
            for (int I=0; I<LargestL*LargestL; I++) {
              long int* sliceStatusAddr = &(PlanExecutionControlDataStructure.SliceStatus[I]);
	      long int sliceStatus = 0;
    	      xSSEObjLoc->xSSE_WriteLongIntToMemAddr_Wrapper(sliceStatusAddr, sliceStatus);
  	    }
            readyTrafo = 1;
            xSSEObjLoc->xSSE_WriteLongIntToMemAddr_Wrapper(readyTrafoAddr, readyTrafo);
          }
	}
      }
   
      xSSEObjLoc->xSSE_PerformNopLoop_Wrapper(ThreadID+ExtremeFFT4D_ThreadControlLoopWaitCycles);
    }
    
    if (ThreadID == 0) {
      while(true) {
        bool allThreadsHaveLeft = true;
        for (int I=1; I<numberOfThreads; I++) {
  	  long int* threadWaitAddr = &(PlanExecutionControlDataStructure.WaitLoopEntered[I]);
          long int threadWait = 0;
  	  xSSEObjLoc->xSSE_ReadLongIntFromMemAddr_Wrapper(threadWaitAddr, threadWait);
          if (threadWait == TrafoType) allThreadsHaveLeft = false;
        }
	
	if (allThreadsHaveLeft) break;
        xSSEObjLoc->xSSE_PerformNopLoop_Wrapper(ExtremeFFT4D_ThreadControlLoopWaitCycles);	
      }
    }

    if (doTiming) {
      double currentCPUCycleCount = getCPUCycleCounter();
      localPerformance_FFTstepCycles[2*localFFTstepCount+0] += currentCPUCycleCount - lastCPUCycleCount;
      lastCPUCycleCount = currentCPUCycleCount;
    }

    readyTrafoAddr = &(PlanExecutionControlDataStructure.readyTrafoB);
    inputSlice2DDistanceInBytes = plan->inputSlice2DDistanceInBytesB;
    OuterLoopInd1 = plan->LindA1;
    OuterLoopInd2 = plan->LindA2;
    OuterLoopCombineFactor = plan->OuterLoopCombineFactorB2;
    blockSimReadRead = plan->blockSimReadReadB;
    blockSimReadWrite = plan->blockSimReadWriteB;
    blockSimWriteWrite = plan->blockSimWriteWriteB;
    localFFTstepCount = plan->FFTstepCountB;
    localFFTsteps = plan->FFTstepsB[ThreadID];
    localPerformance_FFTstepCycles = plan->performance_FFTstepCyclesB[ThreadID];
    input = output;
    WaitLoopEntered = 0;
    xSSEObjLoc->xSSE_WriteLongIntToMemAddr_Wrapper(WaitLoopEnteredAddr, WaitLoopEntered);        
  }
}


void ExtremeFFT4D::executePlan(Complex* input, Complex* output, FFTplanType* plan, int timingIterations, bool forward) {
  startThreads();
  int extraMeasurementMax = 100;

  bool doTiming = (timingIterations>0);
  if (timingIterations<=0) timingIterations = 1;
  double** bestPerformance_FFTstepCyclesA = NULL;
  double** bestPerformance_FFTstepCyclesB = NULL;

  if (doTiming) {
    if (plan->performance_FFTstepCyclesA == NULL) {
      plan->performance_FFTstepCyclesA = new double*[numberOfThreads];
      for (int I=0; I<numberOfThreads; I++) {
        plan->performance_FFTstepCyclesA[I] = new double[2*plan->FFTstepCountA+1];
      }
    }
    if (plan->performance_FFTstepCyclesB == NULL) {
      plan->performance_FFTstepCyclesB = new double*[numberOfThreads];
      for (int I=0; I<numberOfThreads; I++) {
        plan->performance_FFTstepCyclesB[I] = new double[2*plan->FFTstepCountB+1];
      }
    }
    
    bestPerformance_FFTstepCyclesA = new double*[numberOfThreads];
    bestPerformance_FFTstepCyclesB = new double*[numberOfThreads];
    for (int I=0; I<numberOfThreads; I++) {
      bestPerformance_FFTstepCyclesA[I] = new double[2*plan->FFTstepCountA+1];
      for (int I2=0; I2<2*plan->FFTstepCountA+1; I2++) {
        bestPerformance_FFTstepCyclesA[I][I2] = NaN;
	plan->performance_FFTstepCyclesA[I][I2] = 0;
      }
      bestPerformance_FFTstepCyclesB[I] = new double[2*plan->FFTstepCountB+1];
      for (int I2=0; I2<2*plan->FFTstepCountB+1; I2++) {
        bestPerformance_FFTstepCyclesB[I][I2] = NaN;
	plan->performance_FFTstepCyclesB[I][I2] = 0;
      }
    }
  }

  for (int iter=0; iter<timingIterations; iter++) {
    for (int I=0; I<LargestL*LargestL; I++) {
      PlanExecutionControlDataStructure.SliceStatus[I] = 0;
    }
    for (int I=0; I<numberOfThreads; I++) {
      PlanExecutionControlDataStructure.WaitLoopEntered[I] = 0;
    }
    PlanExecutionControlDataStructure.ReadFlag = 0;
    PlanExecutionControlDataStructure.WriteFlag = 0;
    PlanExecutionControlDataStructure.readyTrafoA = 0;
    PlanExecutionControlDataStructure.readyTrafoB = 0;
    PlanExecutionControlDataStructure.input = input;
    PlanExecutionControlDataStructure.output = output;
    PlanExecutionControlDataStructure.plan = plan;
    PlanExecutionControlDataStructure.forward = forward;
    PlanExecutionControlDataStructure.doTiming = doTiming;

    for (int I=1; I<numberOfThreads; I++) {
      long int* commandAddr = &(ThreadControlDataStructure[I].command);
      long int command = ExtremeFFT4D_ThreadControlCommand_PlanExecution;
      xSSEObj->xSSE_WriteLongIntToMemAddr_Wrapper(commandAddr, command);
    }

    threadedExecutionOfPlan(0, xSSEObj);

    if (doTiming) {
      double perfCyclesA = 0;
      double perfCyclesB = 0;
      double bestPerfCyclesA = 0;
      double bestPerfCyclesB = 0;
     
      for (int I2=0; I2<2*plan->FFTstepCountA+1; I2++) {
        perfCyclesA += plan->performance_FFTstepCyclesA[0][I2];
	bestPerfCyclesA += bestPerformance_FFTstepCyclesA[0][I2];
      }
      for (int I2=0; I2<2*plan->FFTstepCountB+1; I2++) {
        perfCyclesB += plan->performance_FFTstepCyclesB[0][I2];
	bestPerfCyclesB += bestPerformance_FFTstepCyclesB[0][I2];
      }
      
      if ((isNaN(bestPerfCyclesA) || (perfCyclesA<bestPerfCyclesA)) && (perfCyclesA>=0)) {
        for (int I=0; I<numberOfThreads; I++) {
          for (int I2=0; I2<2*plan->FFTstepCountA+1; I2++) {
 	    bestPerformance_FFTstepCyclesA[I][I2] = plan->performance_FFTstepCyclesA[I][I2];
	  }
	}
      }
      if ((isNaN(bestPerfCyclesB) || (perfCyclesB<bestPerfCyclesB)) && (perfCyclesB>=0)) {
        for (int I=0; I<numberOfThreads; I++) {
          for (int I2=0; I2<2*plan->FFTstepCountB+1; I2++) {
 	    bestPerformance_FFTstepCyclesB[I][I2] = plan->performance_FFTstepCyclesB[I][I2];
	  }
	}
      }
      if ((perfCyclesA<0) || (perfCyclesB<0)) {
        if (extraMeasurementMax>0) {
	  extraMeasurementMax--;
          if (LogLevel>2) printf("...repeating measurement!\n");
          iter--;
	}
      }
     
      for (int I=0; I<numberOfThreads; I++) {
        for (int I2=0; I2<2*plan->FFTstepCountA+1; I2++) {
	  plan->performance_FFTstepCyclesA[I][I2] = 0;	  
	}
        for (int I2=0; I2<2*plan->FFTstepCountB+1; I2++) {
	  plan->performance_FFTstepCyclesB[I][I2] = 0;
	}
      }
    }
  }

  if (doTiming) {
    for (int I=0; I<numberOfThreads; I++) {
      for (int I2=0; I2<2*plan->FFTstepCountA+1; I2++) plan->performance_FFTstepCyclesA[I][I2] = bestPerformance_FFTstepCyclesA[I][I2];
      for (int I2=0; I2<2*plan->FFTstepCountB+1; I2++) plan->performance_FFTstepCyclesB[I][I2] = bestPerformance_FFTstepCyclesB[I][I2];
      delete[] bestPerformance_FFTstepCyclesA[I];
      delete[] bestPerformance_FFTstepCyclesB[I];      
    }

    delete[] bestPerformance_FFTstepCyclesA;
    delete[] bestPerformance_FFTstepCyclesB;   

    plan->performance_TotalCycles = 0;
    for (int I2=0; I2<2*plan->FFTstepCountA+1; I2++) {
      plan->performance_TotalCycles += plan->performance_FFTstepCyclesA[0][I2];
    }
    for (int I2=0; I2<2*plan->FFTstepCountB+1; I2++) {
      plan->performance_TotalCycles += plan->performance_FFTstepCyclesB[0][I2];
    }
  }

  if (!runThreadsContinously) stopThreadsAndWaitForTermination();
}


void ExtremeFFT4D::ThreadControlLoopRoutine(int ThreadID) {
  long int* runningAddr = &(ThreadControlDataStructure[ThreadID].running);
  long int running;
  long int* commandAddr = &(ThreadControlDataStructure[ThreadID].command);
  long int command = ExtremeFFT4D_ThreadControlCommand_NoCommand;
  long int NoCommand = ExtremeFFT4D_ThreadControlCommand_NoCommand;

  xSSE* xSSEObjLoc = new xSSE();

  running = 1;
  xSSEObjLoc->xSSE_WriteLongIntToMemAddr_Wrapper(runningAddr, running);
  
  xSSEObjLoc->xSSE_ReadLongIntFromMemAddr_Wrapper(runningAddr, running);
  xSSEObjLoc->xSSE_ReadLongIntFromMemAddr_Wrapper(commandAddr, command);
  
  long int mask = CoreMaskForThreads[ThreadID];
  setCurrentThreadAffinityMask(mask);
  
  while (command != ExtremeFFT4D_ThreadControlCommand_StopThread) {
    xSSEObjLoc->xSSE_PerformNopLoop_Wrapper(ExtremeFFT4D_ThreadControlLoopWaitCycles);
    xSSEObjLoc->xSSE_ReadLongIntFromMemAddr_Wrapper(commandAddr, command);
    if (command != NoCommand) xSSEObjLoc->xSSE_WriteLongIntToMemAddr_Wrapper(commandAddr, NoCommand);      
      
    if (command == ExtremeFFT4D_ThreadControlCommand_StopThread) break;
    if (command == ExtremeFFT4D_ThreadControlCommand_PlanExecution) {
      command = ExtremeFFT4D_ThreadControlCommand_NoCommand;
      threadedExecutionOfPlan(ThreadID, xSSEObjLoc);
    }  
  }
  
  running = 0;
  xSSEObjLoc->xSSE_WriteLongIntToMemAddr_Wrapper(runningAddr, running);
  delete xSSEObjLoc;
}


void* ExtremeFFT4D_ThreadControlLoopRoutine(void* para) {
  long int ThreadID = ((long int*)para)[0];
  ExtremeFFT4D* xFFTObj = (ExtremeFFT4D*)(((long int*)para)[1]);
  xFFTObj->ThreadControlLoopRoutine(ThreadID);
  return NULL;
}


void ExtremeFFT4D::startThreads() {
  for (int I=1; I<numberOfThreads; I++) {
    long int* runningAddr = &(ThreadControlDataStructure[I].running);
    long int running = 0;
    xSSEObj->xSSE_ReadLongIntFromMemAddr_Wrapper(runningAddr, running);
    if (running==0) {
      long int* commandAddr = &(ThreadControlDataStructure[I].command);
      long int command = ExtremeFFT4D_ThreadControlCommand_NoCommand;
      xSSEObj->xSSE_WriteLongIntToMemAddr_Wrapper(commandAddr, command);
      
      ThreadControlParameterHolder[I][0] = I;
      ThreadControlParameterHolder[I][1] = (long int) this;
      void* parameterAddr = (void*) &(ThreadControlParameterHolder[I][0]);
      pthread_t thread;
      if (pthread_create(&thread, NULL, ExtremeFFT4D_ThreadControlLoopRoutine, parameterAddr) != 0) {
        printf("ERROR in ExtremeFFT4D: Could not create thread nr %d\n",I);
        exit(0);
      }    
    }
  }
}


void ExtremeFFT4D::stopThreads() {
  for (int I=1; I<numberOfThreads; I++) {
    long int* runningAddr = &(ThreadControlDataStructure[I].running);
    long int running = 0;
    xSSEObj->xSSE_ReadLongIntFromMemAddr_Wrapper(runningAddr, running);
    if (running==1) {
      long int* commandAddr = &(ThreadControlDataStructure[I].command);
      long int command = ExtremeFFT4D_ThreadControlCommand_StopThread;
      xSSEObj->xSSE_WriteLongIntToMemAddr_Wrapper(commandAddr, command);
    }
  }
}


void ExtremeFFT4D::stopThreadsAndWaitForTermination() {
  bool StillRunning = true;
  while (StillRunning) {
    StillRunning = false;
    stopThreads();
    xSSEObj->xSSE_PerformNopLoop_Wrapper(ExtremeFFT4D_ThreadControlLoopWaitCycles);

    for (int I=1; I<numberOfThreads; I++) {
      long int* runningAddr = &(ThreadControlDataStructure[I].running);
      long int running = 0;
      xSSEObj->xSSE_ReadLongIntFromMemAddr_Wrapper(runningAddr, running);
      if (running==1) {
        StillRunning = true;
      }
    }
  }
}


int ExtremeFFT4D::getL0() {
  return L[0];
}


int ExtremeFFT4D::getL1() {
  return L[1];
}

 
int ExtremeFFT4D::getL2() {
  return L[2];
}

 
int ExtremeFFT4D::getL3() {
  return L[3];
}

 
int ExtremeFFT4D::getLocalIndexCount() {
  return localIndexCount;
}


int ExtremeFFT4D::getXtraSize1() {
  return xtrSize1;
}


int ExtremeFFT4D::getXtraSize2() {
  return xtrSize2;
}


int ExtremeFFT4D::getXtraSize3() {
  return xtrSize3;
}


int ExtremeFFT4D::getNumberOfThreads() {
  return numberOfThreads;
}


char* ExtremeFFT4D::getSelectedPlanDescriptor() {
  char* descriptor = new char[10000];
  char* dummy = new char[10000];
  snprintf(descriptor, 10000, "Lind(%d,%d,%d,%d) Iind(%d,%d) OComb(%d,%d) PComb(%d,%d,%d,%d) Iembed(%ld,%d,%d,%ld,%d,%d) block(%d,%d,%d,%d,%d,%d) Scount(%d,%d)",
   selectedFFTplan->LindA1,selectedFFTplan->LindA2,selectedFFTplan->LindB1,selectedFFTplan->LindB2,selectedFFTplan->InnerLoopCountA,selectedFFTplan->InnerLoopCountB,
   selectedFFTplan->OuterLoopCombineFactorA2,selectedFFTplan->OuterLoopCombineFactorB2,
   selectedFFTplan->PrefetchCombineA1,selectedFFTplan->PrefetchCombineA2,selectedFFTplan->PrefetchCombineB1,selectedFFTplan->PrefetchCombineB2,
   selectedFFTplan->inputSlice2DDistanceInBytesA,selectedFFTplan->InternalEmbeddingAOne,selectedFFTplan->InternalEmbeddingATwo,
   selectedFFTplan->inputSlice2DDistanceInBytesB,selectedFFTplan->InternalEmbeddingBOne,selectedFFTplan->InternalEmbeddingBTwo,
   selectedFFTplan->blockSimReadReadA,selectedFFTplan->blockSimReadWriteA,selectedFFTplan->blockSimWriteWriteA,
   selectedFFTplan->blockSimReadReadB,selectedFFTplan->blockSimReadWriteB,selectedFFTplan->blockSimWriteWriteB,
   selectedFFTplan->FFTstepCountA,selectedFFTplan->FFTstepCountB);

  if ((!isNaN(selectedFFTplan->performance_TotalCycles)) && (selectedFFTplan->performance_FFTstepCyclesA!=NULL) && (selectedFFTplan->performance_FFTstepCyclesB!=NULL)) {
    snprintf(dummy, 10000, "%s", descriptor);
    snprintf(descriptor, 10000, "%s Cycles(%1.2f:", dummy, selectedFFTplan->performance_TotalCycles/1E6);

    for (int I=0; I<numberOfThreads; I++) {
      snprintf(dummy, 10000, "%s", descriptor);
      snprintf(descriptor, 10000, "%s T%d((%1.2f,%1.2f)", dummy, I, selectedFFTplan->performance_FFTstepCyclesA[I][1]/1E6, selectedFFTplan->performance_FFTstepCyclesA[I][0]/1E6);
      for (int I2=1; I2<selectedFFTplan->FFTstepCountA; I2++) {
        snprintf(dummy, 10000, "%s", descriptor);
        snprintf(descriptor, 10000, "%s (%1.2f,%1.2f)", dummy, selectedFFTplan->performance_FFTstepCyclesA[I][1+2*I2]/1E6, selectedFFTplan->performance_FFTstepCyclesA[I][0+2*I2]/1E6);
      }
      snprintf(dummy, 10000, "%s", descriptor);
      snprintf(descriptor, 10000, "%s + %1.2f /", dummy, selectedFFTplan->performance_FFTstepCyclesA[I][0+2*selectedFFTplan->FFTstepCountA]/1E6);
      for (int I2=0; I2<selectedFFTplan->FFTstepCountB; I2++) {
        snprintf(dummy, 10000, "%s", descriptor);
        snprintf(descriptor, 10000, "%s (%1.2f,%1.2f)", dummy, selectedFFTplan->performance_FFTstepCyclesB[I][1+2*I2]/1E6, selectedFFTplan->performance_FFTstepCyclesB[I][0+2*I2]/1E6);
      }
      snprintf(dummy, 10000, "%s", descriptor);
      snprintf(descriptor, 10000, "%s + %1.2f)", dummy, selectedFFTplan->performance_FFTstepCyclesB[I][0+2*selectedFFTplan->FFTstepCountB]/1E6);
    }
    
    snprintf(dummy, 10000, "%s", descriptor);
    snprintf(descriptor, 10000, "%s)", dummy);
  }

   
  snprintf(dummy, 10000, "%s", descriptor);
  snprintf(descriptor, 10000, "%s FFTsteps(%s", dummy, availFFTsteps[selectedFFTplan->FFTstepNrA[0]]->getName());
  for (int I=1; I<selectedFFTplan->FFTstepCountA; I++) {
    snprintf(dummy, 10000, "%s", descriptor);
    snprintf(descriptor, 10000, "%s %s", dummy, availFFTsteps[selectedFFTplan->FFTstepNrA[I]]->getName());
  }   
  for (int I=0; I<selectedFFTplan->FFTstepCountB; I++) {
    snprintf(dummy, 10000, "%s", descriptor);
    snprintf(descriptor, 10000, "%s %s", dummy, availFFTsteps[selectedFFTplan->FFTstepNrB[I]]->getName());
  }
  snprintf(dummy, 10000, "%s", descriptor);
  snprintf(descriptor, 10000, "%s )", dummy);

  delete dummy;
  return descriptor;
}


bool ExtremeFFT4D::setSelectedPlanFromDescriptor(char* descriptor) {
  int LindA1, LindA2, LindB1, LindB2;
  int InnerLoopCountA, InnerLoopCountB;
  int OuterLoopCombineFactorA2, OuterLoopCombineFactorB2; 
  int PrefetchCombineA1, PrefetchCombineA2, PrefetchCombineB1, PrefetchCombineB2; 
  long int inputSlice2DDistanceInBytesA, inputSlice2DDistanceInBytesB; 
  int InternalEmbeddingAOne, InternalEmbeddingATwo, InternalEmbeddingBOne, InternalEmbeddingBTwo;          
  int blockSimReadReadA, blockSimReadWriteA, blockSimWriteWriteA;    
  int blockSimReadReadB, blockSimReadWriteB, blockSimWriteWriteB;        
  int FFTstepCountA, FFTstepCountB;
  int FFTstepNrA[2*ExtremeFFT4D_PrimeFactorsMAX+2];
  int FFTstepNrB[2*ExtremeFFT4D_PrimeFactorsMAX+2];
  
  bool planFound = false;
  
  if (26 != sscanf(descriptor, "Lind(%d,%d,%d,%d) Iind(%d,%d) OComb(%d,%d) PComb(%d,%d,%d,%d) Iembed(%ld,%d,%d,%ld,%d,%d) block(%d,%d,%d,%d,%d,%d) Scount(%d,%d)",
   &LindA1,&LindA2,&LindB1,&LindB2,&InnerLoopCountA,&InnerLoopCountB,
   &OuterLoopCombineFactorA2,&OuterLoopCombineFactorB2,
   &PrefetchCombineA1,&PrefetchCombineA2,&PrefetchCombineB1,&PrefetchCombineB2,
   &inputSlice2DDistanceInBytesA,&InternalEmbeddingAOne,&InternalEmbeddingATwo,
   &inputSlice2DDistanceInBytesB,&InternalEmbeddingBOne,&InternalEmbeddingBTwo,
   &blockSimReadReadA,&blockSimReadWriteA,&blockSimWriteWriteA,
   &blockSimReadReadB,&blockSimReadWriteB,&blockSimWriteWriteB,
   &FFTstepCountA,&FFTstepCountB)) {
    return false;
  }
  
  int startPos = 0;
  while ((startPos<((int)strlen(descriptor))) && (descriptor[startPos]!='F')) {
    startPos++;
  }
  startPos += 9;  
  if (startPos>=((int)strlen(descriptor))) return false;
  
  for (int I=0; I<FFTstepCountA; I++) {
    char* dummy = &(descriptor[startPos]);
    char* dummy2 = new char[1000];
    sscanf(dummy, "%s ", dummy2);
    
    bool found = false;
    for (int I2=0; I2<availFFTstepCount; I2++) {
      if (strlen(dummy2) == strlen(availFFTsteps[I2]->getName())) {
        bool diff = false;
        for (int I3=0; I3<((int)strlen(dummy2)); I3++) {
	  if (dummy2[I3]!=availFFTsteps[I2]->getName()[I3]) diff = true;
	}
	if (!diff) {
	  FFTstepNrA[I] = I2;
	  found = true;
	}
      }
    }
    
    startPos += 1 + strlen(dummy2);  
    delete[] dummy2;
    if (!found) return false;
  }
  
  for (int I=0; I<FFTstepCountB; I++) {
    char* dummy = &(descriptor[startPos]);
    char* dummy2 = new char[1000];
    sscanf(dummy, "%s ", dummy2);
    
    bool found = false;
    for (int I2=0; I2<availFFTstepCount; I2++) {
      if (strlen(dummy2) == strlen(availFFTsteps[I2]->getName())) {
        bool diff = false;
        for (int I3=0; I3<((int)strlen(dummy2)); I3++) {
	  if (dummy2[I3]!=availFFTsteps[I2]->getName()[I3]) diff = true;
	}
	if (!diff) {
	  FFTstepNrB[I] = I2;
	  found = true;
	}
      }
    }
    
    startPos += 1 + strlen(dummy2);  
    delete[] dummy2;
    if (!found) return false;
  }
  
  FFTplanType* plan = possibleFFTplanList;
  while (plan != NULL) {
    bool diff = false;
    if (!plan->complete) diff = true;
    if (plan->LindA1 != LindA1) diff = true;
    if (plan->LindA2 != LindA2) diff = true;
    if (plan->LindB1 != LindB1) diff = true;
    if (plan->LindB2 != LindB2) diff = true;
    if (plan->InnerLoopCountA != InnerLoopCountA) diff = true;
    if (plan->InnerLoopCountB != InnerLoopCountB) diff = true;
    if (plan->OuterLoopCombineFactorA2 != OuterLoopCombineFactorA2) diff = true;
    if (plan->OuterLoopCombineFactorB2 != OuterLoopCombineFactorB2) diff = true;
    if (plan->FFTstepCountA != FFTstepCountA) diff = true;
    if (plan->FFTstepCountB != FFTstepCountB) diff = true;
  
    for (int I=0; I<FFTstepCountA; I++) {
      if (FFTstepNrA[I] != plan->FFTstepNrA[I]) diff = true;
    }
    for (int I=0; I<FFTstepCountB; I++) {
      if (FFTstepNrB[I] != plan->FFTstepNrB[I]) diff = true;
    }
  
    if (!diff) {
      dismissPlanDetails(selectedFFTplan, false);
      selectedFFTplan = plan;
      
      planFound = true;
      break;
    }
  
    plan = plan->next;
  }
  
  if (planFound) {
    elaboratePlanDetails(selectedFFTplan);
    
    selectedFFTplan->PrefetchCombineA1 = PrefetchCombineA1;
    selectedFFTplan->PrefetchCombineA2 = PrefetchCombineA2;
    selectedFFTplan->PrefetchCombineB1 = PrefetchCombineB1;
    selectedFFTplan->PrefetchCombineB2 = PrefetchCombineB2;
    selectedFFTplan->inputSlice2DDistanceInBytesA = inputSlice2DDistanceInBytesA;
    selectedFFTplan->inputSlice2DDistanceInBytesB = inputSlice2DDistanceInBytesB;
    selectedFFTplan->InternalEmbeddingAOne = InternalEmbeddingAOne;
    selectedFFTplan->InternalEmbeddingATwo = InternalEmbeddingATwo;
    selectedFFTplan->InternalEmbeddingBOne = InternalEmbeddingBOne;
    selectedFFTplan->InternalEmbeddingBTwo = InternalEmbeddingBTwo;
    selectedFFTplan->blockSimReadReadA = blockSimReadReadA;
    selectedFFTplan->blockSimReadWriteA = blockSimReadWriteA;
    selectedFFTplan->blockSimWriteWriteA = blockSimWriteWriteA;
    selectedFFTplan->blockSimReadReadB = blockSimReadReadB;
    selectedFFTplan->blockSimReadWriteB = blockSimReadWriteB;
    selectedFFTplan->blockSimWriteWriteB = blockSimWriteWriteB;
    
    getPlanReadyForExecution(selectedFFTplan); 

    if (LogLevel>2) {
      printf("\nSelected Plan from descriptor: %s\n", descriptor);
      writeFFTplan(selectedFFTplan, false);
    }
  } else {
    if (LogLevel>2) {
      printf("\nNo Plan selected from descriptor: %s (Keeping old plan)\n", descriptor);
      writeFFTplan(selectedFFTplan, false);
    }
  }
  
  return planFound;
}
