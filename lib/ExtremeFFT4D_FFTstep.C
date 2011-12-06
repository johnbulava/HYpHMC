#include "ExtremeFFT4D_FFTstep.h"

ExtremeFFT4D_FFTstep::ExtremeFFT4D_FFTstep() {
  readFromInput = false; 
  writeToOutput = false; 
  primeFacCount1 = 0; 
  primeFacCount2 = 0;
  for (int I=0; I<ExtremeFFT4D_FFTstep_PrimeFactorsMAX; I++) {
    primeFacs1[I] = 0;
    primeFacs2[I] = 0;
  }
  mustBeFirstTrafo1 = false;
  mustBeFirstTrafo2 = false;        
  readyForExecution = false;
  xSSEObj = NULL;
  L1Cache = new CacheSimulator(L1CacheSizePerWayInBytes*L1Ways, L1Ways, L1CacheLineSizeInBytes);
  Size1 = 0;
  Size2 = 0;
  localIndexCount = 0;  
  InternalEmbeddingValueOne = 0;
  InternalEmbeddingValueTwo = 0;
  readPrefetchCombineFac = 1;
  writePrefetchCombineFac = 1;
  name = new char[1000];
  snprintf(name, 1000, "Generic");

  readPrefetchBaseAddresses = NULL;
  writePrefetchBaseAddresses = NULL;
  readAddresses = NULL;
  writeAddresses = NULL;
  readPrefetchAddresses = NULL;
  writePrefetchAddresses = NULL;
  prefetchDataSpaceHolder = NULL;
  prefetchDataGenerationHelperData1 = NULL;
  prefetchDataGenerationHelperData2 = NULL;
  
  readPrefetchOffsets = NULL;
  writePrefetchOffsets = NULL;
  readPrefetchOffsetList = NULL;
  writePrefetchOffsetList = NULL;
  useReadPrefetching = false;
  useWritePrefetching = false;
}


void ExtremeFFT4D_FFTstep::desini() {
  delete[] readPrefetchBaseAddresses;
  delete[] writePrefetchBaseAddresses;
  delete[] readAddresses;
  delete[] writeAddresses;
  delete[] prefetchDataGenerationHelperData1;
  delete[] prefetchDataGenerationHelperData2;
  delete[] readPrefetchAddresses;
  delete[] writePrefetchAddresses;
  delete[] readPrefetchOffsets;
  delete[] writePrefetchOffsets;
  delete[] readPrefetchOffsetList;
  delete[] writePrefetchOffsetList;
  
  destroySuperAlignedInt(prefetchDataSpaceHolder);  
}



void ExtremeFFT4D_FFTstep::ini() {
  desini();
  if (getNumberOfReadLines()*getNumberOfReadPrefetchesPerTimeIndex()*getNumberOfReadTimeIndicesPerLine()>localIndexCount*Size1*Size2) {
    printf("ERROR in ExtremeFFT4D_FFTstep::ini: Number of time indexes impossible (Read)!\n");
    exit(0);
  }
  if (getNumberOfWriteLines()*getNumberOfWritePrefetchesPerTimeIndex()*getNumberOfWriteTimeIndicesPerLine()>localIndexCount*Size1*Size2) {
    printf("ERROR in ExtremeFFT4D_FFTstep::ini: Number of time indexes impossible (Write)!\n");
    exit(0);
  }

  readPrefetchBaseAddresses = new long int[localIndexCount*Size1*Size2];
  for (int I=0; I<localIndexCount*Size1*Size2; I++) readPrefetchBaseAddresses[I] = 0;
  writePrefetchBaseAddresses = new long int[localIndexCount*Size1*Size2];
  for (int I=0; I<localIndexCount*Size1*Size2; I++) writePrefetchBaseAddresses[I] = 0;
  readAddresses = new long int[2*localIndexCount*Size1*Size2];
  for (int I=0; I<2*localIndexCount*Size1*Size2; I++) readAddresses[I] = 0;
  writeAddresses = new long int[2*localIndexCount*Size1*Size2];
  for (int I=0; I<2*localIndexCount*Size1*Size2; I++) writeAddresses[I] = 0;
  readPrefetchAddresses = new long int[2*localIndexCount*Size1*Size2];
  for (int I=0; I<2*localIndexCount*Size1*Size2; I++) readPrefetchAddresses[I] = 0;
  writePrefetchAddresses = new long int[2*localIndexCount*Size1*Size2];
  for (int I=0; I<2*localIndexCount*Size1*Size2; I++) writePrefetchAddresses[I] = 0;
  prefetchDataSpaceHolder = (int*) createSuperAlignedComplex(localIndexCount*Size1*Size2);
  for (int I=0; I<2*localIndexCount*Size1*Size2; I++) prefetchDataSpaceHolder[I] = 0;
  prefetchDataGenerationHelperData1 = new long int[localIndexCount*Size1*Size2];
  for (int I=0; I<localIndexCount*Size1*Size2; I++) prefetchDataGenerationHelperData1[I] = 0;
  prefetchDataGenerationHelperData2 = new long int[localIndexCount*Size1*Size2];
  for (int I=0; I<localIndexCount*Size1*Size2; I++) prefetchDataGenerationHelperData2[I] = 0;
  
  readPrefetchOffsets = new int*[Size1];
  for (int I=0; I<Size1; I++) readPrefetchOffsets[I] = NULL;
  writePrefetchOffsets = new int*[Size1];
  for (int I=0; I<Size1; I++) writePrefetchOffsets[I] = NULL;
  readPrefetchOffsetList = new int*[Size1];
  for (int I=0; I<Size1; I++) readPrefetchOffsetList[I] = NULL;
  writePrefetchOffsetList = new int*[Size1];
  for (int I=0; I<Size1; I++) writePrefetchOffsetList[I] = NULL;
  
  int index = 0;
  for (int I=0; I<Size1; I++) {
    readPrefetchOffsets[I] = & (prefetchDataSpaceHolder[index]);
    index += getNumberOfReadPrefetchesPerTimeIndex()*getNumberOfReadTimeIndicesPerLine();
    writePrefetchOffsets[I] = & (prefetchDataSpaceHolder[index]);
    index += getNumberOfWritePrefetchesPerTimeIndex()*getNumberOfWriteTimeIndicesPerLine();    
    readPrefetchOffsetList[I] = readPrefetchOffsets[I];
    writePrefetchOffsetList[I] = writePrefetchOffsets[I];
  }
  
  calcReadAndReadPrefetchBaseAddresses();
  calcWriteAndWritePrefetchBaseAddresses();
  generatePrefetchData();  
}


ExtremeFFT4D_FFTstep::~ExtremeFFT4D_FFTstep() {
  desini();
  delete xSSEObj;
  delete[] name;
  delete L1Cache;
}


bool ExtremeFFT4D_FFTstep::getReadFromInput() {
  return readFromInput;
}


bool ExtremeFFT4D_FFTstep::getWriteToOutput() {
  return writeToOutput;
}

 
int ExtremeFFT4D_FFTstep::getPrimeFacCount1()  {
  return primeFacCount1;
}


int ExtremeFFT4D_FFTstep::getPrimeFacCount2() {
  return primeFacCount2;
}


int* ExtremeFFT4D_FFTstep::getPrimeFacs1() {
  return primeFacs1;
}


int* ExtremeFFT4D_FFTstep::getPrimeFacs2() {
  return primeFacs2;
}


bool ExtremeFFT4D_FFTstep::getMustBeFirstTrafo1() {
  return mustBeFirstTrafo1;
}


bool ExtremeFFT4D_FFTstep::getMustBeFirstTrafo2()  {
  return mustBeFirstTrafo2;
}


bool ExtremeFFT4D_FFTstep::getReadyForExecution() {
  return readyForExecution;
}

       
void ExtremeFFT4D_FFTstep::generatePrefetchData(int PrefetchCombineFac, int loadsPerTimeIndex, int TimeIndicesPerLine, int lines, long int* touchedAddresses, long int* prefetchAddresses, long int* prefetchBaseAddresses, int** prefetchOffsets, int** prefetchOffsetsList) {
  //Determine max. number of different cache-lines per TimeIndex
  //+ find all different cache-lines
  int maxNrOfDiffCacheLinesPerTimeIndex = 0;
  int nrOfDiffCacheLines = 0;  
  for (int I=0; I<localIndexCount*Size1*Size2; I++) prefetchDataGenerationHelperData1[I] = 0;
  for (int I=0; I<localIndexCount*Size1*Size2; I++) prefetchDataGenerationHelperData2[I] = -1;

  for (int I=0; I<localIndexCount*Size1*Size2; I++) {
   int timeIndex = touchedAddresses[2*I+0];
    int addr = touchedAddresses[2*I+1];
    addr /= L1CacheLineSizeInBytes;
    addr *= L1CacheLineSizeInBytes;
   
    bool found = false;
    for (int I2=nrOfDiffCacheLines-1; I2>=0; I2--) {
      if (prefetchDataGenerationHelperData2[I2] == addr) {
        found = true;
	break;
      }
    }

    if (!found) {
      prefetchDataGenerationHelperData2[nrOfDiffCacheLines] = addr;
      nrOfDiffCacheLines++;
      prefetchDataGenerationHelperData1[timeIndex]++;
    }
  }
  
  for (int I=0; I<localIndexCount*Size1*Size2; I++) {
    if (prefetchDataGenerationHelperData1[I] > maxNrOfDiffCacheLinesPerTimeIndex) {
      maxNrOfDiffCacheLinesPerTimeIndex = prefetchDataGenerationHelperData1[I];
    }
  }

  //Sort touched cache-lines within blocks of size PrefetchCombineBufferSize
  int PrefetchCombineBufferSize = PrefetchCombineFac*maxNrOfDiffCacheLinesPerTimeIndex;
  if (PrefetchCombineBufferSize>0) {
    for (int I=0; I+PrefetchCombineBufferSize<=Size1*Size2*localIndexCount; I+=PrefetchCombineBufferSize) {
      bool changed = true;
      while (changed) {
        changed = false;
        for (int I2=I; I2<I+PrefetchCombineBufferSize-1; I2++) {
          if (((prefetchDataGenerationHelperData2[I2]>prefetchDataGenerationHelperData2[I2+1]) && (prefetchDataGenerationHelperData2[I2+1]!=-1)) || ((prefetchDataGenerationHelperData2[I2]==-1) && (prefetchDataGenerationHelperData2[I2+1]!=-1))) {
	    long int dummy = prefetchDataGenerationHelperData2[I2];
	    prefetchDataGenerationHelperData2[I2] = prefetchDataGenerationHelperData2[I2+1];
  	    prefetchDataGenerationHelperData2[I2+1] = dummy;
	    changed = true;
          }      
        }
      }
    }   
  }

  //Select Addresses for prefetching and construct prefetching plan
  int LoadCapacity = lines*TimeIndicesPerLine*loadsPerTimeIndex;
  bool fractionLoad = (nrOfDiffCacheLines>LoadCapacity);
  bool overLoad= (nrOfDiffCacheLines<LoadCapacity);
  if (LogLevel>1) {
    if (fractionLoad) printf("WARNING: ExtremeFFT4D in fractional Prefetch-Mode (%s)!!!\n",name);
    if (overLoad) printf("WARNING: ExtremeFFT4D in overload Prefetch-Mode (%s)!!!\n",name);
  }
  int timeIndex = 0;
  int count = 0;
  for (int I=0; I<lines*TimeIndicesPerLine; I++) {
    for (int I2=0; I2<loadsPerTimeIndex; I2++) {
      prefetchAddresses[count] = timeIndex;
      prefetchAddresses[count+1] = -1;
      count += 2;
    }
    timeIndex++;
  }
  if (!fractionLoad) {
    for (int I=0; I<nrOfDiffCacheLines-PrefetchCombineBufferSize; I++) {
      prefetchAddresses[2*I+1] = prefetchDataGenerationHelperData2[I+PrefetchCombineBufferSize];
    }
  } else {
    int loadCount = (PrefetchCombineBufferSize*LoadCapacity)/nrOfDiffCacheLines;
    int p = 0;
    for (int I=PrefetchCombineBufferSize; I<nrOfDiffCacheLines; I+=PrefetchCombineBufferSize) {
      for (int I2=0; I2<loadCount; I2++) {
        if (p<LoadCapacity-PrefetchCombineBufferSize) {
	  prefetchAddresses[2*p+1] = prefetchDataGenerationHelperData2[I+I2];
	  p++;
	}
      }
    }
  }
  
  //Calculate prefetch-offsets
  count = 0;
  for (int I=0; I<lines; I++) {
    int count2 = 0;
    for (int I2=0; I2<TimeIndicesPerLine; I2++) {
      for (int I3=0; I3<loadsPerTimeIndex; I3++) {
        long int pBaseAddr = prefetchBaseAddresses[prefetchAddresses[2*count+0]];
	long int addr = prefetchAddresses[2*count+1];
	
	if (addr>=0) {
	  prefetchOffsets[I][count2] = addr-pBaseAddr;	  
	} else {
	  prefetchOffsets[I][count2] = 0;
	}
	
	count++;
	count2++;
      }
    }
  }
  
  //Find different prefetch-data and define offset-list
  int differentPrefetchDataCount=0;
  for (int I=0; I<lines; I++) {
    int sameIndex = -1;
    for (int I2=0; I2<differentPrefetchDataCount; I2++) {
      bool different = false;
      for (int I3=0; I3<TimeIndicesPerLine*loadsPerTimeIndex; I3++) {
        if (prefetchOffsets[I][I3] != prefetchOffsetsList[I2][I3]) {
	  different = true;
	  break;
	}
      }
      if (!different) {
        sameIndex = I2;
	break;
      }
    }
    if (sameIndex>=0) {
      prefetchOffsetsList[I] = prefetchOffsetsList[sameIndex];
    } else {
      prefetchOffsetsList[I] = prefetchOffsets[I];
      differentPrefetchDataCount++;
    }
  }
}


void ExtremeFFT4D_FFTstep::generatePrefetchData() {
  if (useReadPrefetching) generatePrefetchData(readPrefetchCombineFac, getNumberOfReadPrefetchesPerTimeIndex(), getNumberOfReadTimeIndicesPerLine(), getNumberOfReadLines(), readAddresses, readPrefetchAddresses, readPrefetchBaseAddresses, readPrefetchOffsets, readPrefetchOffsetList);
  if (useWritePrefetching) generatePrefetchData(writePrefetchCombineFac, getNumberOfWritePrefetchesPerTimeIndex(), getNumberOfWriteTimeIndicesPerLine(), getNumberOfWriteLines(), writeAddresses, writePrefetchAddresses, writePrefetchBaseAddresses, writePrefetchOffsets, writePrefetchOffsetList);
  generateLocalPrefetchData();
}


void ExtremeFFT4D_FFTstep::setInternalEmbedding(int intEmbeddOne, int intEmbeddTwo) {
  InternalEmbeddingValueOne = intEmbeddOne;
  InternalEmbeddingValueTwo = intEmbeddTwo;
  changedInternalEmbedding();
  if (!readFromInput) calcReadAndReadPrefetchBaseAddresses();
  if (!writeToOutput) calcWriteAndWritePrefetchBaseAddresses();
  generatePrefetchData();  
}


void ExtremeFFT4D_FFTstep::setPrefetchCombineFacs(int readComF, int writeComF) {
  readPrefetchCombineFac = readComF;
  writePrefetchCombineFac = writeComF;
  generatePrefetchData();
}


void ExtremeFFT4D_FFTstep::calcL1ExcessMisses(long int sizeOfWhole4DFieldInBytes, long int inputSlice2DDistanceInBytes, int &minorMisses, int &majorMisses, bool abortWhenResultWorse) {
  L1Cache->reset();
  
  sizeOfWhole4DFieldInBytes /= L1CacheSizePerWayInBytes;
  sizeOfWhole4DFieldInBytes *= L1CacheSizePerWayInBytes;
  sizeOfWhole4DFieldInBytes += L1CacheSizePerWayInBytes;
  
  long int writeOffset = 0;
  if (readFromInput) writeOffset = sizeOfWhole4DFieldInBytes + inputSlice2DDistanceInBytes;
  bool considerWrite = true;
  if (writeToOutput) considerWrite = false;
  int minimalLoads = (localIndexCount*Size1*Size2*16) / L1CacheLineSizeInBytes;
  if (readFromInput) minimalLoads *= 2;
  
  int timeIndex = 0;
  bool dataAccessed = true;
  int readIndex = 0;
  int writeIndex = 0;
  int readPefetchIndex = 0;
  int writePrefetchIndex = 0;  
  while (dataAccessed) {
    dataAccessed = false;
    
    for (int I=readIndex; I<localIndexCount*Size1*Size2; I++) {
      if (readAddresses[2*I+0] == timeIndex) {
        long int addr = readAddresses[2*I+1];
        L1Cache->accessAddress(addr);      
	dataAccessed = true;
      } else {
        readIndex = I;
	break;
      }
    }

    if (useReadPrefetching) {
      for (int I=readPefetchIndex; I<localIndexCount*Size1*Size2; I++) {
        if (readPrefetchAddresses[2*I+0] == timeIndex) {
          long int addr = readPrefetchAddresses[2*I+1];
          if (addr<0) addr = readPrefetchBaseAddresses[timeIndex];
          L1Cache->accessAddress(addr);       
  	  dataAccessed = true;
        } else {
          readPefetchIndex = I;
  	  break;
        }
      }
    }

    if (useWritePrefetching) {
      for (int I=writePrefetchIndex; I<localIndexCount*Size1*Size2; I++) {
        if (writePrefetchAddresses[2*I+0] == timeIndex) {
          long int addr = writePrefetchAddresses[2*I+1];
          if (addr<0) addr = writePrefetchBaseAddresses[timeIndex];
	  addr += writeOffset;
          L1Cache->accessAddress(addr);      
	  dataAccessed = true;
        } else {
          writePrefetchIndex = I; 
	  break;
        }
      }
    }
    
    if (considerWrite) {
      for (int I=writeIndex; I<localIndexCount*Size1*Size2; I++) {
        if (writeAddresses[2*I+0] == timeIndex) {
          long int addr = writeAddresses[2*I+1];
  	  addr += writeOffset;
          L1Cache->accessAddress(addr);      
  	  dataAccessed = true;
        } else {
          writeIndex = I;  
	  break;
        }
      }
    }
    if (abortWhenResultWorse) {
      if ((minorMisses>=0) && (L1Cache->getLastAccessMisses() - minimalLoads > minorMisses)) break;
      if ((majorMisses>=0) && (L1Cache->getMisses() - minimalLoads > majorMisses)) break;
    }
    
    timeIndex++;
  }

  minorMisses = L1Cache->getLastAccessMisses() - minimalLoads;
  majorMisses = L1Cache->getMisses() - minimalLoads;
}


char* ExtremeFFT4D_FFTstep::getName() {
  return name;
}
