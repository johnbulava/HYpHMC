#include "ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2.h"

ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2::ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2() : ExtremeFFT4D_FFTstep() {
  readFromInput = true; 
  writeToOutput = false; 
  primeFacCount1 = 2; 
  primeFacCount2 = 2;
  primeFacs1[0] = 2;
  primeFacs1[1] = 2;
  primeFacs2[0] = 2;
  primeFacs2[1] = 2;
  mustBeFirstTrafo1 = true;
  mustBeFirstTrafo2 = true;        
  AddressAuxData1 = NULL;
  AddressAuxData2 = NULL;  
  addressIncrementsInBytes1 = 0;
  addressIncrementsInBytes2 = 0;
  long int* BitMask_ComplexForward_Pointer = (long int*)&BitMask_ComplexForward;
  long int* BitMask_ComplexBackward_Pointer = (long int*)&BitMask_ComplexBackward;
  BitMask_ComplexBackward_Pointer[0] = (long int) 0x8000000000000000;
  BitMask_ComplexBackward_Pointer[1] = (long int) 0;
  BitMask_ComplexForward_Pointer[0] = (long int) 0;
  BitMask_ComplexForward_Pointer[1] = (long int) 0x8000000000000000;
  snprintf(name, 1000, "Read2DSliceAndDoFirstAndSecond2DFFTStepBase2");
  combinedPrefetchOffsetList = NULL;
  useReadPrefetching = true;
  useWritePrefetching = true;
}


ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2::~ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2() {
  if (AddressAuxData1 != NULL) {
    destroySuperAlignedLongInt(AddressAuxData1);
  }
  if (AddressAuxData2 != NULL) {
    destroySuperAlignedLongInt(AddressAuxData2);  
  }
  delete[] combinedPrefetchOffsetList;
}


ExtremeFFT4D_FFTstep* ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2::deriveSpecificFFTstepInstance(int siz1, int siz2, int locInd, int add1, int add2, int alreadyPerfPrimeFacProd1, int alreadyPerfPrimeFacProd2, long int* invAddrBase1, long int* invAddrBase2) {
  ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2* fftStep = new ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2();
  fftStep->setValuesInDerivedSpecificFFTstepInstance(siz1, siz2, locInd, add1, add2, alreadyPerfPrimeFacProd1, alreadyPerfPrimeFacProd2, invAddrBase1, invAddrBase2);
  return fftStep;
}


void ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2::setValuesInDerivedSpecificFFTstepInstance(int siz1, int siz2, int locInd, int add1, int add2, int alreadyPerfPrimeFacProd1, int alreadyPerfPrimeFacProd2, long int* invAddrBase1, long int* invAddrBase2) {
  Size1 = siz1;
  Size2 = siz2;
  localIndexCount = locInd;  
  addressIncrementsInBytes1 = add1;
  addressIncrementsInBytes2 = add2;
  if (xSSEObj != NULL) delete xSSEObj;
  xSSEObj = new xSSE();
  if (AddressAuxData1 != NULL) {
    destroySuperAlignedLongInt(AddressAuxData1);
  }
  if (AddressAuxData2 != NULL) {
    destroySuperAlignedLongInt(AddressAuxData2);  
  }  
  AddressAuxData1 = (long int*)createSuperAlignedComplex(Size1);
  AddressAuxData2 = (long int*)createSuperAlignedComplex(Size2); 

  delete[] combinedPrefetchOffsetList;
  combinedPrefetchOffsetList = new int*[Size1];

  for (int I=0; I<Size1; I++) {
    AddressAuxData1[I] = invAddrBase1[Size1-1-I];
  }
  int lastIndex = 16*localIndexCount;
  for (int I=0; I<Size2; I++) {
    AddressAuxData2[I] = 0;
  }
  for (int I=0; I<Size2/4; I++) {
    int p = 4*16*localIndexCount*(invAddrBase2[(4*(I+1)) % Size2]/addressIncrementsInBytes2);
    AddressAuxData2[Size2-1-I*4] = p - lastIndex;
    lastIndex = p + 16*localIndexCount;
  } 
  
  ini();
}


void ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2::executeFFTstep(Complex* pIn, Complex* pOut, Complex* slice2D, bool forward) {
  Complex facpmI = BitMask_ComplexForward;
  if (!forward) facpmI = BitMask_ComplexBackward;

  xSSEObj->xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_Wrapper(pIn, slice2D, Size1, Size2, localIndexCount, AddressAuxData1, AddressAuxData2, combinedPrefetchOffsetList, addressIncrementsInBytes2, InternalEmbeddingValueOne, InternalEmbeddingValueTwo, facpmI);
}


void ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2::calcReadAndReadPrefetchBaseAddresses() {
  for (int I=0; I<localIndexCount*Size1*Size2; I++) readPrefetchBaseAddresses[I] = -1;
  for (int I=0; I<2*localIndexCount*Size1*Size2; I++) readAddresses[I] = -1;
  
  int addressCount = 0;  
  int timeIndex = 0;
  int baseStep1 = addressIncrementsInBytes2 * Size2/4;
  int baseStep2 = addressIncrementsInBytes2;
  for (int I=Size1-1; I>=0; I-=4) {
    for (int I2=0; I2<Size2/4; I2++) {
      for (int I3=0; I3<localIndexCount; I3++) {
        readPrefetchBaseAddresses[timeIndex] = AddressAuxData1[I];
        for (int i=0; i<4; i++) {
	  int addrBase = AddressAuxData1[I-i];
	  for (int i2=0; i2<4; i2++) {
	    long int addr = addrBase + i2*baseStep1 + I2*baseStep2 + I3*16;
	    readAddresses[2*addressCount+0] = timeIndex;
	    readAddresses[2*addressCount+1] = addr;
	    addressCount++;
	  }
	}
	timeIndex++;	
      }
    }
  }
}


void ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2::calcWriteAndWritePrefetchBaseAddresses() {
  for (int I=0; I<localIndexCount*Size1*Size2; I++) writePrefetchBaseAddresses[I] = -1;
  for (int I=0; I<2*localIndexCount*Size1*Size2; I++) writeAddresses[I] = -1;

  int addressCount = 0;
  int timeIndex = 0;
  long int rdiBase = 0;
  int baseStep1 = 16 * localIndexCount;
  int baseStep2 = 16 * localIndexCount * Size2 + 64*InternalEmbeddingValueOne;
  for (int I=Size1-1; I>=0; I-=4) {
    for (int I2=0; I2<Size2/4; I2++) {
      for (int I3=0; I3<localIndexCount; I3++) {
        writePrefetchBaseAddresses[timeIndex] = rdiBase;
        for (int i=0; i<4; i++) {
	  for (int i2=0; i2<4; i2++) {
	    long int addr = rdiBase + i2*baseStep1 + i*baseStep2;

	    writeAddresses[2*addressCount+0] = timeIndex;
	    writeAddresses[2*addressCount+1] = addr;
	    addressCount++;
	  }
	}
	timeIndex++;	
	rdiBase += 16;
      }
      rdiBase += AddressAuxData2[Size2-1-4*I2];
    }
    rdiBase += 4*Size2*localIndexCount*16 + 64*InternalEmbeddingValueTwo + 4*64*InternalEmbeddingValueOne;
  }
}


int ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2::getNumberOfReadPrefetchesPerTimeIndex() {
  return 4;
}


int ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2::getNumberOfWritePrefetchesPerTimeIndex()  {
  return 4;
}


int ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2::getNumberOfReadTimeIndicesPerLine() {
  return (localIndexCount*Size2)/4;
} 


int ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2::getNumberOfWriteTimeIndicesPerLine() {
  return (localIndexCount*Size2)/4;
} 


void ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2::generateLocalPrefetchData() {
  for (int I=0; I<Size1; I++) {
    combinedPrefetchOffsetList[I] = NULL;
  }
  for (int I=0; I<Size1/4; I++) {
    combinedPrefetchOffsetList[Size1-1-4*I] = readPrefetchOffsetList[I];
    combinedPrefetchOffsetList[Size1-2-4*I] = writePrefetchOffsetList[I];    
  }  
}


int ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2::getNumberOfReadLines(){
  return Size1/4;
}


int ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2::getNumberOfWriteLines() {
  return Size1/4;
}


void ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2::changedInternalEmbedding() {
}
