#include "ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2.h"

ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2::ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2() : ExtremeFFT4D_FFTstep() {
  readFromInput = true; 
  writeToOutput = false; 
  primeFacCount1 = 0; 
  primeFacCount2 = 3;
  primeFacs2[0] = 2;
  primeFacs2[1] = 2;
  primeFacs2[2] = 2;  
  mustBeFirstTrafo1 = true;
  mustBeFirstTrafo2 = true;        
  AddressAuxData1 = NULL;
  AddressAuxData2 = NULL;  
  LineNumbers1 = NULL;
  FFTcomplexFactorsForward2 = NULL;
  FFTcomplexFactorsBackward2 = NULL;
  addressIncrementsInBytes1 = 0;
  addressIncrementsInBytes2 = 0;
  long int* BitMask_ComplexForward_Pointer = (long int*)&BitMask_ComplexForward;
  long int* BitMask_ComplexBackward_Pointer = (long int*)&BitMask_ComplexBackward;
  BitMask_ComplexBackward_Pointer[0] = (long int) 0x8000000000000000;
  BitMask_ComplexBackward_Pointer[1] = (long int) 0;
  BitMask_ComplexForward_Pointer[0] = (long int) 0;
  BitMask_ComplexForward_Pointer[1] = (long int) 0x8000000000000000;
  snprintf(name, 1000, "Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2");
  combinedPrefetchOffsetList = NULL;
  useReadPrefetching = true;
  useWritePrefetching = true;
}


ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2::~ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2() {
  if (FFTcomplexFactorsForward2 != NULL) destroySuperAlignedComplex(FFTcomplexFactorsForward2);
  if (FFTcomplexFactorsBackward2 != NULL) destroySuperAlignedComplex(FFTcomplexFactorsBackward2);
  if (AddressAuxData1 != NULL) {
    destroySuperAlignedLongInt(AddressAuxData1);
  }
  if (AddressAuxData2 != NULL) {
    destroySuperAlignedLongInt(AddressAuxData2);  
  }
  if (LineNumbers1 != NULL) {
    destroySuperAlignedLongInt(LineNumbers1);
  }
  delete[] combinedPrefetchOffsetList;
}


ExtremeFFT4D_FFTstep* ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2::deriveSpecificFFTstepInstance(int siz1, int siz2, int locInd, int add1, int add2, int alreadyPerfPrimeFacProd1, int alreadyPerfPrimeFacProd2, long int* invAddrBase1, long int* invAddrBase2) {
  ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2* fftStep = new ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2();
  fftStep->setValuesInDerivedSpecificFFTstepInstance(siz1, siz2, locInd, add1, add2, alreadyPerfPrimeFacProd1, alreadyPerfPrimeFacProd2, invAddrBase1, invAddrBase2);
  return fftStep;
}


void ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2::setValuesInDerivedSpecificFFTstepInstance(int siz1, int siz2, int locInd, int add1, int add2, int alreadyPerfPrimeFacProd1, int alreadyPerfPrimeFacProd2, long int* invAddrBase1, long int* invAddrBase2) {
  Size1 = siz1;
  Size2 = siz2;  
  localIndexCount = locInd;  
  addressIncrementsInBytes1 = add1;
  addressIncrementsInBytes2 = add2;
  if (xSSEObj != NULL) delete xSSEObj;
  xSSEObj = new xSSE;
  if (AddressAuxData1 != NULL) {
    destroySuperAlignedLongInt(AddressAuxData1);
  }
  if (AddressAuxData2 != NULL) {
    destroySuperAlignedLongInt(AddressAuxData2);  
  }  
  if (LineNumbers1 != NULL) {
    destroySuperAlignedLongInt(LineNumbers1);
  }
  AddressAuxData1 = (long int*)createSuperAlignedComplex(Size1);
  AddressAuxData2 = (long int*)createSuperAlignedComplex(Size2); 
  LineNumbers1 = (long int*)createSuperAlignedComplex(Size1);

  delete[] combinedPrefetchOffsetList;
  combinedPrefetchOffsetList = new int*[2*Size1];

  for (int I=0; I<Size1; I++) {
    LineNumbers1[I] = invAddrBase1[Size1-1-I] / addressIncrementsInBytes1;
    AddressAuxData1[I] = LineNumbers1[I]*(Size2*16*localIndexCount+InternalEmbeddingValueOne*64);
    AddressAuxData1[I] += (LineNumbers1[I]/4)*InternalEmbeddingValueTwo*64;
  }
  int lastIndex = 16*localIndexCount;
  for (int I=0; I<Size2; I++) {
    AddressAuxData2[I] = 0;
  }
  for (int I=0; I<Size2/8; I++) {
    int p = 8*16*localIndexCount*(invAddrBase2[(8*(I+1)) % Size2]/addressIncrementsInBytes2);
    AddressAuxData2[Size2-1-I*8] = p - lastIndex;
    lastIndex = p + 16*localIndexCount;
  } 
  
  if (FFTcomplexFactorsForward2 != NULL) destroySuperAlignedComplex(FFTcomplexFactorsForward2);
  if (FFTcomplexFactorsBackward2 != NULL) destroySuperAlignedComplex(FFTcomplexFactorsBackward2);  
  FFTcomplexFactorsForward2 = createSuperAlignedComplex(2);
  FFTcomplexFactorsBackward2 = createSuperAlignedComplex(2);
  
  double fac = -1;
  Complex x = Complex(cos(fac*pi/4), sin(fac*pi/4));
  FFTcomplexFactorsForward2[0].x = x.x;
  FFTcomplexFactorsForward2[0].y = x.x;
  FFTcomplexFactorsForward2[1].x = -x.y;
  FFTcomplexFactorsForward2[1].y = x.y;
  FFTcomplexFactorsBackward2[0].x = x.x;
  FFTcomplexFactorsBackward2[0].y = x.x;
  FFTcomplexFactorsBackward2[1].x = x.y;
  FFTcomplexFactorsBackward2[1].y = -x.y;
  
  ini();
}


void ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2::executeFFTstep(Complex* pIn, Complex* pOut, Complex* slice2D, bool forward) {
  Complex facpmI = BitMask_ComplexForward;
  if (!forward) facpmI = BitMask_ComplexBackward;
  Complex* FFTcomplexFacs2 = FFTcomplexFactorsForward2;
  if (!forward) {
    FFTcomplexFacs2 = FFTcomplexFactorsBackward2;
  }

  xSSEObj->xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2_Wrapper(pIn, slice2D, Size1, Size2, localIndexCount, AddressAuxData1, AddressAuxData2, combinedPrefetchOffsetList, addressIncrementsInBytes1, addressIncrementsInBytes2, FFTcomplexFacs2, facpmI);
}


void ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2::calcReadAndReadPrefetchBaseAddresses() {
  for (int I=0; I<localIndexCount*Size1*Size2; I++) readPrefetchBaseAddresses[I] = -1;
  for (int I=0; I<2*localIndexCount*Size1*Size2; I++) readAddresses[I] = -1;
  
  int addressCount = 0;  
  int timeIndex = 0;
  long int rsi = 0;
  long int baseStep2 = addressIncrementsInBytes2 * Size2/8;  
  for (int I=Size1-1; I>=0; I-=1) {
    for (int I2=0; I2<Size2/8; I2++) {
      for (int I3=0; I3<localIndexCount; I3++) {
        readPrefetchBaseAddresses[timeIndex] = rsi;
	long int addrBase = rsi;	
        for (int i=0; i<8; i++) {
	  long int addr = addrBase + i*baseStep2;
	  readAddresses[2*addressCount+0] = timeIndex;
	  readAddresses[2*addressCount+1] = addr;
	  addressCount++;
	}
	timeIndex++;	
	rsi += 16;
      }
      rsi += addressIncrementsInBytes2 - 16*localIndexCount;
    }
    rsi += addressIncrementsInBytes1 - (Size2/8)*addressIncrementsInBytes2;
  }
}


void ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2::calcWriteAndWritePrefetchBaseAddresses() {
  for (int I=0; I<localIndexCount*Size1*Size2; I++) writePrefetchBaseAddresses[I] = -1;
  for (int I=0; I<2*localIndexCount*Size1*Size2; I++) writeAddresses[I] = -1;

  int addressCount = 0;
  int timeIndex = 0;
  long int rdi = 0;
  long int baseStep2 = 16 * localIndexCount;
  for (int I=Size1-1; I>=0; I-=1) {
    rdi = AddressAuxData1[I];  
    for (int I2=0; I2<Size2/8; I2++) {
      for (int I3=0; I3<localIndexCount; I3++) {
        writePrefetchBaseAddresses[timeIndex] = rdi;
        for (int i=0; i<8; i++) {
          long int addr = rdi + i*baseStep2;

          writeAddresses[2*addressCount+0] = timeIndex;
	  writeAddresses[2*addressCount+1] = addr;
	  addressCount++;
	}
	timeIndex++;	
	rdi += 16;
      }
      rdi += AddressAuxData2[Size2-1-8*I2];
    }
  }
}


int ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2::getNumberOfReadPrefetchesPerTimeIndex() {
  return 2;
}


int ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2::getNumberOfWritePrefetchesPerTimeIndex()  {
  return 2;
}


int ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2::getNumberOfReadTimeIndicesPerLine() {
  return (localIndexCount*Size2)/8;
} 


int ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2::getNumberOfWriteTimeIndicesPerLine() {
  return (localIndexCount*Size2)/8;
} 


void ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2::generateLocalPrefetchData() {
  for (int I=0; I<2*Size1; I++) {
    combinedPrefetchOffsetList[I] = NULL;
  }
  for (int I=0; I<Size1; I++) {
    combinedPrefetchOffsetList[2*Size1-2-2*I] = readPrefetchOffsetList[I];
    combinedPrefetchOffsetList[2*Size1-1-2*I] = writePrefetchOffsetList[I];    
  }  
}


int ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2::getNumberOfReadLines(){
  return Size1;
}


int ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2::getNumberOfWriteLines() {
  return Size1;
}


void ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2::changedInternalEmbedding() {
  for (int I=0; I<Size1; I++) {
    AddressAuxData1[I] = LineNumbers1[I]*(Size2*16*localIndexCount+InternalEmbeddingValueOne*64);
    AddressAuxData1[I] += (LineNumbers1[I]/4)*InternalEmbeddingValueTwo*64;
  }
}
