#include "ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2.h"

ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2::ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2() : ExtremeFFT4D_FFTstep() {
  readFromInput = false; 
  writeToOutput = false; 
  primeFacCount1 = 3; 
  primeFacCount2 = 0;
  primeFacs1[0] = 2;
  primeFacs1[1] = 2;
  primeFacs1[2] = 2;
  mustBeFirstTrafo1 = true;
  mustBeFirstTrafo2 = false;        
  FFTcomplexFactorsForward1 = NULL;
  FFTcomplexFactorsBackward1 = NULL;
  alreadyPerformedPrimeFacProducts1 = 0;
  alreadyPerformedPrimeFacProducts2 = 0;
  long int* BitMask_ComplexForward_Pointer = (long int*)&BitMask_ComplexForward;
  long int* BitMask_ComplexBackward_Pointer = (long int*)&BitMask_ComplexBackward;
  BitMask_ComplexBackward_Pointer[0] = (long int) 0x8000000000000000;
  BitMask_ComplexBackward_Pointer[1] = (long int) 0;
  BitMask_ComplexForward_Pointer[0] = (long int) 0;
  BitMask_ComplexForward_Pointer[1] = (long int) 0x8000000000000000;
  snprintf(name, 1000, "DoFirstSecondAndThird1DD1FFTStepBase2");
  combinedPrefetchOffsetList = NULL;
  useReadPrefetching = true;
  useWritePrefetching = false;
}


ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2::~ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2() {
  if (FFTcomplexFactorsForward1 != NULL) destroySuperAlignedComplex(FFTcomplexFactorsForward1);
  if (FFTcomplexFactorsBackward1 != NULL) destroySuperAlignedComplex(FFTcomplexFactorsBackward1);
  delete[] combinedPrefetchOffsetList;
}


ExtremeFFT4D_FFTstep* ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2::deriveSpecificFFTstepInstance(int siz1, int siz2, int locInd, int add1, int add2, int alreadyPerfPrimeFacProd1, int alreadyPerfPrimeFacProd2, long int* invAddrBase1, long int* invAddrBase2) {
  ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2* fftStep = new ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2();
  fftStep->setValuesInDerivedSpecificFFTstepInstance(siz1, siz2, locInd, add1, add2, alreadyPerfPrimeFacProd1, alreadyPerfPrimeFacProd2, invAddrBase1, invAddrBase2);
  return fftStep;
}


void ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2::setValuesInDerivedSpecificFFTstepInstance(int siz1, int siz2, int locInd, int add1, int add2, int alreadyPerfPrimeFacProd1, int alreadyPerfPrimeFacProd2, long int* invAddrBase1, long int* invAddrBase2) {
  Size1 = siz1;
  Size2 = siz2;
  localIndexCount = locInd;  
  xSSEObj = new xSSE;
  alreadyPerformedPrimeFacProducts1 = alreadyPerfPrimeFacProd1;
  alreadyPerformedPrimeFacProducts2 = alreadyPerfPrimeFacProd2;

  delete[] combinedPrefetchOffsetList;
  combinedPrefetchOffsetList = new int*[Size1];

  if (FFTcomplexFactorsForward1 != NULL) destroySuperAlignedComplex(FFTcomplexFactorsForward1);
  if (FFTcomplexFactorsBackward1 != NULL) destroySuperAlignedComplex(FFTcomplexFactorsBackward1);
  FFTcomplexFactorsForward1 = createSuperAlignedComplex(2);
  FFTcomplexFactorsBackward1 = createSuperAlignedComplex(2);

  double fac = -1;
  Complex x = Complex(cos(fac*pi/4), sin(fac*pi/4));
  FFTcomplexFactorsForward1[0].x = x.x;
  FFTcomplexFactorsForward1[0].y = x.x;
  FFTcomplexFactorsForward1[1].x = -x.y;
  FFTcomplexFactorsForward1[1].y = x.y;
  FFTcomplexFactorsBackward1[0].x = x.x;
  FFTcomplexFactorsBackward1[0].y = x.x;
  FFTcomplexFactorsBackward1[1].x = x.y;
  FFTcomplexFactorsBackward1[1].y = -x.y;
  
  ini();
}


void ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2::executeFFTstep(Complex* pIn, Complex* pOut, Complex* slice2D, bool forward) {
  Complex* FFTcomplexFacs1 = FFTcomplexFactorsForward1;
  if (!forward) {
    FFTcomplexFacs1 = FFTcomplexFactorsBackward1;
  }
  Complex facpmI = BitMask_ComplexForward;
  if (!forward) facpmI = BitMask_ComplexBackward;

  xSSEObj->xSSE_ExtremeFFT4D_DoFirstSecondAndThird1DD1FFTStepBase2_Wrapper(slice2D, Size1, Size2, localIndexCount, InternalEmbeddingValueOne, InternalEmbeddingValueTwo, FFTcomplexFacs1, facpmI, combinedPrefetchOffsetList);
}


void ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2::calcAddressesAndPrefetchBaseAddresses(long int* addresses, long int* prefetchBaseAddresses) {
  for (int I=0; I<localIndexCount*Size1*Size2; I++) prefetchBaseAddresses[I] = -1;
  for (int I=0; I<2*localIndexCount*Size1*Size2; I++) addresses[I] = -1;

  int addressCount = 0;  
  int timeIndex = 0;
  long int rsi = 0;
  int baseStep1 = Size2*localIndexCount*16 + InternalEmbeddingValueOne*64;
  int baseStep2 = 64*InternalEmbeddingValueTwo;
  for (int I=0; I<Size1/8; I+=1) {    
    for (int I2=0; I2<Size2*localIndexCount; I2+=1) {      
      prefetchBaseAddresses[timeIndex] = rsi;
      
      long int addr = rsi;
      for (int I3=0; I3<2; I3++) {
        for (int i=0; i<4; i++) {
          addresses[2*addressCount+0] = timeIndex;
	  addresses[2*addressCount+1] = addr;
	  addressCount++;
          addr += baseStep1;
        }
	addr += baseStep2;
      }
      rsi += 16;
      timeIndex++;	
    }
    rsi += 8*baseStep1 + 2*baseStep2 - 16*Size2*localIndexCount;
  }
}


void ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2::calcReadAndReadPrefetchBaseAddresses() {
  calcAddressesAndPrefetchBaseAddresses(readAddresses, readPrefetchBaseAddresses);
}


void ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2::calcWriteAndWritePrefetchBaseAddresses() {
  calcAddressesAndPrefetchBaseAddresses(writeAddresses, writePrefetchBaseAddresses);
}


int ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2::getNumberOfReadPrefetchesPerTimeIndex() {
  return 2;
}


int ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2::getNumberOfWritePrefetchesPerTimeIndex()  {
  return 2;
}


int ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2::getNumberOfReadTimeIndicesPerLine() {
  return localIndexCount*Size2;
} 


int ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2::getNumberOfWriteTimeIndicesPerLine() {
  return localIndexCount*Size2;
} 


void ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2::generateLocalPrefetchData() {
  for (int I=0; I<Size1; I++) {
    combinedPrefetchOffsetList[I] = NULL;
  }
  for (int I=0; I<Size1/8; I++) {
    combinedPrefetchOffsetList[I] = readPrefetchOffsetList[I];
  }  
}


int ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2::getNumberOfReadLines(){
  return Size1/8;
}


int ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2::getNumberOfWriteLines() {
  return Size1/8;
}


void ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2::changedInternalEmbedding() {
}
