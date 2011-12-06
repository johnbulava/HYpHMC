#include "ExtremeFFT4D_FFTstep_Write2DSlice.h"

ExtremeFFT4D_FFTstep_Write2DSlice::ExtremeFFT4D_FFTstep_Write2DSlice() : ExtremeFFT4D_FFTstep() {
  readFromInput = false; 
  writeToOutput = true; 
  primeFacCount1 = 0; 
  primeFacCount2 = 0;
  mustBeFirstTrafo1 = false;
  mustBeFirstTrafo2 = false;       
  snprintf(name, 1000, "Write2DSlice");
  useReadPrefetching = false;
  useWritePrefetching = false;
}


ExtremeFFT4D_FFTstep_Write2DSlice::~ExtremeFFT4D_FFTstep_Write2DSlice() {
}


ExtremeFFT4D_FFTstep* ExtremeFFT4D_FFTstep_Write2DSlice::deriveSpecificFFTstepInstance(int siz1, int siz2, int locInd, int add1, int add2, int alreadyPerfPrimeFacProd1, int alreadyPerfPrimeFacProd2, long int* invAddrBase1, long int* invAddrBase2) {
  ExtremeFFT4D_FFTstep_Write2DSlice* fftStep = new ExtremeFFT4D_FFTstep_Write2DSlice();
  fftStep->setValuesInDerivedSpecificFFTstepInstance(siz1, siz2, locInd, add1, add2, alreadyPerfPrimeFacProd1, alreadyPerfPrimeFacProd2, invAddrBase1, invAddrBase2);
  return fftStep;
}


void ExtremeFFT4D_FFTstep_Write2DSlice::setValuesInDerivedSpecificFFTstepInstance(int siz1, int siz2, int locInd, int add1, int add2, int alreadyPerfPrimeFacProd1, int alreadyPerfPrimeFacProd2, long int* invAddrBase1, long int* invAddrBase2) {
  Size1 = siz1;
  Size2 = siz2;
  localIndexCount = locInd;  
  xSSEObj = new xSSE();
  additionalIncrementsInBytes1 = add1 - Size2*add2;
  additionalIncrementsInBytes2 = add2 - 16*localIndexCount;
  
  ini();
}


void ExtremeFFT4D_FFTstep_Write2DSlice::executeFFTstep(Complex* pIn, Complex* pOut, Complex* slice2D, bool forward) {
  xSSEObj->xSSE_ExtremeFFT4D_Write2DSlice_Wrapper(slice2D, pOut, Size1, Size2, localIndexCount, additionalIncrementsInBytes1, additionalIncrementsInBytes2, InternalEmbeddingValueOne, InternalEmbeddingValueTwo);
}



void ExtremeFFT4D_FFTstep_Write2DSlice::calcReadAndReadPrefetchBaseAddresses() {
}


void ExtremeFFT4D_FFTstep_Write2DSlice::calcWriteAndWritePrefetchBaseAddresses() {
}


int ExtremeFFT4D_FFTstep_Write2DSlice::getNumberOfReadPrefetchesPerTimeIndex() {
  return 4;
}


int ExtremeFFT4D_FFTstep_Write2DSlice::getNumberOfWritePrefetchesPerTimeIndex()  {
  return 4;
}


int ExtremeFFT4D_FFTstep_Write2DSlice::getNumberOfReadTimeIndicesPerLine() {
  return (localIndexCount*Size2)/4;
} 


int ExtremeFFT4D_FFTstep_Write2DSlice::getNumberOfWriteTimeIndicesPerLine() {
  return (localIndexCount*Size2)/4;
} 


void ExtremeFFT4D_FFTstep_Write2DSlice::generateLocalPrefetchData() {
}


int ExtremeFFT4D_FFTstep_Write2DSlice::getNumberOfReadLines(){
  return Size1/4;
}


int ExtremeFFT4D_FFTstep_Write2DSlice::getNumberOfWriteLines() {
  return Size1/4;
}


void ExtremeFFT4D_FFTstep_Write2DSlice::changedInternalEmbedding() {
}
