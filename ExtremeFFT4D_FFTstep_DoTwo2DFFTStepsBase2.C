#include "ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2.h"

ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2::ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2() : ExtremeFFT4D_FFTstep() {
  readFromInput = false; 
  writeToOutput = false; 
  primeFacCount1 = 2; 
  primeFacCount2 = 2;
  primeFacs1[0] = 2;
  primeFacs1[1] = 2;
  primeFacs2[0] = 2;
  primeFacs2[1] = 2;
  mustBeFirstTrafo1 = false;
  mustBeFirstTrafo2 = false;        
  FFTcomplexFactorsForward1 = NULL;
  FFTcomplexFactorsForward2 = NULL;
  FFTcomplexFactorsBackward1 = NULL;
  FFTcomplexFactorsBackward2 = NULL;
  alreadyPerformedPrimeFacProducts1 = 0;
  alreadyPerformedPrimeFacProducts2 = 0;
  snprintf(name, 1000, "DoTwo2DFFTStepsBase2");
  combinedPrefetchOffsetList = NULL;
  useReadPrefetching = true;
  useWritePrefetching = false;
}


ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2::~ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2() {
  if (FFTcomplexFactorsForward1 != NULL) destroySuperAlignedComplex(FFTcomplexFactorsForward1);
  if (FFTcomplexFactorsForward2 != NULL) destroySuperAlignedComplex(FFTcomplexFactorsForward2);
  if (FFTcomplexFactorsBackward1 != NULL) destroySuperAlignedComplex(FFTcomplexFactorsBackward1);
  if (FFTcomplexFactorsBackward2 != NULL) destroySuperAlignedComplex(FFTcomplexFactorsBackward2);
  delete[] combinedPrefetchOffsetList;
}


ExtremeFFT4D_FFTstep* ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2::deriveSpecificFFTstepInstance(int siz1, int siz2, int locInd, int add1, int add2, int alreadyPerfPrimeFacProd1, int alreadyPerfPrimeFacProd2, long int* invAddrBase1, long int* invAddrBase2) {
  ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2* fftStep = new ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2();
  fftStep->setValuesInDerivedSpecificFFTstepInstance(siz1, siz2, locInd, add1, add2, alreadyPerfPrimeFacProd1, alreadyPerfPrimeFacProd2, invAddrBase1, invAddrBase2);
  return fftStep;
}


void ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2::setValuesInDerivedSpecificFFTstepInstance(int siz1, int siz2, int locInd, int add1, int add2, int alreadyPerfPrimeFacProd1, int alreadyPerfPrimeFacProd2, long int* invAddrBase1, long int* invAddrBase2) {
  Size1 = siz1;
  Size2 = siz2;
  localIndexCount = locInd;  
  xSSEObj = new xSSE;
  alreadyPerformedPrimeFacProducts1 = alreadyPerfPrimeFacProd1;
  alreadyPerformedPrimeFacProducts2 = alreadyPerfPrimeFacProd2;

  delete[] combinedPrefetchOffsetList;
  combinedPrefetchOffsetList = new int*[Size1];

  if (FFTcomplexFactorsForward1 != NULL) destroySuperAlignedComplex(FFTcomplexFactorsForward1);
  if (FFTcomplexFactorsForward2 != NULL) destroySuperAlignedComplex(FFTcomplexFactorsForward2);
  if (FFTcomplexFactorsBackward1 != NULL) destroySuperAlignedComplex(FFTcomplexFactorsBackward1);
  if (FFTcomplexFactorsBackward2 != NULL) destroySuperAlignedComplex(FFTcomplexFactorsBackward2);  
  FFTcomplexFactorsForward1 = createSuperAlignedComplex(2*3*alreadyPerfPrimeFacProd1);
  FFTcomplexFactorsForward2 = createSuperAlignedComplex(2*3*alreadyPerfPrimeFacProd2);
  FFTcomplexFactorsBackward1 = createSuperAlignedComplex(2*3*alreadyPerfPrimeFacProd1);
  FFTcomplexFactorsBackward2 = createSuperAlignedComplex(2*3*alreadyPerfPrimeFacProd2);
  
  double fac = -1;
  for (int I=0; I<alreadyPerfPrimeFacProd1; I++) {
    Complex x = Complex(cos(fac*(pi*I)/alreadyPerfPrimeFacProd1), sin(fac*(pi*I)/alreadyPerfPrimeFacProd1));
    FFTcomplexFactorsForward1[2*3*I+0].x = x.x;
    FFTcomplexFactorsForward1[2*3*I+0].y = x.x;
    FFTcomplexFactorsForward1[2*3*I+1].x = -x.y;
    FFTcomplexFactorsForward1[2*3*I+1].y = x.y;
    FFTcomplexFactorsBackward1[2*3*I+0].x = x.x;
    FFTcomplexFactorsBackward1[2*3*I+0].y = x.x;
    FFTcomplexFactorsBackward1[2*3*I+1].x = x.y;
    FFTcomplexFactorsBackward1[2*3*I+1].y = -x.y;
  
    x = Complex(cos(fac*0.5*(pi*I)/alreadyPerfPrimeFacProd1), sin(fac*0.5*(pi*I)/alreadyPerfPrimeFacProd1));
    FFTcomplexFactorsForward1[2*3*I+2].x = x.x;
    FFTcomplexFactorsForward1[2*3*I+2].y = x.x;
    FFTcomplexFactorsForward1[2*3*I+3].x = -x.y;
    FFTcomplexFactorsForward1[2*3*I+3].y = x.y;
    FFTcomplexFactorsBackward1[2*3*I+2].x = x.x;
    FFTcomplexFactorsBackward1[2*3*I+2].y = x.x;
    FFTcomplexFactorsBackward1[2*3*I+3].x = x.y;
    FFTcomplexFactorsBackward1[2*3*I+3].y = -x.y;
    
  
    x = Complex(cos(fac*(0.5*(pi*I)/alreadyPerfPrimeFacProd1+pi/2)), sin(fac*(0.5*(pi*I)/alreadyPerfPrimeFacProd1+pi/2)));
    FFTcomplexFactorsForward1[2*3*I+4].x = x.x;
    FFTcomplexFactorsForward1[2*3*I+4].y = x.x;
    FFTcomplexFactorsForward1[2*3*I+5].x = -x.y;
    FFTcomplexFactorsForward1[2*3*I+5].y = x.y;
    FFTcomplexFactorsBackward1[2*3*I+4].x = x.x;
    FFTcomplexFactorsBackward1[2*3*I+4].y = x.x;
    FFTcomplexFactorsBackward1[2*3*I+5].x = x.y;
    FFTcomplexFactorsBackward1[2*3*I+5].y = -x.y;
  }  
  
  for (int I=0; I<alreadyPerfPrimeFacProd2; I++) {
    Complex x = Complex(cos(fac*(pi*I)/alreadyPerfPrimeFacProd2), sin(fac*(pi*I)/alreadyPerfPrimeFacProd2));
    FFTcomplexFactorsForward2[2*3*I+0].x = x.x;
    FFTcomplexFactorsForward2[2*3*I+0].y = x.x;
    FFTcomplexFactorsForward2[2*3*I+1].x = -x.y;
    FFTcomplexFactorsForward2[2*3*I+1].y = x.y;
    FFTcomplexFactorsBackward2[2*3*I+0].x = x.x;
    FFTcomplexFactorsBackward2[2*3*I+0].y = x.x;
    FFTcomplexFactorsBackward2[2*3*I+1].x = x.y;
    FFTcomplexFactorsBackward2[2*3*I+1].y = -x.y;
  
    x = Complex(cos(fac*0.5*(pi*I)/alreadyPerfPrimeFacProd2), sin(fac*0.5*(pi*I)/alreadyPerfPrimeFacProd2));
    FFTcomplexFactorsForward2[2*3*I+2].x = x.x;
    FFTcomplexFactorsForward2[2*3*I+2].y = x.x;
    FFTcomplexFactorsForward2[2*3*I+3].x = -x.y;
    FFTcomplexFactorsForward2[2*3*I+3].y = x.y;
    FFTcomplexFactorsBackward2[2*3*I+2].x = x.x;
    FFTcomplexFactorsBackward2[2*3*I+2].y = x.x;
    FFTcomplexFactorsBackward2[2*3*I+3].x = x.y;
    FFTcomplexFactorsBackward2[2*3*I+3].y = -x.y;
    
  
    x = Complex(cos(fac*(0.5*(pi*I)/alreadyPerfPrimeFacProd2+pi/2)), sin(fac*(0.5*(pi*I)/alreadyPerfPrimeFacProd2+pi/2)));
    FFTcomplexFactorsForward2[2*3*I+4].x = x.x;
    FFTcomplexFactorsForward2[2*3*I+4].y = x.x;
    FFTcomplexFactorsForward2[2*3*I+5].x = -x.y;
    FFTcomplexFactorsForward2[2*3*I+5].y = x.y;
    FFTcomplexFactorsBackward2[2*3*I+4].x = x.x;
    FFTcomplexFactorsBackward2[2*3*I+4].y = x.x;
    FFTcomplexFactorsBackward2[2*3*I+5].x = x.y;
    FFTcomplexFactorsBackward2[2*3*I+5].y = -x.y;
  }    
  
  ini();
}


void ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2::executeFFTstep(Complex* pIn, Complex* pOut, Complex* slice2D, bool forward) {
  Complex* FFTcomplexFacs1 = FFTcomplexFactorsForward1;
  Complex* FFTcomplexFacs2 = FFTcomplexFactorsForward2;
  if (!forward) {
    FFTcomplexFacs1 = FFTcomplexFactorsBackward1;
    FFTcomplexFacs2 = FFTcomplexFactorsBackward2;
  }

  xSSEObj->xSSE_ExtremeFFT4D_DoTwo2DFFTStepsBase2_Wrapper(slice2D, Size1, Size2, localIndexCount, alreadyPerformedPrimeFacProducts1, alreadyPerformedPrimeFacProducts2, InternalEmbeddingValueOne, InternalEmbeddingValueTwo, FFTcomplexFacs1, FFTcomplexFacs2, combinedPrefetchOffsetList);
}


void ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2::calcAddressesAndPrefetchBaseAddresses(long int* addresses, long int* prefetchBaseAddresses) {
  for (int I=0; I<localIndexCount*Size1*Size2; I++) prefetchBaseAddresses[I] = -1;
  for (int I=0; I<2*localIndexCount*Size1*Size2; I++) addresses[I] = -1;

  int addressCount = 0;  
  int timeIndex = 0;
  int baseStep1 = alreadyPerformedPrimeFacProducts2*localIndexCount*16;
  int baseStep2 = localIndexCount*16;
  for (int I=0; I<Size1; I+=4*alreadyPerformedPrimeFacProducts1) {    
    for (int I4=0; I4<alreadyPerformedPrimeFacProducts1; I4++) {
      for (int I2=0; I2<Size2; I2+=4*alreadyPerformedPrimeFacProducts2) {      
        for (int I5=0; I5<alreadyPerformedPrimeFacProducts2; I5++) {
          for (int I3=0; I3<localIndexCount; I3++) {
	    long int totalLine = I+I4;
   	    long int addrBase = totalLine*(Size2*localIndexCount*16+InternalEmbeddingValueOne*64);
            addrBase += (totalLine/4)*64*InternalEmbeddingValueTwo;             
            prefetchBaseAddresses[timeIndex] = addrBase;
	    
            for (int i=0; i<4; i++) {
	      totalLine = I+I4+i*alreadyPerformedPrimeFacProducts1;
   	      addrBase = totalLine*(Size2*localIndexCount*16+InternalEmbeddingValueOne*64);
              addrBase += (totalLine/4)*64*InternalEmbeddingValueTwo;             
	      for (int i2=0; i2<4; i2++) {
	        long int addr = addrBase + i2*baseStep1 + (I2+I5)*baseStep2 + I3*16;
	        addresses[2*addressCount+0] = timeIndex;
	        addresses[2*addressCount+1] = addr;
	        addressCount++;
	      }
   	    }
  	    timeIndex++;	
	  }
	}
      }
    }
  }
}



void ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2::calcReadAndReadPrefetchBaseAddresses() {
  calcAddressesAndPrefetchBaseAddresses(readAddresses, readPrefetchBaseAddresses);
}


void ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2::calcWriteAndWritePrefetchBaseAddresses() {
  calcAddressesAndPrefetchBaseAddresses(writeAddresses, writePrefetchBaseAddresses);
}


int ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2::getNumberOfReadPrefetchesPerTimeIndex() {
  return 4;
}


int ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2::getNumberOfWritePrefetchesPerTimeIndex()  {
  return 4;
}


int ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2::getNumberOfReadTimeIndicesPerLine() {
  return (localIndexCount*Size2)/4;
} 


int ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2::getNumberOfWriteTimeIndicesPerLine() {
  return (localIndexCount*Size2)/4;
} 


void ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2::generateLocalPrefetchData() {
  for (int I=0; I<Size1; I++) {
    combinedPrefetchOffsetList[I] = NULL;
  }
  for (int I=0; I<Size1/4; I++) {
    combinedPrefetchOffsetList[I] = readPrefetchOffsetList[I];
  }  
}


int ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2::getNumberOfReadLines(){
  return Size1/4;
}


int ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2::getNumberOfWriteLines() {
  return Size1/4;
}


void ExtremeFFT4D_FFTstep_DoTwo2DFFTStepsBase2::changedInternalEmbedding() {
}
