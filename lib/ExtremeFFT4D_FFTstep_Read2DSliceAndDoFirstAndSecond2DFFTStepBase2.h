#ifndef ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_included
#define ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_included

#include <math.h>
#include "Global.h"
#include "Complex.h"
#include "xSSE.h"
#include "ExtremeFFT4D_FFTstep.h"



class ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2 : public ExtremeFFT4D_FFTstep {
private:
  long int* AddressAuxData1;
  long int* AddressAuxData2;  
  int addressIncrementsInBytes1;
  int addressIncrementsInBytes2;
  Complex BitMask_ComplexForward;
  Complex BitMask_ComplexBackward;
  int** combinedPrefetchOffsetList;

  void calcReadAndReadPrefetchBaseAddresses();
  void calcWriteAndWritePrefetchBaseAddresses();
  int getNumberOfReadPrefetchesPerTimeIndex();
  int getNumberOfWritePrefetchesPerTimeIndex();
  int getNumberOfReadTimeIndicesPerLine();
  int getNumberOfWriteTimeIndicesPerLine();
  int getNumberOfReadLines();
  int getNumberOfWriteLines();  
  void generateLocalPrefetchData();
  void changedInternalEmbedding();


public:
  ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2(); 
  ~ExtremeFFT4D_FFTstep_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2(); 
  
  ExtremeFFT4D_FFTstep* deriveSpecificFFTstepInstance(int siz1, int siz2, int locInd, int add1, int add2, int alreadyPerfPrimeFacProd1, int alreadyPerfPrimeFacProd2, long int* invAddrBase1, long int* invAddrBase2); 
  void setValuesInDerivedSpecificFFTstepInstance(int siz1, int siz2, int locInd, int add1, int add2, int alreadyPerfPrimeFacProd1, int alreadyPerfPrimeFacProd2, long int* invAddrBase1, long int* invAddrBase2); 
  void executeFFTstep(Complex* pIn, Complex* pOut, Complex* slice2D, bool forward);
  
};

#endif
