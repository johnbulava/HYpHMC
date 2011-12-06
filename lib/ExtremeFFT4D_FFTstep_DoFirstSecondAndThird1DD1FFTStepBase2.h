#ifndef ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2_included
#define ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2_included

#include <math.h>
#include "Global.h"
#include "Complex.h"
#include "xSSE.h"
#include "ExtremeFFT4D_FFTstep.h"


class ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2 : public ExtremeFFT4D_FFTstep {
private:
  Complex* FFTcomplexFactorsForward1;
  Complex* FFTcomplexFactorsBackward1;
  Complex BitMask_ComplexForward;
  Complex BitMask_ComplexBackward;
  int alreadyPerformedPrimeFacProducts1;
  int alreadyPerformedPrimeFacProducts2;
  int** combinedPrefetchOffsetList;
  
  void calcAddressesAndPrefetchBaseAddresses(long int* addresses, long int* prefetchBaseAddresses);
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
  ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2(); 
  ~ExtremeFFT4D_FFTstep_DoFirstSecondAndThird1DD1FFTStepBase2(); 
  
  ExtremeFFT4D_FFTstep* deriveSpecificFFTstepInstance(int siz1, int siz2, int locInd, int add1, int add2, int alreadyPerfPrimeFacProd1, int alreadyPerfPrimeFacProd2, long int* invAddrBase1, long int* invAddrBase2); 
  void setValuesInDerivedSpecificFFTstepInstance(int siz1, int siz2, int locInd, int add1, int add2, int alreadyPerfPrimeFacProd1, int alreadyPerfPrimeFacProd2, long int* invAddrBase1, long int* invAddrBase2); 
  void executeFFTstep(Complex* pIn, Complex* pOut, Complex* slice2D, bool forward);
  
};

#endif
