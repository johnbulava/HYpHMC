#ifndef ExtremeFFT4D_FFTstep_included
#define ExtremeFFT4D_FFTstep_included

#include <math.h>
#include "Global.h"
#include "Complex.h"
#include "xSSE.h"
#include "Tools.h"
#include "CacheSimulator.h"

#define ExtremeFFT4D_FFTstep_PrimeFactorsMAX 10


class ExtremeFFT4D_FFTstep {
protected:
  bool readFromInput; 
  bool writeToOutput; 
  int primeFacCount1; 
  int primeFacCount2;
  int primeFacs1[ExtremeFFT4D_FFTstep_PrimeFactorsMAX];
  int primeFacs2[ExtremeFFT4D_FFTstep_PrimeFactorsMAX];
  bool mustBeFirstTrafo1;
  bool mustBeFirstTrafo2;  
  bool readyForExecution;
  xSSE* xSSEObj;
  int Size1;
  int Size2;
  int localIndexCount;  
  char* name;
  CacheSimulator* L1Cache;
  
  long int* readPrefetchBaseAddresses;     // Base-Addresses for prefetching for given Time-Index
  long int* writePrefetchBaseAddresses;
  long int* readAddresses;                 // Time-Index + Address
  long int* writeAddresses;
  long int* readPrefetchAddresses;         // Time-Index + Address
  long int* writePrefetchAddresses;        
  long int* prefetchDataGenerationHelperData1;
  long int* prefetchDataGenerationHelperData2;  
  int* prefetchDataSpaceHolder;
  int** readPrefetchOffsets;
  int** writePrefetchOffsets;
  int** readPrefetchOffsetList;
  int** writePrefetchOffsetList;
  bool useReadPrefetching;
  bool useWritePrefetching;
  
  
  
  
  int readPrefetchCombineFac;
  int writePrefetchCombineFac;
  int InternalEmbeddingValueOne;
  int InternalEmbeddingValueTwo;
  
  
  virtual void calcReadAndReadPrefetchBaseAddresses() = 0;
  virtual void calcWriteAndWritePrefetchBaseAddresses() = 0;
  virtual int getNumberOfReadPrefetchesPerTimeIndex() = 0;
  virtual int getNumberOfWritePrefetchesPerTimeIndex() = 0;
  virtual int getNumberOfReadTimeIndicesPerLine() = 0;
  virtual int getNumberOfWriteTimeIndicesPerLine() = 0;
  virtual int getNumberOfReadLines() = 0;
  virtual int getNumberOfWriteLines() = 0;  
  virtual void generateLocalPrefetchData() = 0;
  virtual void changedInternalEmbedding() = 0;
  
  
  void ini();
  void desini();
  void generatePrefetchData(int PrefetchCombineFac, int loadsPerTimeIndex, int TimeIndicesPerLine, int lines, long int* touchedAddresses, long int* prefetchAddresses, long int* prefetchBaseAddresses, int** prefetchOffsets, int** prefetchOffsetsList);
  void generatePrefetchData();
  

public:
  ExtremeFFT4D_FFTstep(); 
  virtual ~ExtremeFFT4D_FFTstep(); 
  
  char* getName();
  bool getReadFromInput(); 
  bool getWriteToOutput(); 
  int getPrimeFacCount1(); 
  int getPrimeFacCount2();
  int* getPrimeFacs1();
  int* getPrimeFacs2();
  bool getMustBeFirstTrafo1();
  bool getMustBeFirstTrafo2();       
  bool getReadyForExecution();

  void setInternalEmbedding(int intEmbeddOne, int intEmbeddTwo);
  void setPrefetchCombineFacs(int readComF, int writeComF);
  void calcL1ExcessMisses(long int sizeOfWhole4DFieldInBytes, long int inputSlice2DDistanceInBytes, int &minorMisses, int &majorMisses, bool abortWhenResultWorse);
  


  virtual ExtremeFFT4D_FFTstep* deriveSpecificFFTstepInstance(int siz1, int siz2, int locInd, int add1, int add2, int alreadyPerfPrimeFacProd1, int alreadyPerfPrimeFacProd2, long int* invAddrBase1, long int* invAddrBase2) = 0; 
  virtual void setValuesInDerivedSpecificFFTstepInstance(int siz1, int siz2, int locInd, int add1, int add2, int alreadyPerfPrimeFacProd1, int alreadyPerfPrimeFacProd2, long int* invAddrBase1, long int* invAddrBase2) = 0; 
  virtual void executeFFTstep(Complex* pIn, Complex* pOut, Complex* slice2D, bool forward) = 0;
  
};

#endif
