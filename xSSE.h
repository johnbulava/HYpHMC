#ifndef xSSE_included
#define xSSE_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>

#include "Global.h"
#include "Complex.h"

#define xSSE_ParameterSpace_Size 1600
#define xSSE_Level1Cache_Size 65536
#define xSSE_WorkSpace_Size 65536



extern "C" {
  void xSSE_ComplexCopy(char* parameters);
  void xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2(char* parameters);  
  void xSSE_ExtremeFFT4D_DoTwo2DFFTStepsBase2(char* parameters);
  void xSSE_ExtremeFFT4D_Write2DSlice(char* parameters);
  void xSSE_ReadLongIntFromMemAddr(char* parameters);
  void xSSE_WriteLongIntToMemAddr(char* parameters);
  void xSSE_PerformNopLoop(char* parameters);
  void xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2(char* parameters);
  void xSSE_ExtremeFFT4D_DoFirstSecondAndThird1DD1FFTStepBase2(char* parameters);
  void xSSE_ComplexVectorAddition(char* parameters);  
}


class xSSE {
private:
  char* xSSE_ParameterSpace;
  char* xSSE_WorkSpace;
  //Pointers to Parameter-Space
  char* xSSE_ParameterSpace_Pointer;
  char** xSSE_ParameterSpace_Pointer_CharPointer;
  Complex* xSSE_ParameterSpace_Pointer_Complex;
  double* xSSE_ParameterSpace_Pointer_Double;
  long int* xSSE_ParameterSpace_Pointer_LongInt;
  //Pointer to WorkSpace
  char* xSSE_WorkSpace_Pointer;


public:
  xSSE();
  ~xSSE();

  void xSSE_ComplexCopy_Wrappper(Complex* input, Complex* output, int blockSize, int L0, int L1, int L2, int L3, int xtrSize1, int xtrSize2, int xtrSize3);
  void xSSE_ReadLongIntFromMemAddr_Wrapper(long int* memAddr, long int& value);
  void xSSE_WriteLongIntToMemAddr_Wrapper(long int* memAddr, long int value);
  void xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_Wrapper(Complex* input, Complex* slice, int Size1, int Size2, int localIndexCount, long int* bitInvertedInputAddressAuxData1, long int* bitInvertedOutputAddressAuxData2, int** prefetchAddressControlList, int add2, int embeddingOne, int embeddingTwo, Complex facpmI);
  void xSSE_ExtremeFFT4D_DoTwo2DFFTStepsBase2_Wrapper(Complex* slice, int Size1, int Size2, int localIndexCount, int smallerHalfBlockSize1, int smallerHalfBlockSize2, int embeddingOne, int embeddingTwo, Complex* FFTcomplexFacs1, Complex* FFTcomplexFacs2, int** prefetchAddressControlList);
  void xSSE_ExtremeFFT4D_Write2DSlice_Wrapper(Complex* slice, Complex* output, int Size1, int Size2, int localIndexCount, int add1, int add2, int embeddingOne, int embeddingTwo);
  void xSSE_PerformNopLoop_Wrapper(long int count);  
  void xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2_Wrapper(Complex* input, Complex* slice, int Size1, int Size2, int localIndexCount, long int* bitInvertedOutputAddressAuxData1, long int* bitInvertedOutputAddressAuxData2, int** prefetchAddressControlList, int add1, int add2, Complex* FFTcomplexFacs, Complex facpmI);
  void xSSE_ExtremeFFT4D_DoFirstSecondAndThird1DD1FFTStepBase2_Wrapper(Complex* slice, int Size1, int Size2, int localIndexCount, int embeddingOne, int embeddingTwo, Complex* FFTcomplexFacs1, Complex facpmI, int** prefetchAddressControlList);
  void xSSE_ComplexVectorAddition_Wrapper(Complex* input, Complex* output, Complex alpha, int blockSize, int L0, int L1, int L2, int L3, int xtrSize1, int xtrSize2, int xtrSize3);  
};

#endif
