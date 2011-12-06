#include "xSSE.h"

xSSE::xSSE() {
  xSSE_ParameterSpace = new char[xSSE_ParameterSpace_Size+128];  
  xSSE_WorkSpace = new char[xSSE_WorkSpace_Size+128];
  for (int I=0; I<xSSE_ParameterSpace_Size+128; I++) xSSE_ParameterSpace[I] = 0;
  for (int I=0; I<xSSE_WorkSpace_Size+128; I++) xSSE_WorkSpace[I] = 0;

  //Pointers to Parameter-Space
  long int offset1 = 128 - (((long int) xSSE_ParameterSpace) % 128);
  xSSE_ParameterSpace_Pointer = (char*) &(xSSE_ParameterSpace[offset1]);
  xSSE_ParameterSpace_Pointer_CharPointer = (char**) xSSE_ParameterSpace_Pointer;
  xSSE_ParameterSpace_Pointer_Complex = (Complex*) xSSE_ParameterSpace_Pointer;
  xSSE_ParameterSpace_Pointer_Double = (double*) xSSE_ParameterSpace_Pointer;
  xSSE_ParameterSpace_Pointer_LongInt = (long int*) xSSE_ParameterSpace_Pointer;
  //Pointer to Work.Space
  long int offset2 = 128 - (((long int) xSSE_WorkSpace) % 128);
  xSSE_WorkSpace_Pointer = (char*) &(xSSE_WorkSpace[offset2]);
  
  //Alignment - Check
  offset1 = (long int) xSSE_ParameterSpace_Pointer;
  offset2 = (long int) xSSE_WorkSpace_Pointer;
  if ((offset1 % 128) != 0) {
    printf("ERROR in xSSE-Constructor: xSSE_ParameterSpace not properly aligned!\n");
    exit(0);  
  }
  if ((offset2 % 128) != 0) {
    printf("ERROR in xSSE-Constructor: xSSE_WorkSpace_Pointer not properly aligned!\n");
    exit(0);  
  }
}


xSSE::~xSSE() {
  delete[] xSSE_ParameterSpace;
  delete[] xSSE_WorkSpace;
}


void xSSE::xSSE_ComplexCopy_Wrappper(Complex* input, Complex* output, int blockSize, int L0, int L1, int L2, int L3, int xtrSize1, int xtrSize2, int xtrSize3) {
  xSSE_ParameterSpace_Pointer_CharPointer[0] = xSSE_WorkSpace_Pointer;
  xSSE_ParameterSpace_Pointer_CharPointer[1] = (char*) input;
  xSSE_ParameterSpace_Pointer_CharPointer[2] = (char*) output;
  
  int L3LineSize = (L3+xtraSize3)*blockSize*16;
  int L3CombineCount = (xSSE_Level1Cache_Size / 2) / L3LineSize;
  if (L3CombineCount>L2) L3CombineCount = L2;
  
  xSSE_ParameterSpace_Pointer_LongInt[3]  = blockSize;
  xSSE_ParameterSpace_Pointer_LongInt[4]  = L0;
  xSSE_ParameterSpace_Pointer_LongInt[5]  = L1;
  xSSE_ParameterSpace_Pointer_LongInt[6]  = L2;
  xSSE_ParameterSpace_Pointer_LongInt[7]  = L3;
  xSSE_ParameterSpace_Pointer_LongInt[8]  = 16*xtraSize1*blockSize*(L3+xtraSize3)*(L2+xtraSize2);
  xSSE_ParameterSpace_Pointer_LongInt[9]  = 16*xtraSize2*blockSize*(L3+xtraSize3);
  xSSE_ParameterSpace_Pointer_LongInt[10] = 16*xtraSize3*blockSize;
  xSSE_ParameterSpace_Pointer_LongInt[11] = L3CombineCount;
  xSSE_ParameterSpace_Pointer_LongInt[12] = 4;
  xSSE_ParameterSpace_Pointer_LongInt[13] = -1;
  
  xSSE_ComplexCopy(xSSE_ParameterSpace_Pointer);
  
  if (xSSE_ParameterSpace_Pointer_LongInt[13] != 0) {
    printf("ERROR in xSSE_ComplexCopy (error code: %ld)!!!\n", xSSE_ParameterSpace_Pointer_LongInt[13]);
    exit(0);
  }
}


void xSSE::xSSE_PerformNopLoop_Wrapper(long int count) {
  xSSE_ParameterSpace_Pointer_LongInt[0]  = count;
  xSSE_ParameterSpace_Pointer_LongInt[1]  = -1;
  
  xSSE_PerformNopLoop(xSSE_ParameterSpace_Pointer);
  
  if (xSSE_ParameterSpace_Pointer_LongInt[1] != 0) {
    printf("ERROR in xSSE_PerformNopLoop (error code: %ld)!!!\n", xSSE_ParameterSpace_Pointer_LongInt[1]);
    exit(0);
  }
}


void xSSE::xSSE_ReadLongIntFromMemAddr_Wrapper(long int* memAddr, long int& value) {
  xSSE_ParameterSpace_Pointer_CharPointer[0] = (char*) memAddr;
  xSSE_ParameterSpace_Pointer_LongInt[1]  = value;
  xSSE_ParameterSpace_Pointer_LongInt[2]  = -1;
  
  xSSE_ReadLongIntFromMemAddr(xSSE_ParameterSpace_Pointer);
  
  value = xSSE_ParameterSpace_Pointer_LongInt[1];
  
  if (xSSE_ParameterSpace_Pointer_LongInt[2] != 0) {
    printf("ERROR in xSSE_ReadLongIntFromMemAddr (error code: %ld)!!!\n", xSSE_ParameterSpace_Pointer_LongInt[2]);
    exit(0);
  }
}


void xSSE::xSSE_WriteLongIntToMemAddr_Wrapper(long int* memAddr, long int value) {
  xSSE_ParameterSpace_Pointer_CharPointer[0] = (char*) memAddr;
  xSSE_ParameterSpace_Pointer_LongInt[1]  = value;
  xSSE_ParameterSpace_Pointer_LongInt[2]  = -1;
  
  xSSE_WriteLongIntToMemAddr(xSSE_ParameterSpace_Pointer);
  
  if (xSSE_ParameterSpace_Pointer_LongInt[2] != 0) {
    printf("ERROR in xSSE_WriteLongIntToMemAddr (error code: %ld)!!!\n", xSSE_ParameterSpace_Pointer_LongInt[2]);
    exit(0);
  }
}

void xSSE::xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_Wrapper(Complex* input, Complex* slice, int Size1, int Size2, int localIndexCount, long int* bitInvertedInputAddressAuxData1, long int* bitInvertedOutputAddressAuxData2, int** prefetchAddressControlList, int add2, int embeddingOne, int embeddingTwo, Complex facpmI) {
  xSSE_ParameterSpace_Pointer_CharPointer[0] = xSSE_WorkSpace_Pointer;
  xSSE_ParameterSpace_Pointer_CharPointer[1] = (char*) input;
  xSSE_ParameterSpace_Pointer_CharPointer[2] = (char*) slice;
  
  xSSE_ParameterSpace_Pointer_LongInt[3]  = Size1;
  xSSE_ParameterSpace_Pointer_LongInt[4]  = Size2;
  xSSE_ParameterSpace_Pointer_LongInt[5]  = localIndexCount;
  xSSE_ParameterSpace_Pointer_CharPointer[6] = (char*) bitInvertedInputAddressAuxData1;
  xSSE_ParameterSpace_Pointer_CharPointer[7] = (char*) bitInvertedOutputAddressAuxData2;  
  xSSE_ParameterSpace_Pointer_CharPointer[8] = (char*) prefetchAddressControlList;  
  xSSE_ParameterSpace_Pointer_LongInt[9]  = add2;
  xSSE_ParameterSpace_Pointer_LongInt[10] = embeddingOne;
  xSSE_ParameterSpace_Pointer_LongInt[11] = embeddingTwo;
  xSSE_ParameterSpace_Pointer_Complex[6] = facpmI;
  xSSE_ParameterSpace_Pointer_Complex[7].x = 2.0;
  xSSE_ParameterSpace_Pointer_Complex[7].y = 2.0;
  xSSE_ParameterSpace_Pointer_LongInt[16] = -1;  //Error-Code
  
  xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2(xSSE_ParameterSpace_Pointer);
  
  if (xSSE_ParameterSpace_Pointer_LongInt[16] != 0) {
    printf("ERROR in xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2 (error code: %ld)!!!\n", xSSE_ParameterSpace_Pointer_LongInt[16]);
    exit(0);
  }
}


void xSSE::xSSE_ExtremeFFT4D_DoTwo2DFFTStepsBase2_Wrapper(Complex* slice, int Size1, int Size2, int localIndexCount, int smallerHalfBlockSize1, int smallerHalfBlockSize2, int embeddingOne, int embeddingTwo, Complex* FFTcomplexFacs1, Complex* FFTcomplexFacs2, int** prefetchAddressControlList) {
  xSSE_ParameterSpace_Pointer_CharPointer[0] = xSSE_WorkSpace_Pointer;
  xSSE_ParameterSpace_Pointer_CharPointer[1] = (char*) slice;
  
  xSSE_ParameterSpace_Pointer_LongInt[2]  = Size1;
  xSSE_ParameterSpace_Pointer_LongInt[3]  = Size2;
  xSSE_ParameterSpace_Pointer_LongInt[4]  = localIndexCount;
  xSSE_ParameterSpace_Pointer_LongInt[5]  = smallerHalfBlockSize1;
  xSSE_ParameterSpace_Pointer_LongInt[6]  = smallerHalfBlockSize2;
  xSSE_ParameterSpace_Pointer_LongInt[7]  = embeddingOne;  
  xSSE_ParameterSpace_Pointer_LongInt[8]  = embeddingTwo;
  xSSE_ParameterSpace_Pointer_CharPointer[9] = (char*) FFTcomplexFacs1;
  xSSE_ParameterSpace_Pointer_CharPointer[10] = (char*) FFTcomplexFacs2;

  xSSE_ParameterSpace_Pointer_CharPointer[11] = (char*) prefetchAddressControlList;  
  xSSE_ParameterSpace_Pointer_Complex[6].x = 2.0;
  xSSE_ParameterSpace_Pointer_Complex[6].y = 2.0;
  xSSE_ParameterSpace_Pointer_LongInt[14] = -1;  //Error-Code
  
  xSSE_ExtremeFFT4D_DoTwo2DFFTStepsBase2(xSSE_ParameterSpace_Pointer);
  
  if (xSSE_ParameterSpace_Pointer_LongInt[14] != 0) {
    printf("ERROR in xSSE_ExtremeFFT4D_DoTwo2DFFTStepsBase2 (error code: %ld)!!!\n", xSSE_ParameterSpace_Pointer_LongInt[14]);
    exit(0);
  }
}


void xSSE::xSSE_ExtremeFFT4D_Write2DSlice_Wrapper(Complex* slice, Complex* output, int Size1, int Size2, int localIndexCount, int add1, int add2, int embeddingOne, int embeddingTwo) {
  xSSE_ParameterSpace_Pointer_CharPointer[0] = xSSE_WorkSpace_Pointer;
  xSSE_ParameterSpace_Pointer_CharPointer[1] = (char*) slice;
  xSSE_ParameterSpace_Pointer_CharPointer[2] = (char*) output;
  
  xSSE_ParameterSpace_Pointer_LongInt[3]  = Size1;
  xSSE_ParameterSpace_Pointer_LongInt[4]  = Size2;
  xSSE_ParameterSpace_Pointer_LongInt[5]  = localIndexCount;
  xSSE_ParameterSpace_Pointer_LongInt[6]  = add1;
  xSSE_ParameterSpace_Pointer_LongInt[7]  = add2;  
  xSSE_ParameterSpace_Pointer_LongInt[8]  = embeddingOne;
  xSSE_ParameterSpace_Pointer_LongInt[9]  = embeddingTwo;  
  xSSE_ParameterSpace_Pointer_LongInt[10] = -1;  //Error-Code
  
  xSSE_ExtremeFFT4D_Write2DSlice(xSSE_ParameterSpace_Pointer);
  
  if (xSSE_ParameterSpace_Pointer_LongInt[10] != 0) {
    printf("ERROR in xSSE_ExtremeFFT4D_Write2DSlice (error code: %ld)!!!\n", xSSE_ParameterSpace_Pointer_LongInt[10]);
    exit(0);
  }
}


void xSSE::xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2_Wrapper(Complex* input, Complex* slice, int Size1, int Size2, int localIndexCount, long int* bitInvertedOutputAddressAuxData1, long int* bitInvertedOutputAddressAuxData2, int** prefetchAddressControlList, int add1, int add2, Complex* FFTcomplexFacs, Complex facpmI) {
  xSSE_ParameterSpace_Pointer_CharPointer[0] = xSSE_WorkSpace_Pointer;
  xSSE_ParameterSpace_Pointer_CharPointer[1] = (char*) input;
  xSSE_ParameterSpace_Pointer_CharPointer[2] = (char*) slice;
  
  xSSE_ParameterSpace_Pointer_LongInt[3]  = Size1;
  xSSE_ParameterSpace_Pointer_LongInt[4]  = Size2; 
  xSSE_ParameterSpace_Pointer_LongInt[5]  = localIndexCount;
  xSSE_ParameterSpace_Pointer_CharPointer[6] = (char*) bitInvertedOutputAddressAuxData1;
  xSSE_ParameterSpace_Pointer_CharPointer[7] = (char*) bitInvertedOutputAddressAuxData2;  
  xSSE_ParameterSpace_Pointer_CharPointer[8] = (char*) prefetchAddressControlList;  
  xSSE_ParameterSpace_Pointer_LongInt[9]  = add1;
  xSSE_ParameterSpace_Pointer_LongInt[10]  = add2;
  xSSE_ParameterSpace_Pointer_CharPointer[11] = (char*) FFTcomplexFacs;
  xSSE_ParameterSpace_Pointer_Complex[6] = facpmI;
  xSSE_ParameterSpace_Pointer_LongInt[14] = -1;  //Error-Code
  
  xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2(xSSE_ParameterSpace_Pointer);
  
  if (xSSE_ParameterSpace_Pointer_LongInt[14] != 0) {
    printf("ERROR in xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2 (error code: %ld)!!!\n", xSSE_ParameterSpace_Pointer_LongInt[14]);
    exit(0);
  }
}


void xSSE::xSSE_ExtremeFFT4D_DoFirstSecondAndThird1DD1FFTStepBase2_Wrapper(Complex* slice, int Size1, int Size2, int localIndexCount, int embeddingOne, int embeddingTwo, Complex* FFTcomplexFacs1, Complex facpmI, int** prefetchAddressControlList) {
  xSSE_ParameterSpace_Pointer_CharPointer[0] = xSSE_WorkSpace_Pointer;
  xSSE_ParameterSpace_Pointer_CharPointer[1] = (char*) slice;
  
  xSSE_ParameterSpace_Pointer_LongInt[2]  = Size1;
  xSSE_ParameterSpace_Pointer_LongInt[3]  = Size2;
  xSSE_ParameterSpace_Pointer_LongInt[4]  = localIndexCount;
  xSSE_ParameterSpace_Pointer_LongInt[5]  = embeddingOne;  
  xSSE_ParameterSpace_Pointer_LongInt[6]  = embeddingTwo;
  xSSE_ParameterSpace_Pointer_CharPointer[7] = (char*) FFTcomplexFacs1;
  xSSE_ParameterSpace_Pointer_Complex[4] = facpmI;
  xSSE_ParameterSpace_Pointer_CharPointer[10] = (char*) prefetchAddressControlList;  
  xSSE_ParameterSpace_Pointer_LongInt[15] = -1;  //Error-Code
  
  xSSE_ExtremeFFT4D_DoFirstSecondAndThird1DD1FFTStepBase2(xSSE_ParameterSpace_Pointer);
  
  if (xSSE_ParameterSpace_Pointer_LongInt[15] != 0) {
    printf("ERROR in xSSE_ExtremeFFT4D_DoFirstSecondAndThird1DD1FFTStepBase2 (error code: %ld)!!!\n", xSSE_ParameterSpace_Pointer_LongInt[15]);
    exit(0);
  }
}


void xSSE::xSSE_ComplexVectorAddition_Wrapper(Complex* input, Complex* output, Complex alpha, int blockSize, int L0, int L1, int L2, int L3, int xtrSize1, int xtrSize2, int xtrSize3) {
  xSSE_ParameterSpace_Pointer_CharPointer[0] = xSSE_WorkSpace_Pointer;
  xSSE_ParameterSpace_Pointer_CharPointer[1] = (char*) input;
  xSSE_ParameterSpace_Pointer_CharPointer[2] = (char*) output;
  
  int L3LineSize = (L3+xtraSize3)*blockSize*16;
  int L3CombineCount = (xSSE_Level1Cache_Size / 2) / L3LineSize;
  if (L3CombineCount>L2) L3CombineCount = L2;

  xSSE_ParameterSpace_Pointer_LongInt[3]  = blockSize;
  xSSE_ParameterSpace_Pointer_Complex[2] = alpha;
  xSSE_ParameterSpace_Pointer_LongInt[6]  = L0;
  xSSE_ParameterSpace_Pointer_LongInt[7]  = L1;
  xSSE_ParameterSpace_Pointer_LongInt[8]  = L2;
  xSSE_ParameterSpace_Pointer_LongInt[9]  = L3;
  xSSE_ParameterSpace_Pointer_LongInt[10]  = 16*xtraSize1*blockSize*(L3+xtraSize3)*(L2+xtraSize2);
  xSSE_ParameterSpace_Pointer_LongInt[11]  = 16*xtraSize2*blockSize*(L3+xtraSize3);
  xSSE_ParameterSpace_Pointer_LongInt[12] = 16*xtraSize3*blockSize;
  xSSE_ParameterSpace_Pointer_LongInt[13] = L3CombineCount;
  xSSE_ParameterSpace_Pointer_LongInt[14] = 4;
  xSSE_ParameterSpace_Pointer_LongInt[15] = xSSE_Level1Cache_Size / 2;  
  xSSE_ParameterSpace_Pointer_LongInt[16] = -1;
  
  xSSE_ComplexVectorAddition(xSSE_ParameterSpace_Pointer);

  if (xSSE_ParameterSpace_Pointer_LongInt[16] != 0) {
    printf("ERROR in xSSE_ComplexVectorAddition (error code: %ld)!!!\n", xSSE_ParameterSpace_Pointer_LongInt[16]);
    exit(0);
  }
}
