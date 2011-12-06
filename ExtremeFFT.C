#include "ExtremeFFT.h"
#include "SSEroutines.h"

ExtremeFFT::ExtremeFFT() {
  FourierTrafo_bitInverse = NULL;
  FourierTrafo_iniForL0 = -1;
  FourierTrafo_iniForL1 = -1;
  FourierTrafo_iniForL2 = -1;
  FourierTrafo_iniForL3 = -1;
  FourierTrafo_iniForN = -1;
  FourierTrafo_Weights = NULL;
  CACHE_L1 = NULL;
  CACHE_L2 = NULL;
  BUFFER = NULL;
}


void ExtremeFFT::SlowRearrange(int oneDimL0, int oneDimL1, int oneDimL2, int oneDimL3, Complex* input, Complex* output) {
  int I1, I2, I3, I4, i;
  int index1, index2;
  int o1 = oneDimL3;
  int o2 = oneDimL2 * o1;
  int o3 = oneDimL1 * o2;
  int l3 = (int)((log(1.0*oneDimL3) / log(2.0))+0.001);
  int l2 = (int)((log(1.0*oneDimL2) / log(2.0))+0.001);
  int l1 = (int)((log(1.0*oneDimL1) / log(2.0))+0.001);
  int l0 = (int)((log(1.0*oneDimL0) / log(2.0))+0.001);
  
  for (I1=0; I1<oneDimL0; I1++) {
    for (I2=0; I2<oneDimL1; I2++) {
      for (I3=0; I3<oneDimL2; I3++) {
        for (I4=0; I4<oneDimL3; I4++) {
          index1 = I1*o3 + I2*o2 + I3*o1 + I4;
          index2 = bitInverter(I1,l0)*o3 + bitInverter(I2,l1)*o2 + bitInverter(I3,l2)*o1 + bitInverter(I4,l3);
          for (i=0; i<8; i++) {
	    output[8*index2+i] = input[8*index1+i];
	  }
	}
      }
    }
  }
}



void ExtremeFFT::SlowFourierStep(int oneDimL0, int oneDimL1, int oneDimL2, int oneDimL3, int stepNr, Complex* input, Complex* output) {
/*  int O1,O2,O3,O4, I1,I2,I3,I4, F1,F2,F3,F4, i;
  int innerL = powINT(2,stepNr+1);
  int outerL = oneDim / innerL;
  int fmax = innerL / 2;
  double w = -2*pi / innerL;
  int offO1 = oneDim / outerL;
  int offO2 = oneDim*oneDim / outerL;
  int offO3 = oneDim*oneDim*oneDim / outerL;
  int offO4 = oneDim*oneDim*oneDim*oneDim / outerL;
  int index;
  
  for (O1=0; O1<outerL; O1++) {
  for (O2=0; O2<outerL; O2++) {
  for (O3=0; O3<outerL; O3++) {
  for (O4=0; O4<outerL; O4++) {
  
    for (I1=0; I1<innerL; I1++) {
    for (I2=0; I2<innerL; I2++) {
    for (I3=0; I3<innerL; I3++) {
    for (I4=0; I4<innerL; I4++) {
  
      for (i=0; i<8; i++) {
        Complex res(0,0);
        for (F1=0; F1<2; F1++) {
        for (F2=0; F2<2; F2++) {
        for (F3=0; F3<2; F3++) {
        for (F4=0; F4<2; F4++) {
	  int fac = I1*F1+I2*F2+I3*F3+I4*F4;
          Complex W(cos(fac*w),sin(fac*w));          
	  index  = (I4%fmax) + (I3%fmax)*oneDim + (I2%fmax)*oneDim*oneDim + (I1%fmax)*oneDim*oneDim*oneDim;
	  index += (F4*fmax) + (F3*fmax)*oneDim + (F2*fmax)*oneDim*oneDim + (F1*fmax)*oneDim*oneDim*oneDim;
          index += O4*offO1+O3*offO2+O2*offO3+O1*offO4;
	  index *= 8;
  	  index += i;
	  res = res + W * input[index];
        }}}}
	index  = (I4) + (I3)*oneDim + (I2)*oneDim*oneDim + (I1)*oneDim*oneDim*oneDim;
        index += O4*offO1+O3*offO2+O2*offO3+O1*offO4;
	index *= 8;
	index += i;
	output[index] = res;
      }
      
    }}}}
  }}}}*/
}



void ExtremeFFT::destroyFourierTrafoAuxData() {
  Complex* c2 = (Complex*) FourierTrafo_bitInverse;
  destroySuperAlignedComplex(c2);
  FourierTrafo_bitInverse = NULL;
  int I;
  int I2 = 0;
  for (I=8; I<=FourierTrafo_iniForN; I*=2) {
    destroySuperAlignedComplex(FourierTrafo_Weights[I2+0]);
    destroySuperAlignedComplex(FourierTrafo_Weights[I2+1]);
    I2+=2;
  }
  delete[] FourierTrafo_Weights;
  destroySuperAlignedComplex(CACHE_L1);
  destroySuperAlignedComplex(CACHE_L2);
  destroySuperAlignedComplex(BUFFER);
}


//N muss Zweier-Potenz sein!!!
void ExtremeFFT::generateFastFourierTrafoAuxData(int L0, int L1, int L2, int L3) {
  if (LogLevel>3) printf("Generating Fast Fourier Trafo Auxilary Data for L0 = %d, L1 = %d, L2 = %d, L3 = %d...\n", L0, L1, L2, L3);

  int N = L3;

  destroyFourierTrafoAuxData();

  CACHE_L1 = createSuperAlignedComplex(4096);
  CACHE_L2 = createSuperAlignedComplex(8*8*(8+xtraCACHESize3)*(8+xtraCACHESize2)*(8+xtraCACHESize1));
  BUFFER = createSuperAlignedComplex(8*N*(N+xtraSize1)*(N+xtraSize2)*(N+xtraSize3));
  int I;
  for (I=0; I<4096; I++) {
    CACHE_L1[I].x = 0;
    CACHE_L1[I].y = 0;
  }
  for (I=0; I<8*8*(8+xtraCACHESize3)*(8+xtraCACHESize2)*(8+xtraCACHESize1); I++) {
    CACHE_L2[I].x = 0;
    CACHE_L2[I].y = 0;
  }
  for (I=0; I<8*N*(N+xtraSize1)*(N+xtraSize2)*(N+xtraSize3); I++) {
    BUFFER[I].x = 0;
    BUFFER[I].y = 0;
  }
  

  int I2 = 0;
  int nr = 0;
  for (I=8; I<=N; I*=2) nr += 2;
  if (nr == 0) {
    FourierTrafo_Weights = NULL;
  } else {
    FourierTrafo_Weights = new Complex*[nr];
  }
  for (I=8; I<=N; I*=2) {
    int count = 4*((I/2)-1);
    FourierTrafo_Weights[I2+0] = createSuperAlignedComplex(2*(count+1));
    FourierTrafo_Weights[I2+1] = createSuperAlignedComplex(2*(count+1));
    int I3;
    for (I3=0; I3<=count; I3++) {
      double w = (2.0*pi*I3) / I;
      FourierTrafo_Weights[I2+0][2*I3+0].x = cos(w);
      FourierTrafo_Weights[I2+0][2*I3+0].y = cos(w);
      FourierTrafo_Weights[I2+0][2*I3+1].x = -sin(w);
      FourierTrafo_Weights[I2+0][2*I3+1].y = sin(w);
      FourierTrafo_Weights[I2+1][2*I3+0].x = cos(w);
      FourierTrafo_Weights[I2+1][2*I3+0].y = cos(w);
      FourierTrafo_Weights[I2+1][2*I3+1].x = sin(w);
      FourierTrafo_Weights[I2+1][2*I3+1].y = -sin(w);
//      FourierTrafo_Weights[0][2*I3+0].print();
//      FourierTrafo_Weights[0][2*I3+1].print();
    }
    I2+=2;
  }
  
  FourierTrafo_iniForN = N;
  int j;
  long int dummy;
  FourierTrafo_bitInverse = (int**)createSuperAlignedComplex(256);
  for (I=0; I<512; I++) FourierTrafo_bitInverse[I] = (int*) 0;
  for (I=0; I<N; I++) {
    j = 1;
    dummy = 0;
    while (j<N) {
      dummy *= 2;
      if ((I & j) > 0) {
        dummy += 1;
      }
      j *= 2;
    }
    FourierTrafo_bitInverse[I+0] = (int*)(16*8*dummy);
    FourierTrafo_bitInverse[I+128] = (int*)(16*8*(N + xtraSize3)*dummy);
    FourierTrafo_bitInverse[I+256] = (int*)(16*8*(N + xtraSize3)*(N + xtraSize2)*dummy);
    FourierTrafo_bitInverse[I+384] = (int*)(16*8*(N + xtraSize3)*(N + xtraSize2)*(N + xtraSize1)*dummy);
    
//    printBits((long int)FourierTrafo_bitInverse[I+0]);
//    printBits((long int)FourierTrafo_bitInverse[I+128]);
//    printBits((long int)FourierTrafo_bitInverse[I+256]);
//    printBits((long int)FourierTrafo_bitInverse[I+384]);
  }
  if (LogLevel>3) printf("ready!!!\n");
}


void ExtremeFFT::FastFourierTrafo(int oneDimL0, int oneDimL1, int oneDimL2, int oneDimL3, Complex* input, Complex* output, bool fft_forward) {
int oneDim = oneDimL3;
  bool performed = false;
  bool targetEqualSource = false;
  long int adr1 = (long int) input;
  long int adr2 = (long int) output;
  long int diffadr = adr1 - adr2;
  if (diffadr<0) diffadr = -diffadr;
  if (diffadr<128*oneDim*(oneDim+xtraSize1)*(oneDim+xtraSize2)*(oneDim+xtraSize3)) targetEqualSource = true;

    
  if (oneDim<=2) {
    printf("xFFT not implemented for L=%d\n",oneDim);
    exit(0);
  }
  if ((FourierTrafo_bitInverse==NULL) || (FourierTrafo_bitInverse==NULL) || (FourierTrafo_iniForN!=oneDim)) {
    generateFastFourierTrafoAuxData(oneDimL0, oneDimL1, oneDimL2, oneDimL3);
  }
  
  Complex* para_Comp = ASMParameterData;
  long int* para_int = (long int*) para_Comp;
  para_int[4] = (long int) FourierTrafo_bitInverse;
  para_int[20] = 0x8000000000000000;
  para_int[21] = 0x8000000000000000;
  if (fft_forward) {
    para_int[22] = 0x8000000000000000;
    para_int[23] = 0;
    para_int[24] = 0;
    para_int[25] = 0x8000000000000000;
  } else {
    para_int[22] = 0;
    para_int[23] = 0x8000000000000000;
    para_int[24] = 0x8000000000000000;
    para_int[25] = 0;
  }
  para_Comp[13] = Complex(2.0, 2.0);
  int I;
  for (I=0; I<16; I++) {
    para_int[50+I] = 0;
  }
  
  int count = 322;
  int L1Size2 = 4;
  int L1Size1 = 16;
  int L1Size0 = 64;  
  int I0, I1, I2, I3;
  for (I0=1; I0>=0; I0--) {
    for (I1=1; I1>=0; I1--) {
      for (I2=1; I2>=0; I2--) {
        for (I3=1; I3>=0; I3--) {
	  para_int[count] = 256*(I3 + I2*L1Size2 + I1*L1Size1 + I0*L1Size0);
	  count++;
	}
      }
    }
  }
  
/*  if (oneDim == 4) {
    para_int[1] = (long int) 128;  
    para_int[2] = (long int) input;
    para_int[3] = (long int) CACHE_L1;

    para_int[28]=0; 
    para_int[32] = para_int[28] + 16;
    para_int[29]=0;
    para_int[33] = para_int[29] + 16;
    para_int[30]=0;
    para_int[34] = para_int[30] + 16;
    para_int[31]=0;
    para_int[35] = para_int[31] + 16;
  
    SSE_FourierTrafoRearrangerWithFirstFFTstep();

    para_int[2] = (long int) CACHE_L1;
    para_int[3] = (long int) output;
    int outputSize = 4;

    para_int[5] = (long int) (128*(outputSize+xtraSize3)*(outputSize+xtraSize2)*(outputSize+xtraSize1) - 128*(outputSize+xtraSize3)*(outputSize+xtraSize2) - 128*(outputSize+xtraSize3) - 128);
    para_int[6] = (long int) (128*(outputSize+xtraSize3)*(outputSize+xtraSize2) - 128*(outputSize+xtraSize3) - 128);
    para_int[7] = (long int) (128*(outputSize+xtraSize3) - 128);
    para_int[8] = (long int) (128);
    para_int[9] = (long int) (2*128*(outputSize+xtraSize3));
    para_int[10] = (long int) (2*128*(outputSize+xtraSize3)*(outputSize+xtraSize2));
    para_int[11] = (long int) (para_int[9] + para_int[10]);
    para_int[12] = (long int) (2*128*(outputSize+xtraSize3)*(outputSize+xtraSize2)*(outputSize+xtraSize1));
  
    SSE_FourierTrafoSecondFFTstep();
    performed = true;
  }*/
  
/*  Complex* CACHE_L2_ADR = output;
  if (targetEqualSource) {
    CACHE_L2_ADR = CACHE_L2;
  }
  if (oneDim == 8) {
    int maxIndex = 8*oneDim;
    para_int[28]=0; 
    while (true) {
      para_int[32] = para_int[28] + 16;
      para_int[29]=0;
      while (true) {
        para_int[33] = para_int[29] + 16;
        para_int[30]=0;
        while (true) {
          para_int[34] = para_int[30] + 16;
	  para_int[31]=0;
          while (true) {
            para_int[35] = para_int[31] + 16;

            para_int[1] = (long int) 128;  
            para_int[2] = (long int) input;
            para_int[3] = (long int) CACHE_L1;
  
            SSE_FourierTrafoRearrangerWithFirstFFTstep();

            int outputSize = 8;
            para_int[2] = (long int) CACHE_L1;
            para_int[3] = (long int) CACHE_L2_ADR;
	    para_int[3] += para_int[31]*16 + para_int[30]*16*(outputSize+xtraCACHESize3) + para_int[29]*16*(outputSize+xtraCACHESize3)*(outputSize+xtraCACHESize2)+ para_int[28]*16*(outputSize+xtraCACHESize3)*(outputSize+xtraCACHESize2)*(outputSize+xtraCACHESize1);
   
            para_int[5] = (long int) (128*(outputSize+xtraCACHESize3)*(outputSize+xtraCACHESize2)*(outputSize+xtraCACHESize1) - 128*(outputSize+xtraCACHESize3)*(outputSize+xtraCACHESize2) - 128*(outputSize+xtraCACHESize3) - 128);
            para_int[6] = (long int) (128*(outputSize+xtraCACHESize3)*(outputSize+xtraCACHESize2) - 128*(outputSize+xtraCACHESize3) - 128);
            para_int[7] = (long int) (128*(outputSize+xtraCACHESize3) - 128);
            para_int[8] = (long int) (128);
            para_int[9] = (long int) (2*128*(outputSize+xtraCACHESize3));
            para_int[10] = (long int) (2*128*(outputSize+xtraCACHESize3)*(outputSize+xtraCACHESize2));
            para_int[11] = (long int) (para_int[9] + para_int[10]);
            para_int[12] = (long int) (2*128*(outputSize+xtraCACHESize3)*(outputSize+xtraCACHESize2)*(outputSize+xtraCACHESize1));
  
            SSE_FourierTrafoSecondFFTstep();

  	    para_int[31]+=32;
            if (para_int[31]>=maxIndex) break;
  	  }
          para_int[30]+=32;
          if (para_int[30]>=maxIndex) break;
        }
        para_int[29]+=32;
        if (para_int[29]>=maxIndex) break;
      }
      para_int[28]+=32;
      if (para_int[28]>=maxIndex) break;
    }

    para_int[1] = (long int) 128;  
    para_int[2] = (long int) CACHE_L2_ADR;
    para_int[3] = (long int) output;

    if (fft_forward) {
      para_int[13] = (long int) (FourierTrafo_Weights[1]);
    } else {
      para_int[13] = (long int) (FourierTrafo_Weights[0]);
    }
    int L2Size = 8;
    para_int[36] = (long int) (128*4*(L2Size+xtraCACHESize3));
    para_int[37] = (long int) (128*4*(L2Size+xtraCACHESize3)*(L2Size+xtraCACHESize2));
    para_int[38] = (long int) para_int[36] + para_int[37];
    para_int[39] = (long int) (128*4*(L2Size+xtraCACHESize3)*(L2Size+xtraCACHESize2)*(L2Size+xtraCACHESize1));

    para_int[14] = (long int) (128*(L2Size+xtraCACHESize3)*(L2Size+xtraCACHESize2)*(L2Size+xtraCACHESize1) - 4*128*(L2Size+xtraCACHESize3)*(L2Size+xtraCACHESize2));
    para_int[15] = (long int) (128*(L2Size+xtraCACHESize3)*(L2Size+xtraCACHESize2) - 4*128*(L2Size+xtraCACHESize3));
    para_int[16] = (long int) (128*(L2Size+xtraCACHESize3) - 512);

    int outSize = 8;
    para_int[40] = (long int) (128*4*(outSize+xtraSize3));
    para_int[41] = (long int) (128*4*(outSize+xtraSize3)*(outSize+xtraSize2));
    para_int[42] = (long int) para_int[40] + para_int[41];
    para_int[43] = (long int) (128*4*(outSize+xtraSize3)*(outSize+xtraSize2)*(outSize+xtraSize1));
    
    para_int[17] = (long int) (128*(outSize+xtraSize3)*(outSize+xtraSize2)*(outSize+xtraSize1) - 4*128*(outSize+xtraSize3)*(outSize+xtraSize2));
    para_int[18] = (long int) (128*(outSize+xtraSize3)*(outSize+xtraSize2) - 4*128*(outSize+xtraSize3));
    para_int[19] = (long int) (128*(outSize+xtraSize3) - 512);
    
    
    SSE_FourierTrafoThirdFFTstep();
    performed = true;
  }*/
  
  
  if (oneDim == 16) {
    int L2Size = 8;
    int maxAindex = oneDim*8;

    int a0, a1, a2 , a3;
    int L2p0, L2p1 = 0, L2p2, L2p3 = 0;

    a0 = 0;
    while (true) {
      a1 = 0;
      while (true) {
        a2 = 0;
        while (true) {
	  a3 = 0;
          while (true) {
	  
  for (I=0; I<16; I++) {
    para_int[50+I] = 0;
  }
  count = 322;
  for (I0=1; I0>=0; I0--) {
    for (I1=1; I1>=0; I1--) {
      for (I2=1; I2>=0; I2--) {
        for (I3=1; I3>=0; I3--) {
	  para_int[count] = 256*(I3 + I2*L1Size2 + I1*L1Size1 + I0*L1Size0);
	  count++;
	}
      }
    }
  }
	  
	  
            para_int[28]=a0; 
            L2p0 = 0;
            while (true) {
              para_int[32] = para_int[28] + 16;
              para_int[29]=a1;
              L2p1 = 0;
              while (true) {
                para_int[33] = para_int[29] + 16;
                para_int[30]=a2;
                L2p2 = 0;
                while (true) {
                  para_int[34] = para_int[30] + 16;
        	  para_int[31]=a3;
                  L2p3 = 0;
                  while (true) {
                    para_int[35] = para_int[31] + 16;

                    para_int[1] = (long int) 128;  
                    para_int[2] = (long int) input;
                    para_int[3] = (long int) CACHE_L1;
  
                    SSE_FourierTrafoRearrangerWithFirstFFTstep();

                    para_int[2] = (long int) CACHE_L1;
                    para_int[3] = (long int) CACHE_L2;
   	            para_int[3] += L2p3*16 + L2p2*16*(L2Size+xtraCACHESize3) + L2p1*16*(L2Size+xtraCACHESize3)*(L2Size+xtraCACHESize2)+ L2p0*16*(L2Size+xtraCACHESize3)*(L2Size+xtraCACHESize2)*(L2Size+xtraCACHESize1);
   
                    para_int[5] = (long int) (128*(L2Size+xtraCACHESize3)*(L2Size+xtraCACHESize2)*(L2Size+xtraCACHESize1) - 128*(L2Size+xtraCACHESize3)*(L2Size+xtraCACHESize2) - 128*(L2Size+xtraCACHESize3) - 128);
                    para_int[6] = (long int) (128*(L2Size+xtraCACHESize3)*(L2Size+xtraCACHESize2) - 128*(L2Size+xtraCACHESize3) - 128);
                    para_int[7] = (long int) (128*(L2Size+xtraCACHESize3) - 128);
                    para_int[8] = (long int) (128);
                    para_int[9] = (long int) (2*128*(L2Size+xtraCACHESize3));
                    para_int[10] = (long int) (2*128*(L2Size+xtraCACHESize3)*(L2Size+xtraCACHESize2));
                    para_int[11] = (long int) (para_int[9] + para_int[10]);
                    para_int[12] = (long int) (2*128*(L2Size+xtraCACHESize3)*(L2Size+xtraCACHESize2)*(L2Size+xtraCACHESize1));
  
  
                    SSE_FourierTrafoSecondFFTstep();

                    L2p3 += 32;
        	    para_int[31]+=32;
                    if (para_int[31]>=a3+64) break;
                  }
                  L2p2 += 32;
                  para_int[30]+=32;
                  if (para_int[30]>=a2+64) break;
                }
                L2p1 += 32;
                para_int[29]+=32;
                if (para_int[29]>=a1+64) break;
              }
              L2p0 += 32;
              para_int[28]+=32;
              if (para_int[28]>=a0+64) break;
            }

            para_int[1] = (long int) 128;  
            para_int[2] = (long int) CACHE_L2;
            para_int[3] = (long int) BUFFER;
	    
            if (fft_forward) {
              para_int[13] = (long int) (FourierTrafo_Weights[1]);
            } else {
              para_int[13] = (long int) (FourierTrafo_Weights[0]);
            }

            para_int[36] = (long int) (128*4*(L2Size+xtraCACHESize3));
            para_int[37] = (long int) (128*4*(L2Size+xtraCACHESize3)*(L2Size+xtraCACHESize2));
            para_int[38] = (long int) para_int[36] + para_int[37];
            para_int[39] = (long int) (128*4*(L2Size+xtraCACHESize3)*(L2Size+xtraCACHESize2)*(L2Size+xtraCACHESize1));

            para_int[14] = (long int) (128*(L2Size+xtraCACHESize3)*(L2Size+xtraCACHESize2)*(L2Size+xtraCACHESize1) - 4*128*(L2Size+xtraCACHESize3)*(L2Size+xtraCACHESize2));
            para_int[15] = (long int) (128*(L2Size+xtraCACHESize3)*(L2Size+xtraCACHESize2) - 4*128*(L2Size+xtraCACHESize3));
            para_int[16] = (long int) (128*(L2Size+xtraCACHESize3));

            int outSize = 16;
	    para_int[3] += a3*16 + a2*16*(outSize+xtraSize3) + a1*16*(outSize+xtraSize3)*(outSize+xtraSize2)+ a0*16*(outSize+xtraSize3)*(outSize+xtraSize2)*(outSize+xtraSize1);
	    
            para_int[40] = (long int) (128*4*(outSize+xtraSize3));
            para_int[41] = (long int) (128*4*(outSize+xtraSize3)*(outSize+xtraSize2));
            para_int[42] = (long int) para_int[40] + para_int[41];
            para_int[43] = (long int) (128*4*(outSize+xtraSize3)*(outSize+xtraSize2)*(outSize+xtraSize1));
    
            para_int[17] = (long int) (128*(outSize+xtraSize3)*(outSize+xtraSize2)*(outSize+xtraSize1) - 4*128*(outSize+xtraSize3)*(outSize+xtraSize2));
            para_int[18] = (long int) (128*(outSize+xtraSize3)*(outSize+xtraSize2) - 4*128*(outSize+xtraSize3));
            para_int[19] = (long int) (128*(outSize+xtraSize3));


    count = 320;
    int pos = 0;
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = 512;
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[36];
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[36]+512;
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[37];
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[37]+512;
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[38];
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[38]+512;
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }

    pos = para_int[39];
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[39]+512;
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[39]+para_int[36];
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[39]+para_int[36]+512;
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[39]+para_int[37];
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[39]+para_int[37]+512;
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[39]+para_int[38];
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[39]+para_int[38]+512;
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }            
    
    pos = 0;
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = 512;
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[40];
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[40]+512;
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[41];
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[41]+512;
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[42];
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[42]+512;
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }

    pos = para_int[43];
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[43]+512;
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[43]+para_int[40];
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[43]+para_int[40]+512;
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[43]+para_int[41];
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[43]+para_int[41]+512;
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[43]+para_int[42];
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[43]+para_int[42]+512;
    for (I=0; I<8; I++) { para_int[count] = pos+I*64; count++; }            
    
    
    
            para_int[50] = 64;
            para_int[0] = (long int) 128;  
	    
	    for (I=0; I<2; I++) {
              SSE_FourierTrafoThirdFFTstep();
              para_int[2] += 64;
              para_int[3] += 64;
              para_int[50] = 0;
	    }
	    
	    a3 += 64;
	    if (a3>=maxAindex) break;
          }
          a2 += 64;
          if (a2>=maxAindex) break;
        }
	a1 += 64;
        if (a1>=maxAindex) break;
      }
      a0 += 64;
      if (a0>=maxAindex) break;
    }

    //Weitere Transformationen...
    para_int[0] = (long int) 256;  
    para_int[1] = (long int) 128;  
    para_int[2] = (long int) BUFFER;
    para_int[3] = (long int) output;

    if (fft_forward) {
      para_int[13] = (long int) (FourierTrafo_Weights[3]);
    } else {
      para_int[13] = (long int) (FourierTrafo_Weights[2]);
    }

    int outSize = 16;
    para_int[36] = (long int) (128*8*(outSize+xtraSize3));
    para_int[37] = (long int) (128*8*(outSize+xtraSize3)*(outSize+xtraSize2));
    para_int[38] = (long int) para_int[36] + para_int[37];
    para_int[39] = (long int) (128*8*(outSize+xtraSize3)*(outSize+xtraSize2)*(outSize+xtraSize1));

    para_int[14] = (long int) (128*(outSize+xtraSize3)*(outSize+xtraSize2)*(outSize+xtraSize1) - 8*128*(outSize+xtraSize3)*(outSize+xtraSize2));
    para_int[15] = (long int) (128*(outSize+xtraSize3)*(outSize+xtraSize2) - 8*128*(outSize+xtraSize3));
    para_int[16] = (long int) (128*(outSize+xtraSize3));
    

    count = 320;
    int pos = 0;
    for (I=0; I<16; I++) { para_int[count] = pos+I*64; count++; }
    pos = 1024;
    for (I=0; I<16; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[36];
    for (I=0; I<16; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[36]+1024;
    for (I=0; I<16; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[37];
    for (I=0; I<16; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[37]+1024;
    for (I=0; I<16; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[38];
    for (I=0; I<16; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[38]+1024;
    for (I=0; I<16; I++) { para_int[count] = pos+I*64; count++; }

    pos = para_int[39];
    for (I=0; I<16; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[39]+1024;
    for (I=0; I<16; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[39]+para_int[36];
    for (I=0; I<16; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[39]+para_int[36]+1024;
    for (I=0; I<16; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[39]+para_int[37];
    for (I=0; I<16; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[39]+para_int[37]+1024;
    for (I=0; I<16; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[39]+para_int[38];
    for (I=0; I<16; I++) { para_int[count] = pos+I*64; count++; }
    pos = para_int[39]+para_int[38]+1024;
    for (I=0; I<16; I++) { para_int[count] = pos+I*64; count++; }


    para_int[17] = 64;   //Interlacing - Einstellung [18] intern verwendet
    for (I=0; I<2; I++) {
      SSE_FourierTrafoForthFFTstep();
      para_int[2] += 64;
      para_int[3] += 64;
      para_int[17] = 0;
    }    
    
    performed = true;
  }  
  
  if (!performed) {
    printf("xFFT not implemented for L=%d\n",oneDim);
    exit(0);
  }
}
