#ifndef ExtremeFFT_included
#define ExtremeFFT_included

#include <math.h>
#include "Tools.h"

#define ExtremeFFT_Forward true
#define ExtremeFFT_Backward false


class ExtremeFFT {
private:
int** FourierTrafo_bitInverse;
int FourierTrafo_iniForL0;
int FourierTrafo_iniForL1;
int FourierTrafo_iniForL2;
int FourierTrafo_iniForL3;
int FourierTrafo_iniForN;
Complex** FourierTrafo_Weights;
Complex* CACHE_L1;
Complex* CACHE_L2;
Complex* BUFFER;



public:
  ExtremeFFT(); 
  
  
  void SlowRearrange(int oneDimL0, int oneDimL1, int oneDimL2, int oneDimL3, Complex* input, Complex* output);
  void SlowFourierStep(int oneDimL0, int oneDimL1, int oneDimL2, int oneDimL3, int stepNr, Complex* input, Complex* output);
  
  void destroyFourierTrafoAuxData();
  void generateFastFourierTrafoAuxData(int L0, int L1, int L2, int L3);
  void FastFourierTrafo(int oneDimL0, int oneDimL1, int oneDimL2, int oneDimL3, Complex* input, Complex* output, bool fft_forward);
  
};


#endif
