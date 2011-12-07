#ifndef AutoCorrelation_included
#define AutoCorrelation_included

#include <stdlib.h>

#include "Global.h"
#include "Complex.h"
#include "Quat.h"
#include "ComplexVector.h"
#include "ComplexMatrix.h"



class AutoCorrelation {
protected:
  void ini(int obsCount, int GWmax);
  void desini();
  
  int observableCount;
  int GammaWMax;
  double*** Gammas;
  double*** GammaSigmas;  
  double* CombinedGamma;
  double* CombinedGammaSigma;
  double* WstepArray;
  double* CombinedCFunction;
  double* CombinedReducedCFunction;
  double* CombinedCFunctionErrors;
  double* CombinedReducedCFunctionErrors;
  double* averages;
  int totalN;
  void calcCombinedGammaFunction(ComplexVector& derivatives);
  double AutoCorrelationTime;
  
    
public:
  AutoCorrelation();
  AutoCorrelation(int obsCount, int GWmax);
  ~AutoCorrelation();
  
  void loadData(int RunCount, int* RunLengths, double* detData, double* measureData);
  double getAverage(int nr);
  void calcCombinedCFunction(ComplexVector& derivatives);
  double* getCombinedGamma();
  double* getCombinedCFunction();
  double* getCombinedReducedCFunction();  
  double* getCombinedCFunctionErrors();
  double* getCombinedReducedCFunctionErrors();
  double* getWstepArray();
  int getGammaWMax();
  int getTotalN();
  double estimateCombinedError(ComplexVector& derivatives);
  double estimateAutoCorrelationTime(ComplexVector& derivatives);
  double estimateAutoCorrelationTime();  
};

#endif
