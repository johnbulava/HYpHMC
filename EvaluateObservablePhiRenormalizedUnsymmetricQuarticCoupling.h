#ifndef EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling_included
#define EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "LatticeMomentumBins.h"
#include "AutoCorrelation.h"


class EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling : public EvaluateObservable {
private:  
  LatticeMomentumBins* latticeBins;
  AutoCorrelation* autoCorr;
  double* pSqr;
  double* avgLamRen;
  double* sigmaLamRen;
  double* sigmaLamRenHelper;  
  double* autoCorrelationTime;
  double* propTemp;  
  int smallL0;
  int smallL1;
  int smallL2;
  int smallL3;
  
  void calcLamRen(int ignoreStart, int ignoreEnd);
  int getIndex(int i0, int i1, int i2, int i3, int locInd);
  LAPsystemPlot* createPlot1(double maxP);
  LAPsystemPlot* createPlot2();
  
protected:

  
public:    
  EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling();

  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
  int getAnalyzerResultsCount();      
};


#include "EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling.C"

#endif
