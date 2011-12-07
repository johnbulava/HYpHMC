#ifndef EvaluateObservableHiggsRenormalizedQuarticCoupling_included
#define EvaluateObservableHiggsRenormalizedQuarticCoupling_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "LatticeMomentumBins.h"
#include "AutoCorrelation.h"


class EvaluateObservableHiggsRenormalizedQuarticCoupling : public EvaluateObservable {
private:  
  LatticeMomentumBins* latticeBins;
  AutoCorrelation* autoCorr;
  double* pSqr;
  double* avgLamRen;
  double* sigmaLamRen;
  double* sigmaLamRenHelper;  
  double* autoCorrelationTimeHiggs4;
  bool considerZeroMomentum;
  
  void calcLamRen(int ignoreStart, int ignoreEnd);
  LAPsystemPlot* createPlot1(double maxP);
  LAPsystemPlot* createPlot2();
  
protected:

  
public:    
  EvaluateObservableHiggsRenormalizedQuarticCoupling(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableHiggsRenormalizedQuarticCoupling();

  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
  int getAnalyzerResultsCount();      
};


#endif
