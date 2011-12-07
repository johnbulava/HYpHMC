#ifndef EvaluateObservableHiggsGoldstoneUnrotatedCorrelator_included
#define EvaluateObservableHiggsGoldstoneUnrotatedCorrelator_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "AutoCorrelation.h"
#include "MassCorrelationMatrixAnalyzer.h"
#include "EvaluateObservableCorrelatorBase.h"


class EvaluateObservableHiggsGoldstoneUnrotatedCorrelator : public EvaluateObservableCorrelatorBase {
private:  
  
public:    
  EvaluateObservableHiggsGoldstoneUnrotatedCorrelator(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableHiggsGoldstoneUnrotatedCorrelator();

};


#endif
