#ifndef EvaluateObservableHiggsCorrelator_included
#define EvaluateObservableHiggsCorrelator_included

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


class EvaluateObservableHiggsCorrelator : public EvaluateObservableCorrelatorBase {
private:  
  
public:    
  EvaluateObservableHiggsCorrelator(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableHiggsCorrelator();

};


#include "EvaluateObservableHiggsCorrelator.C"

#endif
