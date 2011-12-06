#ifndef EvaluateObservableHiggsDefinedOnTimeSliceCorrelator_included
#define EvaluateObservableHiggsDefinedOnTimeSliceCorrelator_included

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


class EvaluateObservableHiggsDefinedOnTimeSliceCorrelator : public EvaluateObservableCorrelatorBase {
private:  
  
public:    
  EvaluateObservableHiggsDefinedOnTimeSliceCorrelator(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableHiggsDefinedOnTimeSliceCorrelator();

};


#include "EvaluateObservableHiggsDefinedOnTimeSliceCorrelator.C"

#endif
