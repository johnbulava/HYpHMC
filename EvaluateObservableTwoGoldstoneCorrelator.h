#ifndef EvaluateObservableTwoGoldstoneCorrelator_included
#define EvaluateObservableTwoGoldstoneCorrelator_included

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


class EvaluateObservableTwoGoldstoneCorrelator : public EvaluateObservableCorrelatorBase {
private:  
  
public:    
	EvaluateObservableTwoGoldstoneCorrelator(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
	~EvaluateObservableTwoGoldstoneCorrelator();

};


#include "EvaluateObservableTwoGoldstoneCorrelator.C"

#endif
