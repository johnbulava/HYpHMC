#ifndef EvaluateObservableHiggsPropagator_included
#define EvaluateObservableHiggsPropagator_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservablePropagatorBase.h"
#include "SimulationParameterSet.h"


class EvaluateObservableHiggsPropagator : public EvaluateObservablePropagatorBase {
private:  
	
public:    
  EvaluateObservableHiggsPropagator(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableHiggsPropagator();
};


#include "EvaluateObservableHiggsPropagator.C"

#endif
