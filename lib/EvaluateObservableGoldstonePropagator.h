#ifndef EvaluateObservableGoldstonePropagator_included
#define EvaluateObservableGoldstonePropagator_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservablePropagatorBase.h"


class EvaluateObservableGoldstonePropagator : public EvaluateObservablePropagatorBase {
private:  
	
public:    
  EvaluateObservableGoldstonePropagator(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableGoldstonePropagator();
};


#endif
