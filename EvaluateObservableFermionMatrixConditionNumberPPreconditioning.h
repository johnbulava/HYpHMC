#ifndef EvaluateObservableFermionMatrixConditionNumberPPreconditioning_included
#define EvaluateObservableFermionMatrixConditionNumberPPreconditioning_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "EvaluateObservableFermionMatrixConditionNumberBase.h"


class EvaluateObservableFermionMatrixConditionNumberPPreconditioning : public EvaluateObservableFermionMatrixConditionNumberBase {
private:  

public:    
  EvaluateObservableFermionMatrixConditionNumberPPreconditioning(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableFermionMatrixConditionNumberPPreconditioning();

};


#include "EvaluateObservableFermionMatrixConditionNumberPPreconditioning.C"

#endif
