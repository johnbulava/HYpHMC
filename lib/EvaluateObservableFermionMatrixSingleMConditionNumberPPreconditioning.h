#ifndef EvaluateObservableFermionMatrixSingleMConditionNumberPPreconditioning_included
#define EvaluateObservableFermionMatrixSingleMConditionNumberPPreconditioning_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "EvaluateObservableFermionMatrixConditionNumberBase.h"


class EvaluateObservableFermionMatrixSingleMConditionNumberPPreconditioning : public EvaluateObservableFermionMatrixConditionNumberBase {
private:  

public:    
  EvaluateObservableFermionMatrixSingleMConditionNumberPPreconditioning(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableFermionMatrixSingleMConditionNumberPPreconditioning();

};


#endif
