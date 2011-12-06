#ifndef EvaluateObservableFermionMatrixSingleMConditionNumberNoPreconditioning_included
#define EvaluateObservableFermionMatrixSingleMConditionNumberNoPreconditioning_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "EvaluateObservableFermionMatrixConditionNumberBase.h"


class EvaluateObservableFermionMatrixSingleMConditionNumberNoPreconditioning : public EvaluateObservableFermionMatrixConditionNumberBase {
private:  

public:    
  EvaluateObservableFermionMatrixSingleMConditionNumberNoPreconditioning(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableFermionMatrixSingleMConditionNumberNoPreconditioning();

};


#include "EvaluateObservableFermionMatrixSingleMConditionNumberNoPreconditioning.C"

#endif
