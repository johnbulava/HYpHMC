#ifndef EvaluateObservableFermionMatrixConditionNumberNoPreconditioning_included
#define EvaluateObservableFermionMatrixConditionNumberNoPreconditioning_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "EvaluateObservableFermionMatrixConditionNumberBase.h"


class EvaluateObservableFermionMatrixConditionNumberNoPreconditioning : public EvaluateObservableFermionMatrixConditionNumberBase {
private:  

public:    
  EvaluateObservableFermionMatrixConditionNumberNoPreconditioning(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableFermionMatrixConditionNumberNoPreconditioning();

};


#endif
