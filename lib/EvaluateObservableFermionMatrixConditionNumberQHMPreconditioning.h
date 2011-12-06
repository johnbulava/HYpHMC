#ifndef EvaluateObservableFermionMatrixConditionNumberQHMPreconditioning_included
#define EvaluateObservableFermionMatrixConditionNumberQHMPreconditioning_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "EvaluateObservableFermionMatrixConditionNumberBase.h"


class EvaluateObservableFermionMatrixConditionNumberQHMPreconditioning : public EvaluateObservableFermionMatrixConditionNumberBase {
private:  

public:    
  EvaluateObservableFermionMatrixConditionNumberQHMPreconditioning(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableFermionMatrixConditionNumberQHMPreconditioning();

};


#include "EvaluateObservableFermionMatrixConditionNumberQHMPreconditioning.C"

#endif
