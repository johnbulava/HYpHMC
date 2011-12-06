#ifndef EvaluateObservableFermionMatrixConditionNumber_included
#define EvaluateObservableFermionMatrixConditionNumber_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "EvaluateObservableFermionMatrixConditionNumberBase.h"


class EvaluateObservableFermionMatrixConditionNumber : public EvaluateObservableFermionMatrixConditionNumberBase {
private:  

public:    
  EvaluateObservableFermionMatrixConditionNumber(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableFermionMatrixConditionNumber();

};


#include "EvaluateObservableFermionMatrixConditionNumber.C"

#endif
