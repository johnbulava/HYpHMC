#include "EvaluateObservableFermionMatrixConditionNumber.h"

EvaluateObservableFermionMatrixConditionNumber::EvaluateObservableFermionMatrixConditionNumber(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservableFermionMatrixConditionNumberBase(aIOcon, sdr, "FermionMatrixConditionNumber", "cond", obsWeight, obsDetSign, relStart, relEnd) {
  drawUpperBound = true;
  drawLowerBound = true;
}


EvaluateObservableFermionMatrixConditionNumber::~EvaluateObservableFermionMatrixConditionNumber() {
}

