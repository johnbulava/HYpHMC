#include "EvaluateObservableFermionMatrixSingleMConditionNumberPPreconditioning.h"

EvaluateObservableFermionMatrixSingleMConditionNumberPPreconditioning::EvaluateObservableFermionMatrixSingleMConditionNumberPPreconditioning(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservableFermionMatrixConditionNumberBase(aIOcon, sdr, "FermionMatrixSingleMConditionNumberPPreconditioning", "condsmpp", obsWeight, obsDetSign, relStart, relEnd) {
  drawUpperBound = false;
  drawLowerBound = false;
}


EvaluateObservableFermionMatrixSingleMConditionNumberPPreconditioning::~EvaluateObservableFermionMatrixSingleMConditionNumberPPreconditioning() {
}

