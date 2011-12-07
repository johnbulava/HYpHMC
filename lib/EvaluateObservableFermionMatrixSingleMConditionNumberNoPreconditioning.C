#include "EvaluateObservableFermionMatrixSingleMConditionNumberNoPreconditioning.h"

EvaluateObservableFermionMatrixSingleMConditionNumberNoPreconditioning::EvaluateObservableFermionMatrixSingleMConditionNumberNoPreconditioning(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservableFermionMatrixConditionNumberBase(aIOcon, sdr, "FermionMatrixSingleMConditionNumberNoPreconditioning", "condsmnp", obsWeight, obsDetSign, relStart, relEnd) {
  drawUpperBound = false;
  drawLowerBound = false;
}


EvaluateObservableFermionMatrixSingleMConditionNumberNoPreconditioning::~EvaluateObservableFermionMatrixSingleMConditionNumberNoPreconditioning() {
}

