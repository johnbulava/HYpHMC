EvaluateObservableFermionMatrixConditionNumberPPreconditioning::EvaluateObservableFermionMatrixConditionNumberPPreconditioning(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservableFermionMatrixConditionNumberBase(aIOcon, sdr, "FermionMatrixConditionNumberPPreconditioning", "condpp", obsWeight, obsDetSign, relStart, relEnd) {
  drawUpperBound = false;
  drawLowerBound = false;
}


EvaluateObservableFermionMatrixConditionNumberPPreconditioning::~EvaluateObservableFermionMatrixConditionNumberPPreconditioning() {
}

