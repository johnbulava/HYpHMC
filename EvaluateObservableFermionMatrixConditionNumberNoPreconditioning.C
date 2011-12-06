EvaluateObservableFermionMatrixConditionNumberNoPreconditioning::EvaluateObservableFermionMatrixConditionNumberNoPreconditioning(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservableFermionMatrixConditionNumberBase(aIOcon, sdr, "FermionMatrixConditionNumberNoPreconditioning", "condnp", obsWeight, obsDetSign, relStart, relEnd) {
  drawUpperBound = false;
  drawLowerBound = false;
}


EvaluateObservableFermionMatrixConditionNumberNoPreconditioning::~EvaluateObservableFermionMatrixConditionNumberNoPreconditioning() {
}

