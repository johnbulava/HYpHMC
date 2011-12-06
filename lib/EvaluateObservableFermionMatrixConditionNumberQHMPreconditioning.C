EvaluateObservableFermionMatrixConditionNumberQHMPreconditioning::EvaluateObservableFermionMatrixConditionNumberQHMPreconditioning(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservableFermionMatrixConditionNumberBase(aIOcon, sdr, "FermionMatrixConditionNumberQHMPreconditioning", "condqhm", obsWeight, obsDetSign, relStart, relEnd) {
  drawUpperBound = false;
  drawLowerBound = false;
}


EvaluateObservableFermionMatrixConditionNumberQHMPreconditioning::~EvaluateObservableFermionMatrixConditionNumberQHMPreconditioning() {
}

