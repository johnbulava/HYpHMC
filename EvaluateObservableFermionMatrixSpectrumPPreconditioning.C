EvaluateObservableFermionMatrixSpectrumPPreconditioning::EvaluateObservableFermionMatrixSpectrumPPreconditioning(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservableFermionMatrixSpectrumBase(aIOcon, sdr, "FermionMatrixSingleMSpectrumPPreconditioning", "specsmpp", obsWeight, obsDetSign, relStart, relEnd) {
  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
}


EvaluateObservableFermionMatrixSpectrumPPreconditioning::~EvaluateObservableFermionMatrixSpectrumPPreconditioning() {
}


int EvaluateObservableFermionMatrixSpectrumPPreconditioning::getAnalyzerResultsCount() {
  return 40;
}
