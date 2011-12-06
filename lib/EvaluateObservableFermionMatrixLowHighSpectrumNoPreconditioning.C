EvaluateObservableFermionMatrixLowHighSpectrumNoPreconditioning::EvaluateObservableFermionMatrixLowHighSpectrumNoPreconditioning(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservableFermionMatrixSpectrumBase(aIOcon, sdr, "FermionMatrixSingleMLowHighSpectrumNoPreconditioning", "lhspecsmnp", obsWeight, obsDetSign, relStart, relEnd) {
  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
}


EvaluateObservableFermionMatrixLowHighSpectrumNoPreconditioning::~EvaluateObservableFermionMatrixLowHighSpectrumNoPreconditioning() {
}


int EvaluateObservableFermionMatrixLowHighSpectrumNoPreconditioning::getAnalyzerResultsCount() {
  return 240;
}
