EvaluateObservableHiggsCorrelator::EvaluateObservableHiggsCorrelator(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservableCorrelatorBase(aIOcon, sdr, obsWeight, obsDetSign, "HiggsCorrelator", "hcorr", relStart, relEnd, 1,1) { 
  doMassCorrMatrixAnalysis = false;
  doSeparateAnalysis = true;
}


EvaluateObservableHiggsCorrelator::~EvaluateObservableHiggsCorrelator() {
}
