EvaluateObservableTwoGoldstoneCorrelator::EvaluateObservableTwoGoldstoneCorrelator(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservableCorrelatorBase(aIOcon, sdr, obsWeight, obsDetSign, "TwoGoldstoneCorrelator", "twogoldstonecorr", relStart, relEnd, 1,1) { 
  doMassCorrMatrixAnalysis = true;
  doSeparateAnalysis = false;
}


EvaluateObservableTwoGoldstoneCorrelator::~EvaluateObservableTwoGoldstoneCorrelator() {
}
