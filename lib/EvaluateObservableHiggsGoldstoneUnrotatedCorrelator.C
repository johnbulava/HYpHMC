EvaluateObservableHiggsGoldstoneUnrotatedCorrelator::EvaluateObservableHiggsGoldstoneUnrotatedCorrelator(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservableCorrelatorBase(aIOcon, sdr, obsWeight, obsDetSign, "HiggsGoldstoneUnrotatedCorrelator", "hgucorr", relStart, relEnd, 4,4) { 
  doMassCorrMatrixAnalysis = true;
  doSeparateAnalysis = false;
}


EvaluateObservableHiggsGoldstoneUnrotatedCorrelator::~EvaluateObservableHiggsGoldstoneUnrotatedCorrelator() {
}
