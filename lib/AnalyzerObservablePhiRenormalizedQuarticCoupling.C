AnalyzerObservablePhiRenormalizedQuarticCoupling::AnalyzerObservablePhiRenormalizedQuarticCoupling(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "PhiRenormalizedQuarticCoupling", "phirlam") { 
  ini(getAnalyzerResultsCount());
}


AnalyzerObservablePhiRenormalizedQuarticCoupling::~AnalyzerObservablePhiRenormalizedQuarticCoupling() {
}


bool AnalyzerObservablePhiRenormalizedQuarticCoupling::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  analyzerResults[0] = 0;

  return true;
}


int AnalyzerObservablePhiRenormalizedQuarticCoupling::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservablePhiRenormalizedQuarticCoupling::getAnalyzerResultsCount() {
  return 1;
}
