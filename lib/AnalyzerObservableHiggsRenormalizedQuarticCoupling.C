AnalyzerObservableHiggsRenormalizedQuarticCoupling::AnalyzerObservableHiggsRenormalizedQuarticCoupling(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "HiggsRenormalizedQuarticCoupling", "hrlam") { 
  ini(getAnalyzerResultsCount());
}


AnalyzerObservableHiggsRenormalizedQuarticCoupling::~AnalyzerObservableHiggsRenormalizedQuarticCoupling() {
}


bool AnalyzerObservableHiggsRenormalizedQuarticCoupling::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  analyzerResults[0] = 0;

  return true;
}


int AnalyzerObservableHiggsRenormalizedQuarticCoupling::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservableHiggsRenormalizedQuarticCoupling::getAnalyzerResultsCount() {
  return 1;
}
