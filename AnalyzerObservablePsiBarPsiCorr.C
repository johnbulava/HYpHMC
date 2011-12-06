AnalyzerObservablePsiBarPsiCorr::AnalyzerObservablePsiBarPsiCorr(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservablePsiBarPsiCorrBase(fOps, aIOcon, SDreader, "PsiBarPsiCorr", "pbpc") { 
  fixGauge = false;
  randomGauge = false;
  projectorSelection = 0;
  multiplyWithPhiMatBSelection = 0;  
}


AnalyzerObservablePsiBarPsiCorr::~AnalyzerObservablePsiBarPsiCorr() {
}
