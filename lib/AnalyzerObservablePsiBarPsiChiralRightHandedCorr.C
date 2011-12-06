AnalyzerObservablePsiBarPsiChiralRightHandedCorr::AnalyzerObservablePsiBarPsiChiralRightHandedCorr(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservablePsiBarPsiCorrBase(fOps, aIOcon, SDreader, "PsiBarPsiCorrChiralRightHanded", "pbpcxr") { 
  fixGauge = false;
  randomGauge = false;
  projectorSelection = 1;
  multiplyWithPhiMatBSelection = 0;  
}


AnalyzerObservablePsiBarPsiChiralRightHandedCorr::~AnalyzerObservablePsiBarPsiChiralRightHandedCorr() {
}
