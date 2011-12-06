AnalyzerObservablePsiBarPsiChiralRightHandedCorrGauged::AnalyzerObservablePsiBarPsiChiralRightHandedCorrGauged(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservablePsiBarPsiCorrBase(fOps, aIOcon, SDreader, "PsiBarPsiCorrChiralRightHandedGauged", "pbpcxrg") { 
  fixGauge = false;
  randomGauge = false;
  projectorSelection = 1;
  multiplyWithPhiMatBSelection = 1;  
}


AnalyzerObservablePsiBarPsiChiralRightHandedCorrGauged::~AnalyzerObservablePsiBarPsiChiralRightHandedCorrGauged() {
}
