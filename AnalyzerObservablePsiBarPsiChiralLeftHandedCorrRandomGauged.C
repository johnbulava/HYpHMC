AnalyzerObservablePsiBarPsiChiralLeftHandedCorrRandomGauged::AnalyzerObservablePsiBarPsiChiralLeftHandedCorrRandomGauged(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservablePsiBarPsiCorrBase(fOps, aIOcon, SDreader, "PsiBarPsiCorrChiralLeftHandedRandomGauged", "pbpcxlrg") { 
  fixGauge = false;
  randomGauge = true;
  projectorSelection = -1;
  multiplyWithPhiMatBSelection = 0;  
}


AnalyzerObservablePsiBarPsiChiralLeftHandedCorrRandomGauged::~AnalyzerObservablePsiBarPsiChiralLeftHandedCorrRandomGauged() {
}
