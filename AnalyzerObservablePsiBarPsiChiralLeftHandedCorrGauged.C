AnalyzerObservablePsiBarPsiChiralLeftHandedCorrGauged::AnalyzerObservablePsiBarPsiChiralLeftHandedCorrGauged(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservablePsiBarPsiCorrBase(fOps, aIOcon, SDreader, "PsiBarPsiCorrChiralLeftHandedGauged", "pbpcxlg") { 
  fixGauge = false;
  randomGauge = false;
  projectorSelection = -1;
  multiplyWithPhiMatBSelection = -1;  
}


AnalyzerObservablePsiBarPsiChiralLeftHandedCorrGauged::~AnalyzerObservablePsiBarPsiChiralLeftHandedCorrGauged() {
}
