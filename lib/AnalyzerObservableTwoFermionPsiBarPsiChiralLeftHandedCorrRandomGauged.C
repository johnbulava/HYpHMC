AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrRandomGauged::AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrRandomGauged(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase(fOps, aIOcon, SDreader, "TwoFermionPsiBarPsiChiralLeftHandedRandomGauged", "2fpbpcxlrg") { 
  fixGauge = false;
  randomGauge = true;
  projectorSelection = -1;
  multiplyWithPhiMatBSelection = false;  
}


AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrRandomGauged::~AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrRandomGauged() {
}
