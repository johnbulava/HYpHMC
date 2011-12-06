AnalyzerObservableBottomBarBottomChiralLeftHandedCorrGauged::AnalyzerObservableBottomBarBottomChiralLeftHandedCorrGauged(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservablePsiBarPsiCorrBase(fOps, aIOcon, SDreader, "BottomBarBottomCorrChiralLeftHandedGauged", "bbbcxlg") { 
  fixGauge = true;
  randomGauge = false;
  projectorSelection = -1;
  multiplyWithPhiMatBSelection = 0;  
  PsiPsiBarMatrixInd1Start = 4;
  PsiPsiBarMatrixInd2Start = 4;
  PsiPsiBarMatrixInd1End = 7;
  PsiPsiBarMatrixInd2End = 7;
  tresSumStartIndex = 4;
  tresSumEndIndex = 7;  
}


AnalyzerObservableBottomBarBottomChiralLeftHandedCorrGauged::~AnalyzerObservableBottomBarBottomChiralLeftHandedCorrGauged() {
}
