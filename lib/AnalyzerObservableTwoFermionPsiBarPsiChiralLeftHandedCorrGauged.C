#include "AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged.h"

AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged::AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase(fOps, aIOcon, SDreader, "TwoFermionPsiBarPsiChiralLeftHandedGauged", "2fpbpcxlg") { 
  fixGauge = true;
  randomGauge = false;
  projectorSelection = -1;
  multiplyWithPhiMatBSelection = false;  
}


AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged::~AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged() {
}
