#include "AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorr.h"

AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorr::AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorr(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase(fOps, aIOcon, SDreader, "TwoFermionPsiBarPsiChiralLeftHanded", "2fpbpcxl") { 
  fixGauge = false;
  randomGauge = false;
  projectorSelection = -1;
  multiplyWithPhiMatBSelection = false;  
}


AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorr::~AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorr() {
}
