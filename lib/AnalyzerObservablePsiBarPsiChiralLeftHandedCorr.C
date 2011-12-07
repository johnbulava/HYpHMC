#include "AnalyzerObservablePsiBarPsiChiralLeftHandedCorr.h"

AnalyzerObservablePsiBarPsiChiralLeftHandedCorr::AnalyzerObservablePsiBarPsiChiralLeftHandedCorr(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservablePsiBarPsiCorrBase(fOps, aIOcon, SDreader, "PsiBarPsiCorrChiralLeftHanded", "pbpcxl") { 
  fixGauge = false;
  randomGauge = false;
  projectorSelection = -1;
  multiplyWithPhiMatBSelection = 0;  
}


AnalyzerObservablePsiBarPsiChiralLeftHandedCorr::~AnalyzerObservablePsiBarPsiChiralLeftHandedCorr() {
}
