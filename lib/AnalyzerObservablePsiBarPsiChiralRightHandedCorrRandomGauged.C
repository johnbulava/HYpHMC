#include "AnalyzerObservablePsiBarPsiChiralRightHandedCorrRandomGauged.h"

AnalyzerObservablePsiBarPsiChiralRightHandedCorrRandomGauged::AnalyzerObservablePsiBarPsiChiralRightHandedCorrRandomGauged(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservablePsiBarPsiCorrBase(fOps, aIOcon, SDreader, "PsiBarPsiCorrChiralRightHandedRandomGauged", "pbpcxrrg") { 
  fixGauge = false;
  randomGauge = true;
  projectorSelection = 1;
  multiplyWithPhiMatBSelection = 0;  
}


AnalyzerObservablePsiBarPsiChiralRightHandedCorrRandomGauged::~AnalyzerObservablePsiBarPsiChiralRightHandedCorrRandomGauged() {
}
