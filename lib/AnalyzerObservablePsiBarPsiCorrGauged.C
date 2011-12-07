#include "AnalyzerObservablePsiBarPsiCorrGauged.h"

AnalyzerObservablePsiBarPsiCorrGauged::AnalyzerObservablePsiBarPsiCorrGauged(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservablePsiBarPsiCorrBase(fOps, aIOcon, SDreader, "PsiBarPsiCorrGauged", "pbpcg") { 
  fixGauge = true;
  randomGauge = false;
  projectorSelection = 0;
  multiplyWithPhiMatBSelection = 0;  
}


AnalyzerObservablePsiBarPsiCorrGauged::~AnalyzerObservablePsiBarPsiCorrGauged() {
}

