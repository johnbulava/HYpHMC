#include "AnalyzerObservablePsiBarPhiPsiChiralLeftHandedCondensate.h"

AnalyzerObservablePsiBarPhiPsiChiralLeftHandedCondensate::AnalyzerObservablePsiBarPhiPsiChiralLeftHandedCondensate(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservablePsiBarPhiPsiCondensateBase(fOps, aIOcon, SDreader, "PsiBarPhiPsiCondensateChiralLeftHanded", "pbppcondxl") { 
  fixGauge = false;
  randomGauge = false;
  projectorSelection = -1;
  multiplyWithPhiMatBSelection = -1;
}


AnalyzerObservablePsiBarPhiPsiChiralLeftHandedCondensate::~AnalyzerObservablePsiBarPhiPsiChiralLeftHandedCondensate() {
}
