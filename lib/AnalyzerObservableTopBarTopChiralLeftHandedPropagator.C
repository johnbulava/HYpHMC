#include "AnalyzerObservableTopBarTopChiralLeftHandedPropagator.h"

AnalyzerObservableTopBarTopChiralLeftHandedPropagator::AnalyzerObservableTopBarTopChiralLeftHandedPropagator(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservablePsiBarPsiPropagatorBase(fOps, aIOcon, SDreader, "TopBarTopPropagatorChiralLeftHandedGauged", "tbtpxlg") { 
  fixGauge = true;
  randomGauge = false;
  projectorSelection = -1;
  multiplyWithPhiMatBSelection = 0;  
  PsiPsiBarMatrixInd1Start = 0;
  PsiPsiBarMatrixInd2Start = 0;
  PsiPsiBarMatrixInd1End = 3;
  PsiPsiBarMatrixInd2End = 3;
  tresSumStartIndex = 0;
  tresSumEndIndex = 3;  
}


AnalyzerObservableTopBarTopChiralLeftHandedPropagator::~AnalyzerObservableTopBarTopChiralLeftHandedPropagator() {
}
