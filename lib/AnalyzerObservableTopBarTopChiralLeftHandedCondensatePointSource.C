#include "AnalyzerObservableTopBarTopChiralLeftHandedCondensatePointSource.h"

AnalyzerObservableTopBarTopChiralLeftHandedCondensatePointSource::AnalyzerObservableTopBarTopChiralLeftHandedCondensatePointSource(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservablePsiBarPsiCondensateBase(fOps, aIOcon, SDreader, "TopBarTopCondensateChiralLeftHandedPointSource", "tbtcondxlps") { 
  fixGauge = false;
  randomGauge = false;
  stochasticalSource = false;
  projectorSelection = -1;
  PsiPsiBarMatrixInd1Start = 0;
  PsiPsiBarMatrixInd2Start = 0;
  PsiPsiBarMatrixInd1End = 3;
  PsiPsiBarMatrixInd2End = 3;
  tresSumStartIndex = 0;
  tresSumEndIndex = 3;  
}


AnalyzerObservableTopBarTopChiralLeftHandedCondensatePointSource::~AnalyzerObservableTopBarTopChiralLeftHandedCondensatePointSource() {
}
