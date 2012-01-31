#include "AnalyzerObservableTopBarTopChiralLeftHandedCondensateStochasticSourceGauged.h"

AnalyzerObservableTopBarTopChiralLeftHandedCondensateStochasticSourceGauged::AnalyzerObservableTopBarTopChiralLeftHandedCondensateStochasticSourceGauged(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservablePsiBarPsiCondensateBase(fOps, aIOcon, SDreader, "TopBarTopCondensateChiralLeftHandedStochSourceGauged", "tbtcondxlssg") { 
  fixGauge = true;
  randomGauge = false;
  stochasticalSource = true;
  projectorSelection = -1;
  PsiPsiBarMatrixInd1Start = 0;
  PsiPsiBarMatrixInd2Start = 0;
  PsiPsiBarMatrixInd1End = 3;
  PsiPsiBarMatrixInd2End = 3;
  tresSumStartIndex = 0;
  tresSumEndIndex = 3;  
}


AnalyzerObservableTopBarTopChiralLeftHandedCondensateStochasticSourceGauged::~AnalyzerObservableTopBarTopChiralLeftHandedCondensateStochasticSourceGauged() {
}
