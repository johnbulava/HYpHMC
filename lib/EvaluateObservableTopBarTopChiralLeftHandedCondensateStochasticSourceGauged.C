#include "EvaluateObservableTopBarTopChiralLeftHandedCondensateStochasticSourceGauged.h"

EvaluateObservableTopBarTopChiralLeftHandedCondensateStochasticSourceGauged::EvaluateObservableTopBarTopChiralLeftHandedCondensateStochasticSourceGauged(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservablePsiBarPsiCondensateBase(aIOcon, sdr, obsWeight, obsDetSign, "TopBarTopCondensateChiralLeftHandedStochSourceGauged", "tbtcondxlssg", relStart, relEnd) { 
  snprintf(displayedFormula, 1000, "\\bar t_R t_L");

  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
}


EvaluateObservableTopBarTopChiralLeftHandedCondensateStochasticSourceGauged::~EvaluateObservableTopBarTopChiralLeftHandedCondensateStochasticSourceGauged() {
}
