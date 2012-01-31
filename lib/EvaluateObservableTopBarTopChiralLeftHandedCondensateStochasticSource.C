#include "EvaluateObservableTopBarTopChiralLeftHandedCondensateStochasticSource.h"

EvaluateObservableTopBarTopChiralLeftHandedCondensateStochasticSource::EvaluateObservableTopBarTopChiralLeftHandedCondensateStochasticSource(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservablePsiBarPsiCondensateBase(aIOcon, sdr, obsWeight, obsDetSign, "TopBarTopCondensateChiralLeftHandedStochSource", "tbtcondxlss", relStart, relEnd) { 
  snprintf(displayedFormula, 1000, "\\bar t_R t_L");

  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
}


EvaluateObservableTopBarTopChiralLeftHandedCondensateStochasticSource::~EvaluateObservableTopBarTopChiralLeftHandedCondensateStochasticSource() {
}
