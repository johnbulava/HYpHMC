#include "EvaluateObservableTopBarTopChiralLeftHandedCondensatePointSource.h"

EvaluateObservableTopBarTopChiralLeftHandedCondensatePointSource::EvaluateObservableTopBarTopChiralLeftHandedCondensatePointSource(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservablePsiBarPsiCondensateBase(aIOcon, sdr, obsWeight, obsDetSign, "TopBarTopCondensateChiralLeftHandedPointSource", "tbtcondxlps", relStart, relEnd) { 
  snprintf(displayedFormula, 1000, "\\bar t_R t_L");

  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
}


EvaluateObservableTopBarTopChiralLeftHandedCondensatePointSource::~EvaluateObservableTopBarTopChiralLeftHandedCondensatePointSource() {
}
