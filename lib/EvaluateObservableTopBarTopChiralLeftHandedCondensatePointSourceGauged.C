#include "EvaluateObservableTopBarTopChiralLeftHandedCondensatePointSourceGauged.h"


EvaluateObservableTopBarTopChiralLeftHandedCondensatePointSourceGauged::EvaluateObservableTopBarTopChiralLeftHandedCondensatePointSourceGauged(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservablePsiBarPsiCondensateBase(aIOcon, sdr, obsWeight, obsDetSign, "TopBarTopCondensateChiralLeftHandedPointSourceGauged", "tbtcondxlpsg", relStart, relEnd) { 
  snprintf(displayedFormula, 1000, "\\bar t_R t_L");

  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
}


EvaluateObservableTopBarTopChiralLeftHandedCondensatePointSourceGauged::~EvaluateObservableTopBarTopChiralLeftHandedCondensatePointSourceGauged() {
}
