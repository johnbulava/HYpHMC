#include "EvaluateObservableTopBarTopPropagatorChiralLeftHandedGauged.h"

EvaluateObservableTopBarTopPropagatorChiralLeftHandedGauged::EvaluateObservableTopBarTopPropagatorChiralLeftHandedGauged(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservablePsiBarPsiPropagatorBase(aIOcon, sdr, obsWeight, obsDetSign, "TopBarTopPropagatorChiralLeftHandedGauged", "tbtpxlg", relStart, relEnd) { 
  fitThreshold_Lin = 0.5;
  kickOutPHatThreshold_Lin = 1.0;
  fitThreshold_FreeAna = 0.5;
  kickOutPHatThreshold_FreeAna = 1.0;  
  RescaleFactor = 0.5;
  
  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
}


EvaluateObservableTopBarTopPropagatorChiralLeftHandedGauged::~EvaluateObservableTopBarTopPropagatorChiralLeftHandedGauged() {
}

