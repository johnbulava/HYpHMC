#include "EvaluateObservableTopBarTopCorrChiralLeftHandedGauged.h"

EvaluateObservableTopBarTopCorrChiralLeftHandedGauged::EvaluateObservableTopBarTopCorrChiralLeftHandedGauged(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservablePsiBarPsiCorrBase(aIOcon, sdr, obsWeight, obsDetSign, "TopBarTopCorrChiralLeftHandedGauged", "tbtcxlg", relStart, relEnd) { 
  int LargestL = getLargestL(SDReader);  
  massAnalyzer = new MassCorrelationMatrixAnalyzer(1, LargestL, true, false, true, false, false, 1, 100000, "Top");  
  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
}


EvaluateObservableTopBarTopCorrChiralLeftHandedGauged::~EvaluateObservableTopBarTopCorrChiralLeftHandedGauged() {
}
