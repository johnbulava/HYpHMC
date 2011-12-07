#include "EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged.h"

EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged::EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservablePsiBarPsiCorrBase(aIOcon, sdr, obsWeight, obsDetSign, "TwoFermionPsiBarPsiChiralLeftHandedGauged", "2fpbpcxlg", relStart, relEnd) { 
  operatorVEVdataAvail = true;
  int LargestL = getLargestL(SDReader);  
  massAnalyzer = new MassCorrelationMatrixAnalyzer(1, LargestL, true, false, true, true, true, 2, 100000, "Top");  
  fitMassCount = 2;  
  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);  
}


EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged::~EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged() {
}
