#include "EvaluateObservablePsiBarPhiPsiChiralLeftHandedCondensate.h"

EvaluateObservablePsiBarPhiPsiChiralLeftHandedCondensate::EvaluateObservablePsiBarPhiPsiChiralLeftHandedCondensate(AnalyzerIOControl* aIOcon, 
StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservablePsiBarPsiCondensateBase(aIOcon, sdr, obsWeight, obsDetSign, "PsiBarPhiPsiCondensateChiralLeftHanded", "pbppcondxl", relStart, relEnd) { 
  numberOfRepeatedMeasurements = 1;
  snprintf(displayedFormula, 1000, "\\bar\\psi_R\\phi\\psi_L");

  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
}


EvaluateObservablePsiBarPhiPsiChiralLeftHandedCondensate::~EvaluateObservablePsiBarPhiPsiChiralLeftHandedCondensate() {
}
