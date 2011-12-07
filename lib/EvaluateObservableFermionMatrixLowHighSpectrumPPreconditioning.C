#include "EvaluateObservableFermionMatrixLowHighSpectrumPPreconditioning.h"

EvaluateObservableFermionMatrixLowHighSpectrumPPreconditioning::EvaluateObservableFermionMatrixLowHighSpectrumPPreconditioning(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservableFermionMatrixSpectrumBase(aIOcon, sdr, "FermionMatrixSingleMLowHighSpectrumPPreconditioning", "lhspecsmpp", obsWeight, obsDetSign, relStart, relEnd) {
  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
}


EvaluateObservableFermionMatrixLowHighSpectrumPPreconditioning::~EvaluateObservableFermionMatrixLowHighSpectrumPPreconditioning() {
}


int EvaluateObservableFermionMatrixLowHighSpectrumPPreconditioning::getAnalyzerResultsCount() {
  return 240;
}
