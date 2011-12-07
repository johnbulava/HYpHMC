#include "EvaluateObservableFermionMatrixSpectrumNoPreconditioning.h"

EvaluateObservableFermionMatrixSpectrumNoPreconditioning::EvaluateObservableFermionMatrixSpectrumNoPreconditioning(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) :  EvaluateObservableFermionMatrixSpectrumBase(aIOcon, sdr, "FermionMatrixSingleMSpectrumNoPreconditioning", "specsmnp", obsWeight, obsDetSign, relStart, relEnd) {
  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
}


EvaluateObservableFermionMatrixSpectrumNoPreconditioning::~EvaluateObservableFermionMatrixSpectrumNoPreconditioning() {
}


int EvaluateObservableFermionMatrixSpectrumNoPreconditioning::getAnalyzerResultsCount() {
  return 40;
}
