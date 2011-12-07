#include "AnalyzerObservableWeight.h"

AnalyzerObservableWeight::AnalyzerObservableWeight(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "Weight", "weight") { 
  ini(getAnalyzerResultsCount());
}


AnalyzerObservableWeight::~AnalyzerObservableWeight() {
}


bool AnalyzerObservableWeight::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  if (phiFieldConf->isWeightAvail()) {
    double weight = phiFieldConf->getWeight();
    analyzerResults[0] = weight;
    return true;
  } else {
    return false;
  }
}


int AnalyzerObservableWeight::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservableWeight::getAnalyzerResultsCount() {
  return 1;
}
