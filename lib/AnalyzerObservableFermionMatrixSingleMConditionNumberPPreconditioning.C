#include "AnalyzerObservableFermionMatrixSingleMConditionNumberPPreconditioning.h"

AnalyzerObservableFermionMatrixSingleMConditionNumberPPreconditioning::AnalyzerObservableFermionMatrixSingleMConditionNumberPPreconditioning(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "FermionMatrixSingleMConditionNumberPPreconditioning", "condsmpp") { 
  if (fermiOps->get1DSizeLargest() >= 32) {
    analyzeEveryXXXconf = 8;
  }
  ini(getAnalyzerResultsCount());
}


AnalyzerObservableFermionMatrixSingleMConditionNumberPPreconditioning::~AnalyzerObservableFermionMatrixSingleMConditionNumberPPreconditioning() {
}


bool AnalyzerObservableFermionMatrixSingleMConditionNumberPPreconditioning::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  double* phiField = phiFieldConf->getPhiFieldCopy(); 
  double eigMin, eigMax, invCond;
  
  double Pm = phiFieldConf->getMagnetizationM();
  double Ps = phiFieldConf->getMagnetizationS();

  fermiOps->setPreconditioner(true, Pm, Ps);
  fermiOps->setQPreconditioner(false, 0.25, 0.25);
  fermiOps->setRPreconditioner(false, 1.0, 1.0);

  fermiOps->exactFermionMatrixConditionNumber(phiField, eigMin, eigMax, invCond, false, 4, false);

  analyzerResults[0] = eigMin;
  analyzerResults[1] = eigMax;
  analyzerResults[2] = invCond;
  
  return true;
}


int AnalyzerObservableFermionMatrixSingleMConditionNumberPPreconditioning::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservableFermionMatrixSingleMConditionNumberPPreconditioning::getAnalyzerResultsCount() {
  return 3;
}
