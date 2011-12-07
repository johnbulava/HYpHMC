#include "AnalyzerObservableFermionMatrixConditionNumberPPreconditioning.h"

AnalyzerObservableFermionMatrixConditionNumberPPreconditioning::AnalyzerObservableFermionMatrixConditionNumberPPreconditioning(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "FermionMatrixConditionNumberPPreconditioning", "condpp") { 
  if (fermiOps->get1DSizeLargest() >= 32) {
    analyzeEveryXXXconf = 8;
  }
  ini(getAnalyzerResultsCount());
}


AnalyzerObservableFermionMatrixConditionNumberPPreconditioning::~AnalyzerObservableFermionMatrixConditionNumberPPreconditioning() {
}


bool AnalyzerObservableFermionMatrixConditionNumberPPreconditioning::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  double* phiField = phiFieldConf->getPhiFieldCopy(); 
  double eigMin, eigMax, invCond;
  
  double Pm = phiFieldConf->getMagnetizationM();
  double Ps = phiFieldConf->getMagnetizationS();
  
  fermiOps->setPreconditioner(true, Pm, Ps);
  fermiOps->setQPreconditioner(false, 1.0, 1.0);
  fermiOps->setRPreconditioner(false, 1.0, 1.0);

  fermiOps->exactFermionMatrixConditionNumber(phiField, eigMin, eigMax, invCond, true, 1, false);

  analyzerResults[0] = eigMin;
  analyzerResults[1] = eigMax;
  analyzerResults[2] = invCond;
  
  return true;
}


int AnalyzerObservableFermionMatrixConditionNumberPPreconditioning::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservableFermionMatrixConditionNumberPPreconditioning::getAnalyzerResultsCount() {
  return 3;
}
