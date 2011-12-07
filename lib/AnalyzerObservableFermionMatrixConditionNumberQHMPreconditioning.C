#include "AnalyzerObservableFermionMatrixConditionNumberQHMPreconditioning.h"

AnalyzerObservableFermionMatrixConditionNumberQHMPreconditioning::AnalyzerObservableFermionMatrixConditionNumberQHMPreconditioning(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "FermionMatrixConditionNumberQHMPreconditioning", "condqhm") { 
  if (fermiOps->get1DSizeLargest() >= 32) {
    analyzeEveryXXXconf = 8;
  }
  ini(getAnalyzerResultsCount());
}


AnalyzerObservableFermionMatrixConditionNumberQHMPreconditioning::~AnalyzerObservableFermionMatrixConditionNumberQHMPreconditioning() {
}


bool AnalyzerObservableFermionMatrixConditionNumberQHMPreconditioning::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  double* phiField = phiFieldConf->getPhiFieldCopy(); 
  double eigMin, eigMax, invCond;
  
  double Rm = phiFieldConf->getMagnetizationM();
  double Rf = 1.0;
  
  fermiOps->setPreconditioner(false, 1.0, 0.0);
  fermiOps->setQPreconditioner(false, 1.0, 1.0);
  fermiOps->setRPreconditioner(true, Rm, Rf);

  fermiOps->exactFermionMatrixConditionNumber(phiField, eigMin, eigMax, invCond, true, 1, true);

  analyzerResults[0] = eigMin;
  analyzerResults[1] = eigMax;
  analyzerResults[2] = invCond;
  
  return true;
}


int AnalyzerObservableFermionMatrixConditionNumberQHMPreconditioning::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservableFermionMatrixConditionNumberQHMPreconditioning::getAnalyzerResultsCount() {
  return 3;
}
