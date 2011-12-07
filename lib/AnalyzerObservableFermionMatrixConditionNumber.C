#include "AnalyzerObservableFermionMatrixConditionNumber.h"

AnalyzerObservableFermionMatrixConditionNumber::AnalyzerObservableFermionMatrixConditionNumber(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "FermionMatrixConditionNumber", "cond") { 
  if (fermiOps->get1DSizeLargest() >= 32) {
    analyzeEveryXXXconf = 8;
  }
  ini(getAnalyzerResultsCount());
}


AnalyzerObservableFermionMatrixConditionNumber::~AnalyzerObservableFermionMatrixConditionNumber() {
}


bool AnalyzerObservableFermionMatrixConditionNumber::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  double* phiField = phiFieldConf->getPhiFieldCopy(); 
  double eigMin, eigMax, invCond;
  
  bool quasiHermiteanMode = SDReader->getUseQHM();
  bool useP = SDReader->getUseP();
  bool useQ = SDReader->getUseQ();
  bool useR = SDReader->getUseR();
  double Pm = SDReader->getPPrecondParameterM();  
  double Ps = SDReader->getPPrecondParameterS();  
  double Qmu = SDReader->getQPrecondParameterMu();
  double Qbeta = SDReader->getQPrecondParameterBeta(); 
  double Rm = SDReader->getRPrecondParameterM(); 
  double Rf = SDReader->getRPrecondParameterF();  
  
  fermiOps->setPreconditioner(useP, Pm, Ps);
  fermiOps->setQPreconditioner(useQ, Qmu, Qbeta);
  fermiOps->setRPreconditioner(useR, Rm, Rf);

  fermiOps->exactFermionMatrixConditionNumber(phiField, eigMin, eigMax, invCond, true, 1, quasiHermiteanMode);

  analyzerResults[0] = eigMin;
  analyzerResults[1] = eigMax;
  analyzerResults[2] = invCond;
  
  return true;
}


int AnalyzerObservableFermionMatrixConditionNumber::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservableFermionMatrixConditionNumber::getAnalyzerResultsCount() {
  return 3;
}
