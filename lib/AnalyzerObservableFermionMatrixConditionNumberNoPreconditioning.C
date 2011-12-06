AnalyzerObservableFermionMatrixConditionNumberNoPreconditioning::AnalyzerObservableFermionMatrixConditionNumberNoPreconditioning(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "FermionMatrixConditionNumberNoPreconditioning", "condnp") { 
  if (fermiOps->get1DSizeLargest() >= 32) {
    analyzeEveryXXXconf = 8;
  }
  ini(getAnalyzerResultsCount());
}


AnalyzerObservableFermionMatrixConditionNumberNoPreconditioning::~AnalyzerObservableFermionMatrixConditionNumberNoPreconditioning() {
}


bool AnalyzerObservableFermionMatrixConditionNumberNoPreconditioning::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  double* phiField = phiFieldConf->getPhiFieldCopy(); 
  double eigMin, eigMax, invCond;
  
  fermiOps->setPreconditioner(false, 1.0, 0.0);
  fermiOps->setQPreconditioner(false, 1.0, 1.0);
  fermiOps->setRPreconditioner(false, 1.0, 1.0);

  fermiOps->exactFermionMatrixConditionNumber(phiField, eigMin, eigMax, invCond, true, 1, false);

  analyzerResults[0] = eigMin;
  analyzerResults[1] = eigMax;
  analyzerResults[2] = invCond;
  
  return true;
}


int AnalyzerObservableFermionMatrixConditionNumberNoPreconditioning::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservableFermionMatrixConditionNumberNoPreconditioning::getAnalyzerResultsCount() {
  return 3;
}
