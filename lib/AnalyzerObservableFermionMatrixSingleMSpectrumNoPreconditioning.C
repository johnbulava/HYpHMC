#include "AnalyzerObservableFermionMatrixSingleMSpectrumNoPreconditioning.h"

AnalyzerObservableFermionMatrixSingleMSpectrumNoPreconditioning::AnalyzerObservableFermionMatrixSingleMSpectrumNoPreconditioning(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "FermionMatrixSingleMSpectrumNoPreconditioning", "specsmnp") { 
  if (fermiOps->get1DSizeLargest() >= 32) {
    analyzeEveryXXXconf = 8;
  }
  ini(getAnalyzerResultsCount());
}


AnalyzerObservableFermionMatrixSingleMSpectrumNoPreconditioning::~AnalyzerObservableFermionMatrixSingleMSpectrumNoPreconditioning() {
}


bool AnalyzerObservableFermionMatrixSingleMSpectrumNoPreconditioning::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  double* phiField = phiFieldConf->getPhiFieldCopy(); 
  
  double rho, r;  
  fermiOps->getDiracParameters(rho, r);
  
  fermiOps->setPreconditioner(false, 1.0, 0.0);
  fermiOps->setQPreconditioner(false, 0.25, 0.25);
  fermiOps->setRPreconditioner(false, 1.0, 1.0);

  Complex* eigenvalues = fermiOps->calcFermionMatrixARPACKEigenValues(0, 20, phiField, 1E-1, false, NULL, false, true);

  for (int I=0; I<20; I++) {
    analyzerResults[2*I+0] = (-1.0/(2*rho)) * eigenvalues[I].x;
    analyzerResults[2*I+1] = (-1.0/(2*rho)) * eigenvalues[I].y;
  }
  
  delete[] eigenvalues;
  
  return true;
}


int AnalyzerObservableFermionMatrixSingleMSpectrumNoPreconditioning::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservableFermionMatrixSingleMSpectrumNoPreconditioning::getAnalyzerResultsCount() {
  return 40;
}
