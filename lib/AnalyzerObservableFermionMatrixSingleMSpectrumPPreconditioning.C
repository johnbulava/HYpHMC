#include "AnalyzerObservableFermionMatrixSingleMSpectrumPPreconditioning.h"

AnalyzerObservableFermionMatrixSingleMSpectrumPPreconditioning::AnalyzerObservableFermionMatrixSingleMSpectrumPPreconditioning(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "FermionMatrixSingleMSpectrumPPreconditioning", "specsmpp") { 
  if (fermiOps->get1DSizeLargest() >= 32) {
    analyzeEveryXXXconf = 8;
  }
  ini(getAnalyzerResultsCount());
}


AnalyzerObservableFermionMatrixSingleMSpectrumPPreconditioning::~AnalyzerObservableFermionMatrixSingleMSpectrumPPreconditioning() {
}


bool AnalyzerObservableFermionMatrixSingleMSpectrumPPreconditioning::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  double* phiField = phiFieldConf->getPhiFieldCopy(); 
  
  double Pm = phiFieldConf->getMagnetizationM();
  double Ps = phiFieldConf->getMagnetizationS();

  fermiOps->setPreconditioner(true, Pm, Ps);
  fermiOps->setQPreconditioner(false, 0.25, 0.25);
  fermiOps->setRPreconditioner(false, 1.0, 1.0);

  Complex* eigenvalues = fermiOps->calcFermionMatrixARPACKEigenValues(1, 20, phiField, 1E-1, false, NULL, false, true);

  for (int I=0; I<20; I++) {
    analyzerResults[2*I+0] = eigenvalues[I].x;
    analyzerResults[2*I+1] = eigenvalues[I].y;
  }
  
  delete[] eigenvalues;
  
  return true;
}


int AnalyzerObservableFermionMatrixSingleMSpectrumPPreconditioning::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservableFermionMatrixSingleMSpectrumPPreconditioning::getAnalyzerResultsCount() {
  return 40;
}
