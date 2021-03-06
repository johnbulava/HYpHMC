#include "AnalyzerObservableFermionMatrixSingleMFullNoXiSpectrumNoPreconditioning.h"

AnalyzerObservableFermionMatrixSingleMFullNoXiSpectrumNoPreconditioning::AnalyzerObservableFermionMatrixSingleMFullNoXiSpectrumNoPreconditioning(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "FermionMatrixSingleMFullNoXiSpectrumNoPreconditioning", "specsmfnxnp") { 
  L0 = fermiOps->get1DSizeL0(); 
  L1 = fermiOps->get1DSizeL1(); 
  L2 = fermiOps->get1DSizeL2(); 
  L3 = fermiOps->get1DSizeL3(); 

  if (L0*L1*L2*L3 >= 6*6*6*6) {
    analyzeEveryXXXconf = 4;
  }

  randomPhiField = false;
  
  ini(getAnalyzerResultsCount());
}


AnalyzerObservableFermionMatrixSingleMFullNoXiSpectrumNoPreconditioning::~AnalyzerObservableFermionMatrixSingleMFullNoXiSpectrumNoPreconditioning() {
}


bool AnalyzerObservableFermionMatrixSingleMFullNoXiSpectrumNoPreconditioning::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  double* phiField = phiFieldConf->getPhiFieldCopy(); 
  ComplexMatrix calcMatrix(1);
  double Parameter_MassSplit = fermiOps->getMassSplitRatio();
  
  if (randomPhiField) {
    if (LogLevel > 1) printf("WARNING: Random Phi-Field!!!\n");
    phiFieldConf->randomizeGaussHiggsField(phiField); 
  }
  
  fermiOps->constructNeubergerWithOutXiFermionMatrix(calcMatrix, (vector4D*)phiField, Parameter_MassSplit);
  bool b = calcMatrix.calcEigenvalues();

  for (int I=0; I<calcMatrix.matrixSize; I++) {
    analyzerResults[2*I+0] = calcMatrix.eigenvalues[I].x;
    analyzerResults[2*I+1] = calcMatrix.eigenvalues[I].y;
  }
  
  return b;
}


int AnalyzerObservableFermionMatrixSingleMFullNoXiSpectrumNoPreconditioning::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservableFermionMatrixSingleMFullNoXiSpectrumNoPreconditioning::getAnalyzerResultsCount() {
  return 16*(L0*L1*L2*L3);
}
