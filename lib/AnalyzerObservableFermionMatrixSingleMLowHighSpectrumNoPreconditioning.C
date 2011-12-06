AnalyzerObservableFermionMatrixSingleMLowHighSpectrumNoPreconditioning::AnalyzerObservableFermionMatrixSingleMLowHighSpectrumNoPreconditioning(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "FermionMatrixSingleMLowHighSpectrumNoPreconditioning", "lhspecsmnp") { 
  if (fermiOps->get1DSizeLargest() >= 32) {
    analyzeEveryXXXconf = 8;
  }
  ini(getAnalyzerResultsCount());
}


AnalyzerObservableFermionMatrixSingleMLowHighSpectrumNoPreconditioning::~AnalyzerObservableFermionMatrixSingleMLowHighSpectrumNoPreconditioning() {
}


bool AnalyzerObservableFermionMatrixSingleMLowHighSpectrumNoPreconditioning::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  double* phiField = phiFieldConf->getPhiFieldCopy(); 
  
  double rho, r;  
  fermiOps->getDiracParameters(rho, r);
  
  fermiOps->setPreconditioner(false, 1.0, 0.0);
  fermiOps->setQPreconditioner(false, 0.25, 0.25);
  fermiOps->setRPreconditioner(false, 1.0, 1.0);
  
  Complex* eigenvalues = NULL;

  for (int I2=0; I2<6; I2++) {
    eigenvalues = fermiOps->calcFermionMatrixARPACKEigenValues(I2, 20, phiField, 0.1, false, NULL, false, true);

    for (int I=0; I<20; I++) {
      analyzerResults[40*I2+2*I+0] = (-1.0/(2*rho)) * eigenvalues[I].x;
      analyzerResults[40*I2+2*I+1] = (-1.0/(2*rho)) * eigenvalues[I].y;
    }  
  
    delete[] eigenvalues;
  }
  
  return true;
}


int AnalyzerObservableFermionMatrixSingleMLowHighSpectrumNoPreconditioning::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservableFermionMatrixSingleMLowHighSpectrumNoPreconditioning::getAnalyzerResultsCount() {
  return 240;
}
