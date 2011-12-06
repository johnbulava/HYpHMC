AnalyzerObservableFermionMatrixSingleMLowHighSpectrumPPreconditioning::AnalyzerObservableFermionMatrixSingleMLowHighSpectrumPPreconditioning(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "FermionMatrixSingleMLowHighSpectrumPPreconditioning", "lhspecsmpp") { 
  if (fermiOps->get1DSizeLargest() >= 32) {
    analyzeEveryXXXconf = 8;
  }
  ini(getAnalyzerResultsCount());
}


AnalyzerObservableFermionMatrixSingleMLowHighSpectrumPPreconditioning::~AnalyzerObservableFermionMatrixSingleMLowHighSpectrumPPreconditioning() {
}


bool AnalyzerObservableFermionMatrixSingleMLowHighSpectrumPPreconditioning::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  double* phiField = phiFieldConf->getPhiFieldCopy(); 
  
  double Pm = phiFieldConf->getMagnetizationM();
  double Ps = phiFieldConf->getMagnetizationS();

  fermiOps->setPreconditioner(true, Pm, Ps);
  fermiOps->setQPreconditioner(false, 0.25, 0.25);
  fermiOps->setRPreconditioner(false, 1.0, 1.0);

  Complex* eigenvalues = NULL;
  
  for (int I2=0; I2<6; I2++) {
    eigenvalues = fermiOps->calcFermionMatrixARPACKEigenValues(I2, 20, phiField, 0.1, false, NULL, false, true);

    for (int I=0; I<20; I++) {
      analyzerResults[I2*40+2*I+0] = eigenvalues[I].x;
      analyzerResults[I2*40+2*I+1] = eigenvalues[I].y;
    }  
    
    delete[] eigenvalues;
  }
  
  return true;
}


int AnalyzerObservableFermionMatrixSingleMLowHighSpectrumPPreconditioning::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservableFermionMatrixSingleMLowHighSpectrumPPreconditioning::getAnalyzerResultsCount() {
  return 240;
}
