AnalyzerObservableMagnetizations::AnalyzerObservableMagnetizations(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "Magnetizations", "mags") { 
  ini(getAnalyzerResultsCount());
}


AnalyzerObservableMagnetizations::~AnalyzerObservableMagnetizations() {
}


bool AnalyzerObservableMagnetizations::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  analyzerResults[0] = phiFieldConf->getMagnetizationM();
  analyzerResults[1] = phiFieldConf->getMagnetizationS();
  analyzerResults[2] = phiFieldConf->getPhiFieldAvgVectorLength();
  analyzerResults[3] = phiFieldConf->getPhiFieldAvgVectorLengthVariation();
  analyzerResults[4] = phiFieldConf->getAvgPhiFieldVectorComponent(0);
  analyzerResults[5] = phiFieldConf->getAvgPhiFieldVectorComponent(1);
  analyzerResults[6] = phiFieldConf->getAvgPhiFieldVectorComponent(2);
  analyzerResults[7] = phiFieldConf->getAvgPhiFieldVectorComponent(3);
  return true;
}


int AnalyzerObservableMagnetizations::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservableMagnetizations::getAnalyzerResultsCount() {
  return 8;
}
