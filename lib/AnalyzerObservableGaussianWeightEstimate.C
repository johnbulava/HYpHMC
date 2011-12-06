AnalyzerObservableGaussianWeightEstimate::AnalyzerObservableGaussianWeightEstimate(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "GaussianWeightEstimate", "gwest") { 
  ini(getAnalyzerResultsCount());
  pHMCProp = NULL;
}


AnalyzerObservableGaussianWeightEstimate::~AnalyzerObservableGaussianWeightEstimate() {
  delete pHMCProp;
}


bool AnalyzerObservableGaussianWeightEstimate::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  if (pHMCProp==NULL) {
    double* polEps = new double[5];
    double* polLam = new double[5];
    int* polDeg = new int[5];
    
    polEps[0] = SDReader->getPolynomLowerBound_P0();
    polEps[1] = SDReader->getPolynomLowerBound_P1();
    polEps[2] = SDReader->getPolynomLowerBound_P2();
    polEps[3] = SDReader->getPolynomLowerBound_P3();
    polEps[4] = SDReader->getPolynomLowerBound_P4();

    polLam[0] = SDReader->getPolynomUpperBound();
    polLam[1] = SDReader->getPolynomUpperBound();
    polLam[2] = SDReader->getPolynomUpperBound();
    polLam[3] = SDReader->getPolynomUpperBound();
    polLam[4] = SDReader->getPolynomUpperBound();
    
    polDeg[0] = SDReader->getPolynomDegree_P0();
    polDeg[1] = SDReader->getPolynomDegree_P1();
    polDeg[2] = SDReader->getPolynomDegree_P2();
    polDeg[3] = SDReader->getPolynomDegree_P3();
    polDeg[4] = SDReader->getPolynomDegree_P4();
    
    pHMCProp = new pHMCPropagator(fermiOps, SDReader->getLambda(), SDReader->getKappa(), SDReader->getExternalCurrent(), 
                                  SDReader->getModelParameterC6(), SDReader->getModelParameterC8(), SDReader->getModelParameterC10(),
                                  SDReader->getModelParameterLambda6(), SDReader->getModelParameterLambda8(), SDReader->getModelParameterLambda10(),
                                  SDReader->getNf(), 1.0, 
                                  SDReader->getSphericalHiggsIntegrationMode(), SDReader->getZetaForHiggsIntegrationMode(),
                                  SDReader->getTheta(),
                                  SDReader->getSubPolynomCount(), polEps, polLam, polDeg, 0, NULL,
                                  SDReader->getPolynomDigits(), SDReader->getAlpha(), SDReader->getMaximumPolynomDegreePerNode(),10);

    bool quasiHermiteanMode = SDReader->getUseQHM();				  
    pHMCProp->setPhiForceFourierType(0, 0);  
    pHMCProp->setOmegaMassAdaptionMode(0);
    pHMCProp->resetExactMMdagInverseSQRTOmegaAction();
    pHMCProp->synchronizedChangeOfQuasiHermiteanMode(quasiHermiteanMode);
    pHMCProp->synchronizedChangeOfModelSelection(1);
    pHMCProp->getNodesReady();
  }


  double* phiField = phiFieldConf->getPhiFieldCopy(); 
  
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

  pHMCForce* force = pHMCProp->getForce(0);
  int avgNcg = 0;
  int avgNmmdagApplications = 0;
  double avgWeight = 0;
  double avgWeightSigma = 0;
  int avgNmmdagApplicationsNewTechnique = 0;
  int rep = 10;
  
  for (int I=0; I<rep; I++) {
    int Ncg = 0;
    int NmmdagApplications = 0;
    force->sampleOmegaFieldsPurelyGaussian();
    double dummy = force->calcGaussianWeightFactor(phiField, 1E-10, Ncg, NmmdagApplications);
    avgWeight += dummy;
    avgWeightSigma += dummy*dummy;    
    avgNcg += Ncg;
    avgNmmdagApplications += NmmdagApplications;

    int neededIter;
    force->calcExactMMdagInverseSQRTomegaAction(phiField, 0.5, neededIter);
    avgNmmdagApplicationsNewTechnique += neededIter;
  }
  avgWeight /= rep;
  avgWeightSigma = sqrt(avgWeightSigma/rep - sqr(avgWeight));  
  avgNcg /= rep;
  avgNmmdagApplications /= rep;
  avgNmmdagApplicationsNewTechnique /= rep;

  analyzerResults[0] = avgWeight;
  analyzerResults[1] = avgWeightSigma;
  analyzerResults[2] = rep;
  analyzerResults[3] = avgNcg;
  analyzerResults[4] = avgNmmdagApplications;
  analyzerResults[5] = avgNmmdagApplicationsNewTechnique;
  
  return true;
}


int AnalyzerObservableGaussianWeightEstimate::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservableGaussianWeightEstimate::getAnalyzerResultsCount() {
  return 6;
}
