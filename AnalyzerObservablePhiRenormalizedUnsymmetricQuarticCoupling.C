AnalyzerObservablePhiRenormalizedUnsymmetricQuarticCoupling::AnalyzerObservablePhiRenormalizedUnsymmetricQuarticCoupling(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "PhiRenormalizedUnsymmetricQuarticCoupling", "phirulam") { 
  int L0 = fOps->get1DSizeL0();
  int L1 = fOps->get1DSizeL1();
  int L2 = fOps->get1DSizeL2();
  int L3 = fOps->get1DSizeL3();
  
  smallL0 = L0/8;
  smallL1 = L1/8;
  smallL2 = L2/8;
  smallL3 = L3/8;
  if (smallL0<=2) smallL0 = 2;
  if (smallL1<=2) smallL1 = 2;
  if (smallL2<=2) smallL2 = 2;
  if (smallL3<=2) smallL3 = 2;

  ini(getAnalyzerResultsCount());
}


AnalyzerObservablePhiRenormalizedUnsymmetricQuarticCoupling::~AnalyzerObservablePhiRenormalizedUnsymmetricQuarticCoupling() {
}


bool AnalyzerObservablePhiRenormalizedUnsymmetricQuarticCoupling::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  double* phiField = phiFieldConf->getPhiFieldCopy(); 
  double rescale = SimParaSet->reparametrize_HiggsField(SimulationParameterSet_ContinuumNotation);
  phiFieldConf->multiplyHiggsFieldWithConst(phiField, rescale);

  Complex* phiMomentumBuffer = phiFieldConf->performFourierTransform(phiField, true, 4);  

  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
   
  //Phi - Modes
  int count = 0;
  double normFac = 1.0 / sqrt(L0*L1*L2*L3);
  for (int I0=1-smallL0; I0<smallL0; I0++) {
    int p0 = (I0+L0) % L0;  
    for (int I1=1-smallL1; I1<smallL1; I1++) {
      int p1 = (I1+L1) % L1;    
      for (int I2=1-smallL2; I2<smallL2; I2++) {
        int p2 = (I2+L2) % L2;      
        for (int I3=1-smallL3; I3<smallL3; I3++) {
          int p3 = (I3+L3) % L3;
	  int index = p3 + p2*L3 + p1*L2*L3 + p0*L1*L2*L3;
	  index *= 4;
	  for (int i=0; i<4; i++) {
  	    analyzerResults[count+0] = normFac*phiMomentumBuffer[index+i].x;
	    analyzerResults[count+1] = normFac*phiMomentumBuffer[index+i].y;
	    count += 2;
	  }
	}
      }
    }
  }

  return true;
}


int AnalyzerObservablePhiRenormalizedUnsymmetricQuarticCoupling::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservablePhiRenormalizedUnsymmetricQuarticCoupling::getAnalyzerResultsCount() {
  return 8*(2*smallL0-1)*(2*smallL1-1)*(2*smallL2-1)*(2*smallL3-1);
}
