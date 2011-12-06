AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator::AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "HiggsGoldstoneUnrotatedTwoParticleK012Correlator", "hgu2pk012corr") { 
  ini(getAnalyzerResultsCount());
}


AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator::~AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator() {
}


bool AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  double* phiField = phiFieldConf->getPhiFieldCopy(); 
  double rescale = SimParaSet->reparametrize_HiggsField(SimulationParameterSet_ContinuumNotation);
  phiFieldConf->multiplyHiggsFieldWithConst(phiField, rescale);
  Complex* phiMomentumBuffer = phiFieldConf->performFourierTransform(phiField, true, 4);  

  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  int LargestL = fermiOps->get1DSizeLargest();
  for (int I=0; I<getAnalyzerResultsCount(); I++) {
    analyzerResults[I] = 0;
  }
  int relOffset[4];
  relOffset[3] = 4;
  relOffset[2] = relOffset[3] * L3;
  relOffset[1] = relOffset[2] * L2;
  relOffset[0] = relOffset[1] * L1;


  int timeDirection = fermiOps->getTimeDirection();
  int dirLoopMax = 1;
  if ((L0==L1) && (L0==L2) && (L0==L3)) {
    timeDirection = 0;
    dirLoopMax = 4;
  }
  
  for (int dirLoopCount=0; dirLoopCount<dirLoopMax; dirLoopCount++) {
    int count = 0;
    Complex* modes[4][3][3];
    double normFac = 1.0 / (L0*L1*L2*L3);
    for (int dir=0; dir<4; dir++) {
      if (dir != timeDirection) {
        for (int type=0; type<4; type++) {
          for (int k=0; k<3; k++) {
  	    int baseIndex = k*relOffset[dir];
   	    modes[type][count][k] = new Complex[LargestL];
            for (int It=0; It<LargestL; It++) {
  	      modes[type][count][k][It].x = 0;
	      modes[type][count][k][It].y = 0;
              for (int Ip=0; Ip<LargestL; Ip++) {
  	        modes[type][count][k][It] = modes[type][count][k][It] + exp((2*It*Ip*pi/LargestL) * ComplexI) * phiMomentumBuffer[type + baseIndex + Ip*relOffset[timeDirection]];
	      }
  	      modes[type][count][k][It] = normFac * modes[type][count][k][It];	    
	    }
 	  }
        }
        count++;
      }
    }  

    count = 0;
    for (int It=0; It<LargestL; It++) {
      for (int I=0; I<4; I++) {
        analyzerResults[dirLoopCount*52*LargestL + count] = modes[I][0][0][It].x;
        count++;
      }
      for (int k=1; k<3; k++) {
        for (int dir=0; dir<3; dir++) {
          for (int type=0; type<4; type++) {
            analyzerResults[dirLoopCount*52*LargestL + count+0] = modes[type][dir][k][It].x;
            analyzerResults[dirLoopCount*52*LargestL + count+1] = modes[type][dir][k][It].y;
  	    count+=2;
          }
        }
      }
    }
  
    for (int k=0; k<3; k++) {
      for (int dir=0; dir<3; dir++) {
        for (int type=0; type<4; type++) {
          delete[] modes[type][dir][k];
        }
      }
    }
    timeDirection++;
  }

  return true;
}


int AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator::getAnalyzerResultsCount() {
  int LargestL = fermiOps->get1DSizeLargest();
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  
  if ((L0==L1) && (L0==L2) && (L0==L3)) {
    return 4*52*LargestL;
  } else {
    return 52*LargestL;
  }
}
