AnalyzerObservableHiggsDefinedOnTimeSliceCorrelator::AnalyzerObservableHiggsDefinedOnTimeSliceCorrelator(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "HiggsDefinedOnTimeSliceCorrelator", "hdtcorr") { 
  ini(getAnalyzerResultsCount());
}


AnalyzerObservableHiggsDefinedOnTimeSliceCorrelator::~AnalyzerObservableHiggsDefinedOnTimeSliceCorrelator() {
}


bool AnalyzerObservableHiggsDefinedOnTimeSliceCorrelator::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {   
  double* phiField = phiFieldConf->getPhiFieldCopy(); 
  double rescale = SimParaSet->reparametrize_HiggsField(SimulationParameterSet_ContinuumNotation);
  phiFieldConf->multiplyHiggsFieldWithConst(phiField, rescale);  
  
  //Higgs-Modes
  int ind[4];
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  int LargestL = fermiOps->get1DSizeLargest();
  for (int I=0; I<getAnalyzerResultsCount(); I++) {
    analyzerResults[I] = 0;
  }

  int timeDirection = fermiOps->getTimeDirection();
  int dirLoopMax = 1;
  if ((L0==L1) && (L0==L2) && (L0==L3)) {
    timeDirection = 0;
    dirLoopMax = 4;
  }
  
  for (int dirLoopCount=0; dirLoopCount<dirLoopMax; dirLoopCount++) {
    vector4D* averageVector;
    averageVector = new vector4D[LargestL];
    for (int I=0; I<LargestL; I++) {
      averageVector[I][0] = 0;
      averageVector[I][1] = 0;
      averageVector[I][2] = 0;
      averageVector[I][3] = 0;
    }

    int count = 0;
    for (ind[0]=0; ind[0]<L0; ind[0]++) {
      for (ind[1]=0; ind[1]<L1; ind[1]++) {
        for (ind[2]=0; ind[2]<L2; ind[2]++) {
          for (ind[3]=0; ind[3]<L3; ind[3]++) {	
            averageVector[ind[timeDirection]][0] += phiField[4*count+0];
            averageVector[ind[timeDirection]][1] += phiField[4*count+1];
            averageVector[ind[timeDirection]][2] += phiField[4*count+2];
            averageVector[ind[timeDirection]][3] += phiField[4*count+3];

            count++;
  	  }
        }
      }
    }
    for (int I=0; I<LargestL; I++) {
      double norm  = sqr(averageVector[I][0]);
      norm += sqr(averageVector[I][1]);
      norm += sqr(averageVector[I][2]);
      norm += sqr(averageVector[I][3]);
      norm = sqrt(norm);
    
      averageVector[I][0] /= norm;
      averageVector[I][1] /= norm;
      averageVector[I][2] /= norm;
      averageVector[I][3] /= norm;  
    }
  
    count = 0;
    for (ind[0]=0; ind[0]<L0; ind[0]++) {
      for (ind[1]=0; ind[1]<L1; ind[1]++) {
        for (ind[2]=0; ind[2]<L2; ind[2]++) {
          for (ind[3]=0; ind[3]<L3; ind[3]++) {
            double scalar = 0;
            scalar += phiField[4*count+0]*averageVector[ind[timeDirection]][0];
            scalar += phiField[4*count+1]*averageVector[ind[timeDirection]][1];
            scalar += phiField[4*count+2]*averageVector[ind[timeDirection]][2];
            scalar += phiField[4*count+3]*averageVector[ind[timeDirection]][3];
    
            analyzerResults[dirLoopCount*LargestL + ind[timeDirection]] += scalar;

            count++;
	  }
        }
      }
    }
    timeDirection++;
    delete[] averageVector;    
  }

  for (int I=0; I<getAnalyzerResultsCount(); I++) {
    analyzerResults[I] /= L0*L1*L2*L3;
    analyzerResults[I] *= LargestL;
  }

  return true;  
}


int AnalyzerObservableHiggsDefinedOnTimeSliceCorrelator::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservableHiggsDefinedOnTimeSliceCorrelator::getAnalyzerResultsCount() {
  int LargestL = fermiOps->get1DSizeLargest();
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  
  if ((L0==L1) && (L0==L2) && (L0==L3)) {
    return 4*LargestL;
  } else {
    return LargestL;
  }
}
