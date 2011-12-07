#include "AnalyzerObservableHiggsGoldstoneUnrotatedCorrelator.h"


AnalyzerObservableHiggsGoldstoneUnrotatedCorrelator::AnalyzerObservableHiggsGoldstoneUnrotatedCorrelator(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "HiggsGoldstoneUnrotatedCorrelator", "hgucorr") { 
  ini(getAnalyzerResultsCount());
}


AnalyzerObservableHiggsGoldstoneUnrotatedCorrelator::~AnalyzerObservableHiggsGoldstoneUnrotatedCorrelator() {
}


bool AnalyzerObservableHiggsGoldstoneUnrotatedCorrelator::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  double* phiField = phiFieldConf->getPhiFieldCopy(); 
  double rescale = SimParaSet->reparametrize_HiggsField(SimulationParameterSet_ContinuumNotation);
  phiFieldConf->multiplyHiggsFieldWithConst(phiField, rescale);
  
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
    int count = 0;
    for (ind[0]=0; ind[0]<L0; ind[0]++) {
      for (ind[1]=0; ind[1]<L1; ind[1]++) {
        for (ind[2]=0; ind[2]<L2; ind[2]++) {
          for (ind[3]=0; ind[3]<L3; ind[3]++) {
            analyzerResults[dirLoopCount*4*LargestL + 4*ind[timeDirection]+0] += phiField[4*count+0];
            analyzerResults[dirLoopCount*4*LargestL + 4*ind[timeDirection]+1] += phiField[4*count+1];
            analyzerResults[dirLoopCount*4*LargestL + 4*ind[timeDirection]+2] += phiField[4*count+2];
            analyzerResults[dirLoopCount*4*LargestL + 4*ind[timeDirection]+3] += phiField[4*count+3];

            count++; 
	  }
        }
      }
    }
    timeDirection++;
  }

  for (int I=0; I<getAnalyzerResultsCount(); I++) {
    analyzerResults[I] /= L0*L1*L2*L3;
    analyzerResults[I] *= LargestL;
  }

  return true;
}


int AnalyzerObservableHiggsGoldstoneUnrotatedCorrelator::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservableHiggsGoldstoneUnrotatedCorrelator::getAnalyzerResultsCount() {
  int LargestL = fermiOps->get1DSizeLargest();
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  
  if ((L0==L1) && (L0==L2) && (L0==L3)) {
    return 4*4*LargestL;
  } else {
    return 4*LargestL;
  }
}
