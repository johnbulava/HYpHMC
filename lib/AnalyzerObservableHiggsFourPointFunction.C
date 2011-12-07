#include "AnalyzerObservableHiggsFourPointFunction.h"

AnalyzerObservableHiggsFourPointFunction::AnalyzerObservableHiggsFourPointFunction(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "HiggsPropagator", "hprop") { 
  ini(getAnalyzerResultsCount());
}


AnalyzerObservableHiggsFourPointFunction::~AnalyzerObservableHiggsFourPointFunction() {
}


bool AnalyzerObservableHiggsFourPointFunction::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  vector4D averageVector;
  averageVector[0] = phiFieldConf->getAvgPhiFieldVectorComponent(0) / phiFieldConf->getAvgPhiFieldVectorLength();
  averageVector[1] = phiFieldConf->getAvgPhiFieldVectorComponent(1) / phiFieldConf->getAvgPhiFieldVectorLength();
  averageVector[2] = phiFieldConf->getAvgPhiFieldVectorComponent(2) / phiFieldConf->getAvgPhiFieldVectorLength();
  averageVector[3] = phiFieldConf->getAvgPhiFieldVectorComponent(3) / phiFieldConf->getAvgPhiFieldVectorLength();
    
  double* phiField = phiFieldConf->getPhiFieldCopy(); 
  double mag = phiFieldConf->getMagnetizationM();
  double rescale = SimParaSet->reparametrize_HiggsField(SimulationParameterSet_ContinuumNotation);
  phiFieldConf->multiplyHiggsFieldWithConst(phiField, rescale);
  double magRescaled = rescale*mag;
  
  
  //Higgs-Modes
  int ind[4];
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  int count = 0;
  for (ind[0]=0; ind[0]<L0; ind[0]++) {
    for (ind[1]=0; ind[1]<L1; ind[1]++) {
      for (ind[2]=0; ind[2]<L2; ind[2]++) {
        for (ind[3]=0; ind[3]<L3; ind[3]++) {
          double scalar = 0;
          scalar += phiField[4*count+0]*averageVector[0];
          scalar += phiField[4*count+1]*averageVector[1];
          scalar += phiField[4*count+2]*averageVector[2];
          scalar += phiField[4*count+3]*averageVector[3];
    
          phiField[4*count+0] = scalar - magRescaled;
          phiField[4*count+1] = NaN;
          phiField[4*count+2] = NaN;
          phiField[4*count+3] = NaN;

          count++;
	}
      }
    }
  }
  

  Complex* phiMomentumBuffer = phiFieldConf->performFourierTransform(phiField, true, 1);  
 
  
  //Higgs - General 4-Point Function for selected values of momentum p
  double normFac = 1.0 / (L0*L1*L2*L3);
  double* data = new double[L0*L1*L2*L3];
  for (int I=0; I<L0*L1*L2*L3; I++) {
    data[I] = sqr(phiMomentumBuffer[I].x) + sqr(phiMomentumBuffer[I].y);
    data[I] *= normFac;
  }  
  delete[] data;
  data = NULL;
  

  for (int I=0; I<latticeBins->getMomentumSqrSlotCount(); I++) {
    analyzerResults[I] = data[I];
  }

  return true;
}


int AnalyzerObservableHiggsFourPointFunction::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservableHiggsFourPointFunction::getAnalyzerResultsCount() {
  return 2;
}
