#ifndef AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator_included
#define AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "Global.h"
#include "AnalyzerIOControl.h"
#include "FermionMatrixOperations.h"
#include "StateDescriptorReader.h"
#include "AnalyzerObservable.h"
#include "AnalyzerPhiFieldConfiguration.h"
#include "LatticeMomentumBins.h"


class AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator : public AnalyzerObservable {
private:


public:
  AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#include "AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator.C"

#endif
