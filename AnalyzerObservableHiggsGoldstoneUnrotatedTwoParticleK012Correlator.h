#ifndef AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator_included
#define AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator_included

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


class AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator : public AnalyzerObservable {
private:


public:
  AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#include "AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator.C"

#endif
