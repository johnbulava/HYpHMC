#ifndef AnalyzerObservableHiggsGoldstoneUnrotatedCorrelator_included
#define AnalyzerObservableHiggsGoldstoneUnrotatedCorrelator_included

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


class AnalyzerObservableHiggsGoldstoneUnrotatedCorrelator : public AnalyzerObservable {
private:


public:
  AnalyzerObservableHiggsGoldstoneUnrotatedCorrelator(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableHiggsGoldstoneUnrotatedCorrelator();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#endif
