#ifndef AnalyzerObservableHiggsRenormalizedQuarticCoupling_included
#define AnalyzerObservableHiggsRenormalizedQuarticCoupling_included

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


class AnalyzerObservableHiggsRenormalizedQuarticCoupling : public AnalyzerObservable {
private:

public:
  AnalyzerObservableHiggsRenormalizedQuarticCoupling(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableHiggsRenormalizedQuarticCoupling();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#endif
