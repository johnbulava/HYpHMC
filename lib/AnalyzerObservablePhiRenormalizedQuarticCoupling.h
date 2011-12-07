#ifndef AnalyzerObservablePhiRenormalizedQuarticCoupling_included
#define AnalyzerObservablePhiRenormalizedQuarticCoupling_included

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


class AnalyzerObservablePhiRenormalizedQuarticCoupling : public AnalyzerObservable {
private:

public:
  AnalyzerObservablePhiRenormalizedQuarticCoupling(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservablePhiRenormalizedQuarticCoupling();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#endif
