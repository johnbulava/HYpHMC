#ifndef AnalyzerObservablePhiRenormalizedUnsymmetricQuarticCoupling_included
#define AnalyzerObservablePhiRenormalizedUnsymmetricQuarticCoupling_included

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


class AnalyzerObservablePhiRenormalizedUnsymmetricQuarticCoupling : public AnalyzerObservable {
private:
  int smallL0;
  int smallL1;
  int smallL2;
  int smallL3;

public:
  AnalyzerObservablePhiRenormalizedUnsymmetricQuarticCoupling(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservablePhiRenormalizedUnsymmetricQuarticCoupling();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#include "AnalyzerObservablePhiRenormalizedUnsymmetricQuarticCoupling.C"

#endif
