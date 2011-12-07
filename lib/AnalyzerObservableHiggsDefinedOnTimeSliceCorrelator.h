#ifndef AnalyzerObservableHiggsDefinedOnTimeSliceCorrelator_included
#define AnalyzerObservableHiggsDefinedOnTimeSliceCorrelator_included

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


class AnalyzerObservableHiggsDefinedOnTimeSliceCorrelator : public AnalyzerObservable {
private:


public:
  AnalyzerObservableHiggsDefinedOnTimeSliceCorrelator(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableHiggsDefinedOnTimeSliceCorrelator();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#endif
