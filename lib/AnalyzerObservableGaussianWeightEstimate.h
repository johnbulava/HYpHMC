#ifndef AnalyzerObservableGaussianWeightEstimate_included
#define AnalyzerObservableGaussianWeightEstimate_included

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
#include "pHMCPropagator.h"


class AnalyzerObservableGaussianWeightEstimate : public AnalyzerObservable {
private:
  pHMCPropagator* pHMCProp;


public:
  AnalyzerObservableGaussianWeightEstimate(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableGaussianWeightEstimate();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#include "AnalyzerObservableGaussianWeightEstimate.C"

#endif
