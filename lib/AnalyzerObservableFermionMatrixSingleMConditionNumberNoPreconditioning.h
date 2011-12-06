#ifndef AnalyzerObservableFermionMatrixSingleMConditionNumberNoPreconditioning_included
#define AnalyzerObservableFermionMatrixSingleMConditionNumberNoPreconditioning_included

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


class AnalyzerObservableFermionMatrixSingleMConditionNumberNoPreconditioning : public AnalyzerObservable {
private:


public:
  AnalyzerObservableFermionMatrixSingleMConditionNumberNoPreconditioning(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableFermionMatrixSingleMConditionNumberNoPreconditioning();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#include "AnalyzerObservableFermionMatrixSingleMConditionNumberNoPreconditioning.C"

#endif
