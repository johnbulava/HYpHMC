#ifndef AnalyzerObservableFermionMatrixSingleMConditionNumberPPreconditioning_included
#define AnalyzerObservableFermionMatrixSingleMConditionNumberPPreconditioning_included

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


class AnalyzerObservableFermionMatrixSingleMConditionNumberPPreconditioning : public AnalyzerObservable {
private:


public:
  AnalyzerObservableFermionMatrixSingleMConditionNumberPPreconditioning(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableFermionMatrixSingleMConditionNumberPPreconditioning();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#include "AnalyzerObservableFermionMatrixSingleMConditionNumberPPreconditioning.C"

#endif
