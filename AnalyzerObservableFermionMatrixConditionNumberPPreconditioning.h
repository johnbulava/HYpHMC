#ifndef AnalyzerObservableFermionMatrixConditionNumberPPreconditioning_included
#define AnalyzerObservableFermionMatrixConditionNumberPPreconditioning_included

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


class AnalyzerObservableFermionMatrixConditionNumberPPreconditioning : public AnalyzerObservable {
private:


public:
  AnalyzerObservableFermionMatrixConditionNumberPPreconditioning(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableFermionMatrixConditionNumberPPreconditioning();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#include "AnalyzerObservableFermionMatrixConditionNumberPPreconditioning.C"

#endif
