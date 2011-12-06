#ifndef AnalyzerObservableFermionMatrixConditionNumberNoPreconditioning_included
#define AnalyzerObservableFermionMatrixConditionNumberNoPreconditioning_included

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


class AnalyzerObservableFermionMatrixConditionNumberNoPreconditioning : public AnalyzerObservable {
private:


public:
  AnalyzerObservableFermionMatrixConditionNumberNoPreconditioning(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableFermionMatrixConditionNumberNoPreconditioning();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#include "AnalyzerObservableFermionMatrixConditionNumberNoPreconditioning.C"

#endif
