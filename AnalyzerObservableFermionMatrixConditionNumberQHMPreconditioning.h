#ifndef AnalyzerObservableFermionMatrixConditionNumberQHMPreconditioning_included
#define AnalyzerObservableFermionMatrixConditionNumberQHMPreconditioning_included

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


class AnalyzerObservableFermionMatrixConditionNumberQHMPreconditioning : public AnalyzerObservable {
private:


public:
  AnalyzerObservableFermionMatrixConditionNumberQHMPreconditioning(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableFermionMatrixConditionNumberQHMPreconditioning();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#include "AnalyzerObservableFermionMatrixConditionNumberQHMPreconditioning.C"

#endif
