#ifndef AnalyzerObservableFermionMatrixSingleMSpectrumNoPreconditioning_included
#define AnalyzerObservableFermionMatrixSingleMSpectrumNoPreconditioning_included

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


class AnalyzerObservableFermionMatrixSingleMSpectrumNoPreconditioning : public AnalyzerObservable {
private:


public:
  AnalyzerObservableFermionMatrixSingleMSpectrumNoPreconditioning(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableFermionMatrixSingleMSpectrumNoPreconditioning();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#endif
