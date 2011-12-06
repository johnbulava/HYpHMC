#ifndef AnalyzerObservableFermionMatrixSingleMLowHighSpectrumPPreconditioning_included
#define AnalyzerObservableFermionMatrixSingleMLowHighSpectrumPPreconditioning_included

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


class AnalyzerObservableFermionMatrixSingleMLowHighSpectrumPPreconditioning : public AnalyzerObservable {
private:


public:
  AnalyzerObservableFermionMatrixSingleMLowHighSpectrumPPreconditioning(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableFermionMatrixSingleMLowHighSpectrumPPreconditioning();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#include "AnalyzerObservableFermionMatrixSingleMLowHighSpectrumPPreconditioning.C"

#endif
