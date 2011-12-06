#ifndef AnalyzerObservableFermionMatrixSingleMLowHighSpectrumNoPreconditioning_included
#define AnalyzerObservableFermionMatrixSingleMLowHighSpectrumNoPreconditioning_included

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


class AnalyzerObservableFermionMatrixSingleMLowHighSpectrumNoPreconditioning : public AnalyzerObservable {
private:


public:
  AnalyzerObservableFermionMatrixSingleMLowHighSpectrumNoPreconditioning(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableFermionMatrixSingleMLowHighSpectrumNoPreconditioning();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#include "AnalyzerObservableFermionMatrixSingleMLowHighSpectrumNoPreconditioning.C"

#endif
