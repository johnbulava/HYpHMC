#ifndef AnalyzerObservableFermionMatrixSingleMFullSpectrumNoPreconditioning_included
#define AnalyzerObservableFermionMatrixSingleMFullSpectrumNoPreconditioning_included

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


class AnalyzerObservableFermionMatrixSingleMFullSpectrumNoPreconditioning : public AnalyzerObservable {
private:
  int L0;
  int L1;
  int L2;
  int L3;
  bool randomPhiField;


public:
  AnalyzerObservableFermionMatrixSingleMFullSpectrumNoPreconditioning(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableFermionMatrixSingleMFullSpectrumNoPreconditioning();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#endif
