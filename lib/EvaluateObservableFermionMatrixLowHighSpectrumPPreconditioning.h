#ifndef EvaluateObservableFermionMatrixLowHighSpectrumPPreconditioning_included
#define EvaluateObservableFermionMatrixLowHighSpectrumPPreconditioning_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "EvaluateObservableFermionMatrixSpectrumBase.h"


class EvaluateObservableFermionMatrixLowHighSpectrumPPreconditioning : public EvaluateObservableFermionMatrixSpectrumBase {
private:  

public:    
  EvaluateObservableFermionMatrixLowHighSpectrumPPreconditioning(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableFermionMatrixLowHighSpectrumPPreconditioning();

  int getAnalyzerResultsCount();
};


#include "EvaluateObservableFermionMatrixLowHighSpectrumPPreconditioning.C"

#endif
