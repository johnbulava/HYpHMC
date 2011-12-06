#ifndef EvaluateObservableFermionMatrixLowHighSpectrumNoPreconditioning_included
#define EvaluateObservableFermionMatrixLowHighSpectrumNoPreconditioning_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "EvaluateObservableFermionMatrixSpectrumBase.h"


class EvaluateObservableFermionMatrixLowHighSpectrumNoPreconditioning : public EvaluateObservableFermionMatrixSpectrumBase {
private:  

public:    
  EvaluateObservableFermionMatrixLowHighSpectrumNoPreconditioning(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableFermionMatrixLowHighSpectrumNoPreconditioning();

  int getAnalyzerResultsCount();
};


#include "EvaluateObservableFermionMatrixLowHighSpectrumNoPreconditioning.C"

#endif
