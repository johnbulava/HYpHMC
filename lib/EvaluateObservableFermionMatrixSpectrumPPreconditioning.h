#ifndef EvaluateObservableFermionMatrixSpectrumPPreconditioning_included
#define EvaluateObservableFermionMatrixSpectrumPPreconditioning_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "EvaluateObservableFermionMatrixSpectrumBase.h"


class EvaluateObservableFermionMatrixSpectrumPPreconditioning : public EvaluateObservableFermionMatrixSpectrumBase {
private:  

public:    
  EvaluateObservableFermionMatrixSpectrumPPreconditioning(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableFermionMatrixSpectrumPPreconditioning();

  int getAnalyzerResultsCount();
};


#include "EvaluateObservableFermionMatrixSpectrumPPreconditioning.C"

#endif
