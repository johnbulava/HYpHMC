#ifndef EvaluateObservableFermionMatrixSpectrumNoPreconditioning_included
#define EvaluateObservableFermionMatrixSpectrumNoPreconditioning_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "EvaluateObservableFermionMatrixSpectrumBase.h"


class EvaluateObservableFermionMatrixSpectrumNoPreconditioning : public EvaluateObservableFermionMatrixSpectrumBase {
private:  

public:    
  EvaluateObservableFermionMatrixSpectrumNoPreconditioning(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableFermionMatrixSpectrumNoPreconditioning();


  int getAnalyzerResultsCount();
};


#include "EvaluateObservableFermionMatrixSpectrumNoPreconditioning.C"

#endif
