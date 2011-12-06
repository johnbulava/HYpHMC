#ifndef EvaluateObservablePsiBarPsiCorrChiralLeftHanded_included
#define EvaluateObservablePsiBarPsiCorrChiralLeftHanded_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservablePsiBarPsiCorrBase.h"


class EvaluateObservablePsiBarPsiCorrChiralLeftHanded : public EvaluateObservablePsiBarPsiCorrBase {
private:  
	
public:    
  EvaluateObservablePsiBarPsiCorrChiralLeftHanded(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservablePsiBarPsiCorrChiralLeftHanded();
};


#include "EvaluateObservablePsiBarPsiCorrChiralLeftHanded.C"

#endif
