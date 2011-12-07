#ifndef EvaluateObservablePsiBarPsiCorrChiralRightHanded_included
#define EvaluateObservablePsiBarPsiCorrChiralRightHanded_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservablePsiBarPsiCorrBase.h"


class EvaluateObservablePsiBarPsiCorrChiralRightHanded : public EvaluateObservablePsiBarPsiCorrBase {
private:  
	
public:    
  EvaluateObservablePsiBarPsiCorrChiralRightHanded(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservablePsiBarPsiCorrChiralRightHanded();
};


#endif
