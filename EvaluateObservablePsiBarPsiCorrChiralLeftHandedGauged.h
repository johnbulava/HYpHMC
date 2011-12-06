#ifndef EvaluateObservablePsiBarPsiCorrChiralLeftHandedGauged_included
#define EvaluateObservablePsiBarPsiCorrChiralLeftHandedGauged_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservablePsiBarPsiCorrBase.h"


class EvaluateObservablePsiBarPsiCorrChiralLeftHandedGauged : public EvaluateObservablePsiBarPsiCorrBase {
private:  
	
public:    
  EvaluateObservablePsiBarPsiCorrChiralLeftHandedGauged(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservablePsiBarPsiCorrChiralLeftHandedGauged();
};


#include "EvaluateObservablePsiBarPsiCorrChiralLeftHandedGauged.C"

#endif
