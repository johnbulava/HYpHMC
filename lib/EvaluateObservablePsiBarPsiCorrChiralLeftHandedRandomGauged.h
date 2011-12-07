#ifndef EvaluateObservablePsiBarPsiCorrChiralLeftHandedRandomGauged_included
#define EvaluateObservablePsiBarPsiCorrChiralLeftHandedRandomGauged_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservablePsiBarPsiCorrBase.h"


class EvaluateObservablePsiBarPsiCorrChiralLeftHandedRandomGauged : public EvaluateObservablePsiBarPsiCorrBase {
private:  
	
public:    
  EvaluateObservablePsiBarPsiCorrChiralLeftHandedRandomGauged(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservablePsiBarPsiCorrChiralLeftHandedRandomGauged();
};


#endif
