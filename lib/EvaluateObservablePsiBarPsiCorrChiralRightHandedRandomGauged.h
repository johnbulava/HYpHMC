#ifndef EvaluateObservablePsiBarPsiCorrChiralRightHandedRandomGauged_included
#define EvaluateObservablePsiBarPsiCorrChiralRightHandedRandomGauged_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservablePsiBarPsiCorrBase.h"


class EvaluateObservablePsiBarPsiCorrChiralRightHandedRandomGauged : public EvaluateObservablePsiBarPsiCorrBase {
private:  
	
public:    
  EvaluateObservablePsiBarPsiCorrChiralRightHandedRandomGauged(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservablePsiBarPsiCorrChiralRightHandedRandomGauged();
};


#include "EvaluateObservablePsiBarPsiCorrChiralRightHandedRandomGauged.C"

#endif
