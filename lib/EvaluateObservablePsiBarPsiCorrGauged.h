#ifndef EvaluateObservablePsiBarPsiCorrGauged_included
#define EvaluateObservablePsiBarPsiCorrGauged_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservablePsiBarPsiCorrBase.h"


class EvaluateObservablePsiBarPsiCorrGauged : public EvaluateObservablePsiBarPsiCorrBase {
private:  
	
public:    
  EvaluateObservablePsiBarPsiCorrGauged(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservablePsiBarPsiCorrGauged();
};


#endif
