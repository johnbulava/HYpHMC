#ifndef EvaluateObservablePsiBarPsiCorr_included
#define EvaluateObservablePsiBarPsiCorr_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservablePsiBarPsiCorrBase.h"


class EvaluateObservablePsiBarPsiCorr : public EvaluateObservablePsiBarPsiCorrBase {
private:  
	
public:    
  EvaluateObservablePsiBarPsiCorr(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservablePsiBarPsiCorr();
};


#include "EvaluateObservablePsiBarPsiCorr.C"

#endif
