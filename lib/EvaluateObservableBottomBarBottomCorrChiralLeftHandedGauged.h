#ifndef EvaluateObservableBottomBarBottomCorrChiralLeftHandedGauged_included
#define EvaluateObservableBottomBarBottomCorrChiralLeftHandedGauged_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservablePsiBarPsiCorrBase.h"


class EvaluateObservableBottomBarBottomCorrChiralLeftHandedGauged : public EvaluateObservablePsiBarPsiCorrBase {
private:  
	
public:    
  EvaluateObservableBottomBarBottomCorrChiralLeftHandedGauged(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableBottomBarBottomCorrChiralLeftHandedGauged();
};


#include "EvaluateObservableBottomBarBottomCorrChiralLeftHandedGauged.C"

#endif
