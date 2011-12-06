#ifndef EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorrRandomGauged_included
#define EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorrRandomGauged_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservablePsiBarPsiCorrBase.h"


class EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorrRandomGauged : public EvaluateObservablePsiBarPsiCorrBase {
private:  
	
public:    
  EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorrRandomGauged(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorrRandomGauged();
};


#include "EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorrRandomGauged.C"

#endif
