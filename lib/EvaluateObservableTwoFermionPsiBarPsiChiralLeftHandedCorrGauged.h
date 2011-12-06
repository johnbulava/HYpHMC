#ifndef EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged_included
#define EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservablePsiBarPsiCorrBase.h"


class EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged : public EvaluateObservablePsiBarPsiCorrBase {
private:  
	
public:    
  EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged();
};


#include "EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged.C"

#endif
