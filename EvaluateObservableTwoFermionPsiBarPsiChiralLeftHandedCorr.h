#ifndef EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorr_included
#define EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorr_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservablePsiBarPsiCorrBase.h"


class EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorr : public EvaluateObservablePsiBarPsiCorrBase {
private:  
	
public:    
  EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorr(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorr();
};


#include "EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorr.C"

#endif
