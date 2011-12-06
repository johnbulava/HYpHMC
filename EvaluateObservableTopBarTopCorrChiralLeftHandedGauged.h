#ifndef EvaluateObservableTopBarTopCorrChiralLeftHandedGauged_included
#define EvaluateObservableTopBarTopCorrChiralLeftHandedGauged_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservablePsiBarPsiCorrBase.h"


class EvaluateObservableTopBarTopCorrChiralLeftHandedGauged : public EvaluateObservablePsiBarPsiCorrBase {
private:  
	
public:    
  EvaluateObservableTopBarTopCorrChiralLeftHandedGauged(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableTopBarTopCorrChiralLeftHandedGauged();
};


#include "EvaluateObservableTopBarTopCorrChiralLeftHandedGauged.C"

#endif
