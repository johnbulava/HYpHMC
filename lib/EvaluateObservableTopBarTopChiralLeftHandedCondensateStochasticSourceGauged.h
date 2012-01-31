#ifndef EvaluateObservableTopBarTopChiralLeftHandedCondensateStochasticSourceGauged_included
#define EvaluateObservableTopBarTopChiralLeftHandedCondensateStochasticSourceGauged_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservablePsiBarPsiCondensateBase.h"


class EvaluateObservableTopBarTopChiralLeftHandedCondensateStochasticSourceGauged : public EvaluateObservablePsiBarPsiCondensateBase {
private:  
	
public:    
  EvaluateObservableTopBarTopChiralLeftHandedCondensateStochasticSourceGauged(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableTopBarTopChiralLeftHandedCondensateStochasticSourceGauged();
};


#endif
