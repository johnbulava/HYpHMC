#ifndef EvaluateObservableTopBarTopChiralLeftHandedCondensateStochasticSource_included
#define EvaluateObservableTopBarTopChiralLeftHandedCondensateStochasticSource_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservablePsiBarPsiCondensateBase.h"


class EvaluateObservableTopBarTopChiralLeftHandedCondensateStochasticSource : public EvaluateObservablePsiBarPsiCondensateBase {
private:  
	
public:    
  EvaluateObservableTopBarTopChiralLeftHandedCondensateStochasticSource(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableTopBarTopChiralLeftHandedCondensateStochasticSource();
};


#endif
