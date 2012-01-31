#ifndef EvaluateObservableTopBarTopChiralLeftHandedCondensatePointSourceGauged_included
#define EvaluateObservableTopBarTopChiralLeftHandedCondensatePointSourceGauged_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservablePsiBarPsiCondensateBase.h"


class EvaluateObservableTopBarTopChiralLeftHandedCondensatePointSourceGauged : public EvaluateObservablePsiBarPsiCondensateBase {
private:  
	
public:    
  EvaluateObservableTopBarTopChiralLeftHandedCondensatePointSourceGauged(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableTopBarTopChiralLeftHandedCondensatePointSourceGauged();
};


#endif
