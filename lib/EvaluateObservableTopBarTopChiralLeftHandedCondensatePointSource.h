#ifndef EvaluateObservableTopBarTopChiralLeftHandedCondensatePointSource_included
#define EvaluateObservableTopBarTopChiralLeftHandedCondensatePointSource_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservablePsiBarPsiCondensateBase.h"


class EvaluateObservableTopBarTopChiralLeftHandedCondensatePointSource : public EvaluateObservablePsiBarPsiCondensateBase {
private:  
	
public:    
  EvaluateObservableTopBarTopChiralLeftHandedCondensatePointSource(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableTopBarTopChiralLeftHandedCondensatePointSource();
};


#endif
