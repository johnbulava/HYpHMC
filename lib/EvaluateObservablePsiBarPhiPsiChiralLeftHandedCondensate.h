#ifndef EvaluateObservablePsiBarPhiPsiChiralLeftHandedCondensate_included
#define EvaluateObservablePsiBarPhiPsiChiralLeftHandedCondensate_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservablePsiBarPsiCondensateBase.h"


class EvaluateObservablePsiBarPhiPsiChiralLeftHandedCondensate : public EvaluateObservablePsiBarPsiCondensateBase {
private:  
	
public:    
  EvaluateObservablePsiBarPhiPsiChiralLeftHandedCondensate(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservablePsiBarPhiPsiChiralLeftHandedCondensate();
};


#endif
