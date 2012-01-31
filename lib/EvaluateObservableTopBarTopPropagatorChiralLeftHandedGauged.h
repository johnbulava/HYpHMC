#ifndef EvaluateObservableTopBarTopPropagatorChiralLeftHandedGauged_included
#define EvaluateObservableTopBarTopPropagatorChiralLeftHandedGauged_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservablePsiBarPsiPropagatorBase.h"
#include "NeubergerMatrix.h"


class EvaluateObservableTopBarTopPropagatorChiralLeftHandedGauged : public EvaluateObservablePsiBarPsiPropagatorBase {
private:  
	
public:    
  EvaluateObservableTopBarTopPropagatorChiralLeftHandedGauged(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableTopBarTopPropagatorChiralLeftHandedGauged();
    
};


#endif
