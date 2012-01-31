#ifndef EvaluateObservableHiggsTopBarTopChiralLeftHandedVertexGauged_included
#define EvaluateObservableHiggsTopBarTopChiralLeftHandedVertexGauged_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservableScalarPsiBarPsiVertexBase.h"


class EvaluateObservableHiggsTopBarTopChiralLeftHandedVertexGauged : public EvaluateObservableScalarPsiBarPsiVertexBase {
private:  
	
public:    
  EvaluateObservableHiggsTopBarTopChiralLeftHandedVertexGauged(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableHiggsTopBarTopChiralLeftHandedVertexGauged();
};


#endif
