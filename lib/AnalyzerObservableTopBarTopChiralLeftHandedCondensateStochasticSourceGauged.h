#ifndef AnalyzerObservableTopBarTopChiralLeftHandedCondensateStochasticSourceGauged_included
#define AnalyzerObservableTopBarTopChiralLeftHandedCondensateStochasticSourceGauged_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "Global.h"
#include "AnalyzerIOControl.h"
#include "FermionMatrixOperations.h"
#include "StateDescriptorReader.h"
#include "AnalyzerObservable.h"
#include "AnalyzerPhiFieldConfiguration.h"
#include "AnalyzerObservablePsiBarPsiCondensateBase.h"


class AnalyzerObservableTopBarTopChiralLeftHandedCondensateStochasticSourceGauged : public AnalyzerObservablePsiBarPsiCondensateBase {
private:
  

public:
  AnalyzerObservableTopBarTopChiralLeftHandedCondensateStochasticSourceGauged(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableTopBarTopChiralLeftHandedCondensateStochasticSourceGauged();
};


#endif
