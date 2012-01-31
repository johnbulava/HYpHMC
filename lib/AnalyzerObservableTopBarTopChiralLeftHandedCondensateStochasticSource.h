#ifndef AnalyzerObservableTopBarTopChiralLeftHandedCondensateStochasticSource_included
#define AnalyzerObservableTopBarTopChiralLeftHandedCondensateStochasticSource_included

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


class AnalyzerObservableTopBarTopChiralLeftHandedCondensateStochasticSource : public AnalyzerObservablePsiBarPsiCondensateBase {
private:
  

public:
  AnalyzerObservableTopBarTopChiralLeftHandedCondensateStochasticSource(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableTopBarTopChiralLeftHandedCondensateStochasticSource();
};


#endif
