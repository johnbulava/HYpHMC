#ifndef AnalyzerObservableTopBarTopChiralLeftHandedCondensatePointSourceGauged_included
#define AnalyzerObservableTopBarTopChiralLeftHandedCondensatePointSourceGauged_included

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


class AnalyzerObservableTopBarTopChiralLeftHandedCondensatePointSourceGauged : public AnalyzerObservablePsiBarPsiCondensateBase {
private:
  

public:
  AnalyzerObservableTopBarTopChiralLeftHandedCondensatePointSourceGauged(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableTopBarTopChiralLeftHandedCondensatePointSourceGauged();
};


#endif
