#ifndef AnalyzerObservableTopBarTopChiralLeftHandedCondensatePointSource_included
#define AnalyzerObservableTopBarTopChiralLeftHandedCondensatePointSource_included

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


class AnalyzerObservableTopBarTopChiralLeftHandedCondensatePointSource : public AnalyzerObservablePsiBarPsiCondensateBase {
private:
  

public:
  AnalyzerObservableTopBarTopChiralLeftHandedCondensatePointSource(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableTopBarTopChiralLeftHandedCondensatePointSource();
};


#endif
