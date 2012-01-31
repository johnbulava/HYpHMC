#ifndef AnalyzerObservableTopBarTopChiralLeftHandedPropagator_included
#define AnalyzerObservableTopBarTopChiralLeftHandedPropagator_included

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
#include "AnalyzerObservablePsiBarPsiPropagatorBase.h"


class AnalyzerObservableTopBarTopChiralLeftHandedPropagator : public AnalyzerObservablePsiBarPsiPropagatorBase {
private:
  

public:
  AnalyzerObservableTopBarTopChiralLeftHandedPropagator(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableTopBarTopChiralLeftHandedPropagator();
};


#endif
