#ifndef AnalyzerObservableTopBarTopChiralLeftHandedCorrGauged_included
#define AnalyzerObservableTopBarTopChiralLeftHandedCorrGauged_included

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
#include "AnalyzerObservablePsiBarPsiCorrBase.h"


class AnalyzerObservableTopBarTopChiralLeftHandedCorrGauged : public AnalyzerObservablePsiBarPsiCorrBase {
private:
  

public:
  AnalyzerObservableTopBarTopChiralLeftHandedCorrGauged(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableTopBarTopChiralLeftHandedCorrGauged();
};


#endif
