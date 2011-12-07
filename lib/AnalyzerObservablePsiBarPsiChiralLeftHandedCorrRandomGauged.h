#ifndef AnalyzerObservablePsiBarPsiChiralLeftHandedCorrRandomGauged_included
#define AnalyzerObservablePsiBarPsiChiralLeftHandedCorrRandomGauged_included

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


class AnalyzerObservablePsiBarPsiChiralLeftHandedCorrRandomGauged : public AnalyzerObservablePsiBarPsiCorrBase {
private:
  

public:
  AnalyzerObservablePsiBarPsiChiralLeftHandedCorrRandomGauged(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservablePsiBarPsiChiralLeftHandedCorrRandomGauged();
};


#endif
