#ifndef AnalyzerObservablePsiBarPsiChiralRightHandedCorrRandomGauged_included
#define AnalyzerObservablePsiBarPsiChiralRightHandedCorrRandomGauged_included

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


class AnalyzerObservablePsiBarPsiChiralRightHandedCorrRandomGauged : public AnalyzerObservablePsiBarPsiCorrBase {
private:
  

public:
  AnalyzerObservablePsiBarPsiChiralRightHandedCorrRandomGauged(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservablePsiBarPsiChiralRightHandedCorrRandomGauged();
};


#include "AnalyzerObservablePsiBarPsiChiralRightHandedCorrRandomGauged.C"

#endif
