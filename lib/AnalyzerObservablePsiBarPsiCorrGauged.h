#ifndef AnalyzerObservablePsiBarPsiCorrGauged_included
#define AnalyzerObservablePsiBarPsiCorrGauged_included

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


class AnalyzerObservablePsiBarPsiCorrGauged : public AnalyzerObservablePsiBarPsiCorrBase {
private:
  

public:
  AnalyzerObservablePsiBarPsiCorrGauged(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservablePsiBarPsiCorrGauged();
};


#endif
