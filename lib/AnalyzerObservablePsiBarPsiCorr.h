#ifndef AnalyzerObservablePsiBarPsiCorr_included
#define AnalyzerObservablePsiBarPsiCorr_included

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


class AnalyzerObservablePsiBarPsiCorr : public AnalyzerObservablePsiBarPsiCorrBase {
private:
  

public:
  AnalyzerObservablePsiBarPsiCorr(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservablePsiBarPsiCorr();
};


#endif
