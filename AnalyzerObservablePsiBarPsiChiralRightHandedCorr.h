#ifndef AnalyzerObservablePsiBarPsiChiralRightHandedCorr_included
#define AnalyzerObservablePsiBarPsiChiralRightHandedCorr_included

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


class AnalyzerObservablePsiBarPsiChiralRightHandedCorr : public AnalyzerObservablePsiBarPsiCorrBase {
private:
  

public:
  AnalyzerObservablePsiBarPsiChiralRightHandedCorr(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservablePsiBarPsiChiralRightHandedCorr();
};


#include "AnalyzerObservablePsiBarPsiChiralRightHandedCorr.C"

#endif
