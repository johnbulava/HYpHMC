#ifndef AnalyzerObservablePsiBarPsiChiralLeftHandedCorr_included
#define AnalyzerObservablePsiBarPsiChiralLeftHandedCorr_included

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


class AnalyzerObservablePsiBarPsiChiralLeftHandedCorr : public AnalyzerObservablePsiBarPsiCorrBase {
private:
  

public:
  AnalyzerObservablePsiBarPsiChiralLeftHandedCorr(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservablePsiBarPsiChiralLeftHandedCorr();
};


#include "AnalyzerObservablePsiBarPsiChiralLeftHandedCorr.C"

#endif
