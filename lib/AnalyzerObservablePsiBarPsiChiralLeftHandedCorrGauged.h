#ifndef AnalyzerObservablePsiBarPsiChiralLeftHandedCorrGauged_included
#define AnalyzerObservablePsiBarPsiChiralLeftHandedCorrGauged_included

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


class AnalyzerObservablePsiBarPsiChiralLeftHandedCorrGauged : public AnalyzerObservablePsiBarPsiCorrBase {
private:
  

public:
  AnalyzerObservablePsiBarPsiChiralLeftHandedCorrGauged(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservablePsiBarPsiChiralLeftHandedCorrGauged();
};


#endif
