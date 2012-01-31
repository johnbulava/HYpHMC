#ifndef AnalyzerObservablePsiBarPhiPsiChiralLeftHandedCondensate_included
#define AnalyzerObservablePsiBarPhiPsiChiralLeftHandedCondensate_included

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
#include "AnalyzerObservablePsiBarPhiPsiCondensateBase.h"


class AnalyzerObservablePsiBarPhiPsiChiralLeftHandedCondensate : public  AnalyzerObservablePsiBarPhiPsiCondensateBase {
private:
  

public:
  AnalyzerObservablePsiBarPhiPsiChiralLeftHandedCondensate(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservablePsiBarPhiPsiChiralLeftHandedCondensate();
};


#endif
