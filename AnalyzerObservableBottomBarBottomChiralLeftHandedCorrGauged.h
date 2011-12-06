#ifndef AnalyzerObservableBottomBarBottomChiralLeftHandedCorrGauged_included
#define AnalyzerObservableBottomBarBottomChiralLeftHandedCorrGauged_included

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


class AnalyzerObservableBottomBarBottomChiralLeftHandedCorrGauged : public AnalyzerObservablePsiBarPsiCorrBase {
private:
  

public:
  AnalyzerObservableBottomBarBottomChiralLeftHandedCorrGauged(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableBottomBarBottomChiralLeftHandedCorrGauged();
};


#include "AnalyzerObservableBottomBarBottomChiralLeftHandedCorrGauged.C"

#endif
