#ifndef AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorr_included
#define AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorr_included

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
#include "AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase.h"


class AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorr : public AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase {
private:

protected:

public:
  AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorr(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorr();

};


#include "AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorr.C"

#endif
