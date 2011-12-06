#ifndef AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrRandomGauged_included
#define AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrRandomGauged_included

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


class AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrRandomGauged : public AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase {
private:

protected:

public:
  AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrRandomGauged(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrRandomGauged();

};


#include "AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrRandomGauged.C"

#endif
