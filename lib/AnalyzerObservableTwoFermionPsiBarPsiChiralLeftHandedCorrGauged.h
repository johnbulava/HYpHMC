#ifndef AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged_included
#define AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged_included

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


class AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged : public AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase {
private:

protected:

public:
  AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged();

};


#endif
