#ifndef AnalyzerObservableHiggsTopBarTopChiralLeftHandedVertexGauged_included
#define AnalyzerObservableHiggsTopBarTopChiralLeftHandedVertexGauged_included

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
#include "AnalyzerObservableScalarPsiBarPsiVertexBase.h"


class AnalyzerObservableHiggsTopBarTopChiralLeftHandedVertexGauged : public AnalyzerObservableScalarPsiBarPsiVertexBase {
private:
  

public:
  AnalyzerObservableHiggsTopBarTopChiralLeftHandedVertexGauged(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableHiggsTopBarTopChiralLeftHandedVertexGauged();
};


#endif
