#include "AnalyzerObservableHiggsTopBarTopChiralLeftHandedVertexGauged.h"

AnalyzerObservableHiggsTopBarTopChiralLeftHandedVertexGauged::AnalyzerObservableHiggsTopBarTopChiralLeftHandedVertexGauged(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservableScalarPsiBarPsiVertexBase(fOps, aIOcon, SDreader, "HiggsTopBarTopChiralLeftHandedVertexGauged", "htbtcxlvg") { 
  fixGauge = true;
  projectorSelection = -1;
  PsiPsiBarMatrixInd1Start = 0;
  PsiPsiBarMatrixInd2Start = 0;
  PsiPsiBarMatrixInd1End = 3;
  PsiPsiBarMatrixInd2End = 3;
  tresSumStartIndex = 0;
  tresSumEndIndex = 3;  
}


AnalyzerObservableHiggsTopBarTopChiralLeftHandedVertexGauged::~AnalyzerObservableHiggsTopBarTopChiralLeftHandedVertexGauged() {
}
