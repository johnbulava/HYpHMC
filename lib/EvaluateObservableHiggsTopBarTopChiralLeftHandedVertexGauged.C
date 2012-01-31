#include "EvaluateObservableHiggsTopBarTopChiralLeftHandedVertexGauged.h"

EvaluateObservableHiggsTopBarTopChiralLeftHandedVertexGauged::EvaluateObservableHiggsTopBarTopChiralLeftHandedVertexGauged(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservableScalarPsiBarPsiVertexBase(aIOcon, sdr, obsWeight, obsDetSign, "HiggsTopBarTopChiralLeftHandedVertexGauged", "htbtcxlvg", relStart, relEnd) { 
  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
}


EvaluateObservableHiggsTopBarTopChiralLeftHandedVertexGauged::~EvaluateObservableHiggsTopBarTopChiralLeftHandedVertexGauged() {
}
