#include "EvaluateObservableHiggsDefinedOnTimeSliceCorrelator.h"

EvaluateObservableHiggsDefinedOnTimeSliceCorrelator::EvaluateObservableHiggsDefinedOnTimeSliceCorrelator(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservableCorrelatorBase(aIOcon, sdr, obsWeight, obsDetSign, "HiggsDefinedOnTimeSliceCorrelator", "hdtcorr", relStart, relEnd, 1,1) { 
  doMassCorrMatrixAnalysis = false;
  doSeparateAnalysis = true;
}


EvaluateObservableHiggsDefinedOnTimeSliceCorrelator::~EvaluateObservableHiggsDefinedOnTimeSliceCorrelator() {
}
