#include "EvaluateObservableHiggsPropagator.h"

EvaluateObservableHiggsPropagator::EvaluateObservableHiggsPropagator(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservablePropagatorBase(aIOcon, sdr, obsWeight, obsDetSign, "HiggsPropagator", "hprop", relStart, relEnd) { 
  considerZeroMomentum = true;
  
  selfEnergySubLamFac = 12;
  doArcTanhFit = true;  
  gamma = 4.0;
  HiggsPropFourParameterFit = true;
}


EvaluateObservableHiggsPropagator::~EvaluateObservableHiggsPropagator() {
}
