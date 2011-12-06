EvaluateObservableGoldstonePropagator::EvaluateObservableGoldstonePropagator(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservablePropagatorBase(aIOcon, sdr, obsWeight, obsDetSign, "GoldstonePropagator", "gprop", relStart, relEnd) { 
  considerZeroMomentum = false;
  selfEnergySubLamFac = NaN;
  doArcTanhFit = true;
  gamma = 4.0;
}


EvaluateObservableGoldstonePropagator::~EvaluateObservableGoldstonePropagator() {
}
