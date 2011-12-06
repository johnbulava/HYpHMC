#ifndef EvaluateObservableGaussianWeightEstimate_included
#define EvaluateObservableGaussianWeightEstimate_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "AutoCorrelation.h"
#include "MassCorrelationMatrixAnalyzer.h"


class EvaluateObservableGaussianWeightEstimate : public EvaluateObservable {
private:  
  LAPsystemPlot* createPlot1(double rescale);
  LAPsystemPlot* createPlot2(bool gauss);
  double averageWeight;
  double averageNCG;
  double averageNCGSigma;
  double averageMMdagApplic;
  double averageMMdagApplicSigma;
  double averageMMdagApplicNewApproach;
  double averageMMdagApplicNewApproachSigma;
  
public:    
  EvaluateObservableGaussianWeightEstimate(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableGaussianWeightEstimate();

  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
  int getAnalyzerResultsCount();      
  
  double getAverageWeight();
};


#include "EvaluateObservableGaussianWeightEstimate.C"

#endif
