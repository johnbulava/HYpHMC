#ifndef EvaluateObservableWeight_included
#define EvaluateObservableWeight_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "AutoCorrelation.h"
#include "MassCorrelationMatrixAnalyzer.h"


class EvaluateObservableWeight : public EvaluateObservable {
private:  
  LAPsystemPlot* createPlot1(double rescale);
  double averageWeight;
  double sigmaWeight;
  double avgLogWeight;
  double sigmaLogWeight;
  
public:    
  EvaluateObservableWeight(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableWeight();

  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
  int getAnalyzerResultsCount();      
  
  double getAverageWeight();
};


#endif
