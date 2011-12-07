#ifndef EvaluateObservableMultipleTimeScaleIntegration_included
#define EvaluateObservableMultipleTimeScaleIntegration_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "AutoCorrelation.h"
#include "MassCorrelationMatrixAnalyzer.h"


class EvaluateObservableMultipleTimeScaleIntegration : public EvaluateObservable {
private:  
  bool MultiplePolynomFlag;
  LAPsystemPlot* createPlot1(int startInd, int indCount, char* tag, char* des);
  double trajectoryLength;
  
public:    
  EvaluateObservableMultipleTimeScaleIntegration(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableMultipleTimeScaleIntegration();

  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
  int getAnalyzerResultsCount();      
};


#endif
