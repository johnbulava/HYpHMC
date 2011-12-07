#ifndef EvaluateObservableMultipleTimeScaleIntegrationTest_included
#define EvaluateObservableMultipleTimeScaleIntegrationTest_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "AutoCorrelation.h"
#include "MassCorrelationMatrixAnalyzer.h"


class EvaluateObservableMultipleTimeScaleIntegrationTest : public EvaluateObservable {
private:  
  bool MultiplePolynomFlag;
  LAPsystemPlot* createPlot1(int startInd, int indCount, char* tag, char* des);
  double trajectoryLength;
  
public:    
  EvaluateObservableMultipleTimeScaleIntegrationTest(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableMultipleTimeScaleIntegrationTest();

  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
  int getAnalyzerResultsCount();      
};


#endif
