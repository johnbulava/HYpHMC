#ifndef EvaluateObservableMultipleTimeScaleIntegration3_included
#define EvaluateObservableMultipleTimeScaleIntegration3_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "AutoCorrelation.h"
#include "MassCorrelationMatrixAnalyzer.h"


class EvaluateObservableMultipleTimeScaleIntegration3 : public EvaluateObservable {
private:  
  bool MultiplePolynomFlag;
  LAPsystemPlot* createPlot1(int startInd, int indCount, const char* tag, const char* des);
  double trajectoryLength;
  
public:    
  EvaluateObservableMultipleTimeScaleIntegration3(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableMultipleTimeScaleIntegration3();

  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
  int getAnalyzerResultsCount();      
};


#endif
