#ifndef EvaluateObservableMultipleTimeScaleIntegration2_included
#define EvaluateObservableMultipleTimeScaleIntegration2_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "AutoCorrelation.h"
#include "MassCorrelationMatrixAnalyzer.h"


class EvaluateObservableMultipleTimeScaleIntegration2 : public EvaluateObservable {
private:  
  bool MultiplePolynomFlag;
  LAPsystemPlot* createPlot1(int startInd, int indCount, const char* tag, const char* des);
  double trajectoryLength;
  
public:    
  EvaluateObservableMultipleTimeScaleIntegration2(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableMultipleTimeScaleIntegration2();

  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
  int getAnalyzerResultsCount();      
};


#endif
