#ifndef EvaluateObservableMultipleTimeScaleIntegration4_included
#define EvaluateObservableMultipleTimeScaleIntegration4_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "AutoCorrelation.h"
#include "MassCorrelationMatrixAnalyzer.h"


class EvaluateObservableMultipleTimeScaleIntegration4 : public EvaluateObservable {
private:  
  bool MultiplePolynomFlag;
  LAPsystemPlot* createPlot1(int startInd, int indCount, char* tag, char* des);
  double trajectoryLength;
  
public:    
  EvaluateObservableMultipleTimeScaleIntegration4(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableMultipleTimeScaleIntegration4();

  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
  int getAnalyzerResultsCount();      
};


#include "EvaluateObservableMultipleTimeScaleIntegration4.C"

#endif
