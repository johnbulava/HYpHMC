#ifndef EvaluateObservableScalarCondensate_included
#define EvaluateObservableScalarCondensate_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "AutoCorrelation.h"
#include "MassCorrelationMatrixAnalyzer.h"


class EvaluateObservableScalarCondensate : public EvaluateObservable {
private:  
  AutoCorrelation* autoCorrS;

  double sCond;
  double sCondError;
  
  
  
public:    
  EvaluateObservableScalarCondensate(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableScalarCondensate();

  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
  int getAnalyzerResultsCount();      
  
};


#endif
