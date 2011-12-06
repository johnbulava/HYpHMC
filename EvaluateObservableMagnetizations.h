#ifndef EvaluateObservableMagnetizations_included
#define EvaluateObservableMagnetizations_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "AutoCorrelation.h"
#include "MassCorrelationMatrixAnalyzer.h"


class EvaluateObservableMagnetizations : public EvaluateObservable {
private:  
  AutoCorrelation* autoCorM;
  AutoCorrelation* autoCorS;
  
  double magM;
  double magS;
  double magMsigma;
  double magSsigma;
  double autoCorMtime;
  double autoCorStime;
  double susM;
  double susS;
  double susMsigma;
  double susSsigma;
  
  
  LAPsystemPlot* createPlot1();
  LAPsystemPlot* createPlot2();
  LAPsystemPlot* createPlot3();
  
public:    
  EvaluateObservableMagnetizations(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableMagnetizations();

  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
  int getAnalyzerResultsCount();      
  
  double getMagnetizationM();
  double getMagnetizationMError();  
  double getMagnetizationS();
  double getMagnetizationSError();    
};


#include "EvaluateObservableMagnetizations.C"

#endif
