#ifndef EvaluateObservableCorrelatorBase_included
#define EvaluateObservableCorrelatorBase_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "AutoCorrelation.h"
#include "MassCorrelationMatrixAnalyzer.h"


class EvaluateObservableCorrelatorBase : public EvaluateObservable {
private:  
  MassCorrelationMatrixAnalyzer* massAnalyzerCombined;
  MassCorrelationMatrixAnalyzer* massAnalyzerCombinedWithConstFit;
  MassCorrelationMatrixAnalyzer** massAnalyzerSeparate;
  MassCorrelationMatrixAnalyzer** massAnalyzerSeparateWithConstFit;
    
  LAPsystemPlot* createPlot1(bool logY, bool combinedAna);
  LAPsystemPlot* createPlot2(bool fitConst, bool combinedAna);

protected:
  int NrOfNestedVariables;
  int NrOfIndependentVariables;
  bool doMassCorrMatrixAnalysis;
  bool doSeparateAnalysis;
  
  virtual double transformData(double* data, int index);
  
public:    
  EvaluateObservableCorrelatorBase(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, const char* oName, const char* nick, double relStart, double relEnd, int nrOfVar, int nrOfIndVars); 
  ~EvaluateObservableCorrelatorBase();

  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
  int getAnalyzerResultsCount();      
};


#endif
