#ifndef EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator_included
#define EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "AutoCorrelation.h"
#include "MassCorrelationMatrixAnalyzer.h"


class EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator : public EvaluateObservable {
private:  
  MassCorrelationMatrixAnalyzer* massAnalyzerCombined;
  MassCorrelationMatrixAnalyzer* massAnalyzerCombinedWithConstFit;

  LAPsystemPlot* createPlot1(bool logY);
  LAPsystemPlot* createPlot2(bool fitConst);
  int NrOfNestedVariables;
  int NrOfIndependentVariables;
  
public:    
  EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator();


  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
  int getAnalyzerResultsCount();      

};


#include "EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator.C"

#endif
