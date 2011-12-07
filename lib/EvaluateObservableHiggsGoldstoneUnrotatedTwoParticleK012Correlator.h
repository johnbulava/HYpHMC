#ifndef EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator_included
#define EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "AutoCorrelation.h"
#include "MassCorrelationMatrixAnalyzer.h"


class EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator : public EvaluateObservable {
private:  
  MassCorrelationMatrixAnalyzer* massAnalyzerCombinedGoldstone;
  MassCorrelationMatrixAnalyzer* massAnalyzerCombinedK01;
  MassCorrelationMatrixAnalyzer* massAnalyzerCombinedK012;

  LAPsystemPlot* createPlot1(bool logY, bool fullK);
  LAPsystemPlot* createPlot2(bool fullK);
  int NrOfNestedVariables;
  double GoldstoneMass;
  double GoldstoneMassError;
  int L0;
  int L1;
  int L2;
  int L3;
  
public:    
  EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator();


  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
  int getAnalyzerResultsCount();      

};


#endif
