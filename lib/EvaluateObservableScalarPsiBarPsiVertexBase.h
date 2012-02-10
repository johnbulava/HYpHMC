#ifndef EvaluateObservableScalarPsiBarPsiVertexBase_included
#define EvaluateObservableScalarPsiBarPsiVertexBase_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "MassCorrelationMatrixAnalyzer.h"
#include "EvaluateObservableTopBarTopPropagatorChiralLeftHandedGauged.h"
#include "EvaluateObservableGoldstonePropagator.h"

class EvaluateObservableScalarPsiBarPsiVertexBase : public EvaluateObservable {
private:  
  LAPsystemPlot* createPlot1();
  double* yr;
  double* yrError;
  
  
  
protected:
  MassCorrelationMatrixAnalyzer* massAnalyzer;

  
public:    
  EvaluateObservableScalarPsiBarPsiVertexBase(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, const char* oName, const char* nick, double relStart, double relEnd); 
  ~EvaluateObservableScalarPsiBarPsiVertexBase();

  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
  int getAnalyzerResultsCount();      
};


#endif
