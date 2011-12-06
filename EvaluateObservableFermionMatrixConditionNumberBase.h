#ifndef EvaluateObservableFermionMatrixConditionNumberBase_included
#define EvaluateObservableFermionMatrixConditionNumberBase_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"


class EvaluateObservableFermionMatrixConditionNumberBase : public EvaluateObservable {
private:      
  double smallestEW;
  double largestEW;  
  double InverseCondNr;  
    
    
  LAPsystemPlot* createPlot1(bool low, bool logY);
  LAPsystemPlot* createPlot2(bool logY);
  
protected:
  bool drawUpperBound;
  bool drawLowerBound;
    
  
public:    
  EvaluateObservableFermionMatrixConditionNumberBase(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, char* oName, char* nick, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableFermionMatrixConditionNumberBase();

  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
  int getAnalyzerResultsCount();      
};


#include "EvaluateObservableFermionMatrixConditionNumberBase.C"

#endif
