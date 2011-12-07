#ifndef EvaluateObservablePsiBarPsiCorrBase_included
#define EvaluateObservablePsiBarPsiCorrBase_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "MassCorrelationMatrixAnalyzer.h"


class EvaluateObservablePsiBarPsiCorrBase : public EvaluateObservable {
private:  
  LAPsystemPlot* createPlot1(bool logY);
  LAPsystemPlot* createPlot2();
  
protected:
  MassCorrelationMatrixAnalyzer* massAnalyzer;

  bool operatorVEVdataAvail; 
  int fitMassCount; 
  
public:    
  EvaluateObservablePsiBarPsiCorrBase(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, char* oName, char* nick, double relStart, double relEnd); 
  ~EvaluateObservablePsiBarPsiCorrBase();

  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
  int getAnalyzerResultsCount();      
};


#endif
