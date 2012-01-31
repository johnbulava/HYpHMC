#ifndef EvaluateObservablePsiBarPsiCondensateBase_included
#define EvaluateObservablePsiBarPsiCondensateBase_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "AutoCorrelation.h"

class EvaluateObservablePsiBarPsiCondensateBase : public EvaluateObservable {
private:  
  LAPsystemPlot* createPlot1(bool realPart);
  AutoCorrelation* autoCorR;
  AutoCorrelation* autoCorI;

  
protected:
  Complex condensate;
  Complex condensateSigma;
  int numberOfRepeatedMeasurements;
  char* displayedFormula;

public:    
  EvaluateObservablePsiBarPsiCondensateBase(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, char* oName, char* nick, double relStart, double relEnd); 
  ~EvaluateObservablePsiBarPsiCondensateBase();

  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
  int getAnalyzerResultsCount();      
};


#endif
