#ifndef EvaluateObservableFermionMatrixSpectrumBase_included
#define EvaluateObservableFermionMatrixSpectrumBase_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"


class EvaluateObservableFermionMatrixSpectrumBase : public EvaluateObservable {
private:      
  LAPsystemPlot* createPlot1();
  
protected:
    
  
public:    
  EvaluateObservableFermionMatrixSpectrumBase(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, const char* oName, const char* nick, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableFermionMatrixSpectrumBase();

  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
};


#endif
