#ifndef EvaluateObservableDetSign_included
#define EvaluateObservableDetSign_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"


class EvaluateObservableDetSign : public EvaluateObservable {
private:  
  int neg_det_count; 
	
public:    
  EvaluateObservableDetSign(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableDetSign();

  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
  int getAnalyzerResultsCount();     
};


#endif
