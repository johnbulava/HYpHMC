#ifndef EvaluateObservableFermionMatrixSingleMFullSpectrumNoPreconditioning_included
#define EvaluateObservableFermionMatrixSingleMFullSpectrumNoPreconditioning_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"


class EvaluateObservableFermionMatrixSingleMFullSpectrumNoPreconditioning : public EvaluateObservable {
private:      
  LAPsystemPlot* createPlot1();
  double MaxPhaseExtension;
  double MaxPhaseExtensionFac;
  int BinCount;
  double* Bins;
  int L0;
  int L1;
  int L2;
  int L3;
  
  
protected:
    
  
public:    
  EvaluateObservableFermionMatrixSingleMFullSpectrumNoPreconditioning(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd); 
  ~EvaluateObservableFermionMatrixSingleMFullSpectrumNoPreconditioning();

  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
  int getAnalyzerResultsCount();      
};


#endif
