#ifndef AnalyzerObservableMultipleTimeScaleIntegration2_included
#define AnalyzerObservableMultipleTimeScaleIntegration2_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "Global.h"
#include "AnalyzerIOControl.h"
#include "FermionMatrixOperations.h"
#include "StateDescriptorReader.h"
#include "AnalyzerObservable.h"
#include "AnalyzerPhiFieldConfiguration.h"
#include "pHMCPropagator.h"


class AnalyzerObservableMultipleTimeScaleIntegration2 : public AnalyzerObservable {
private:
  pHMCPropagator* pHMCProp;
  bool MultiplePolynomFlag;
  
  void readFourierAccelerationData();
  double PerformIntegration(int outerIntSteps, int outerIntType, int HiggsIntSteps, int HiggsIntType, int &Nmmdag, double &epsi, int polyMode);
  int howManyMatrixApplicationsForIntegrator(int* iter, int* Parameter_PolyDegree, int* Parameter_IntegratorType, int polyMode);
  
public:
  AnalyzerObservableMultipleTimeScaleIntegration2(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableMultipleTimeScaleIntegration2();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#include "AnalyzerObservableMultipleTimeScaleIntegration2.C"

#endif
