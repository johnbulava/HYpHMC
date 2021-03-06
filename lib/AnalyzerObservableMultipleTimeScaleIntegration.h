#ifndef AnalyzerObservableMultipleTimeScaleIntegration_included
#define AnalyzerObservableMultipleTimeScaleIntegration_included

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


class AnalyzerObservableMultipleTimeScaleIntegration : public AnalyzerObservable {
private:
  pHMCPropagator* pHMCProp;
  bool MultiplePolynomFlag;
  
  void readFourierAccelerationData();
  double PerformIntegration(int outerIntSteps, int outerIntType, int &Nmmdag, double &epsi, int polyMode);
  int howManyMatrixApplicationsForIntegrator(int* iter, int* Parameter_PolyDegree, int* Parameter_IntegratorType, int polyMode);
  
public:
  AnalyzerObservableMultipleTimeScaleIntegration(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableMultipleTimeScaleIntegration();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#endif
