#ifndef AnalyzerObservableMultipleTimeScaleIntegrationTest_included
#define AnalyzerObservableMultipleTimeScaleIntegrationTest_included

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


class AnalyzerObservableMultipleTimeScaleIntegrationTest : public AnalyzerObservable {
private:
  pHMCPropagator* pHMCProp;
  bool MultiplePolynomFlag;
  
  void readFourierAccelerationData();
  double PerformIntegration(int outerIntSteps, int outerIntType, int innerIntSteps, int innerIntType, int &Nmmdag, double &epsi, int polyMode);
  int howManyMatrixApplicationsForIntegrator(int* iter, int* Parameter_PolyDegree, int* Parameter_IntegratorType, int polyMode);
  
public:
  AnalyzerObservableMultipleTimeScaleIntegrationTest(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableMultipleTimeScaleIntegrationTest();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#include "AnalyzerObservableMultipleTimeScaleIntegrationTest.C"

#endif
