#include "EvaluateObservableScalarCondensate.h"

EvaluateObservableScalarCondensate::EvaluateObservableScalarCondensate(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservable(aIOcon, sdr, "ScalarCondensate", "scond", relStart, relEnd) { 
  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
  autoCorrS = new AutoCorrelation(5, 100);  
  sCond = NaN;
  sCondError = NaN;
}


EvaluateObservableScalarCondensate::~EvaluateObservableScalarCondensate() {
  delete autoCorrS;  
}


int EvaluateObservableScalarCondensate::getAnalyzerResultsCount() {
  return 1;
}


void EvaluateObservableScalarCondensate::defineObsDependencies() { 
}


bool EvaluateObservableScalarCondensate::evaluate() {
  double* measureData = new double[dataAvailCount];
  double* measureDataWeights = new double[dataAvailCount];
  
  sCond = 0;
  sCondError = 0;
  //  int L0 = SDReader->getL0();  
  //  int L1 = SDReader->getL1();  
  //  int L2 = SDReader->getL2();  
  //  int L3 = SDReader->getL3();  

  for (int I2=0; I2<dataAvailCount; I2++) {
    measureData[I2] = dataAvail[I2][0];
    measureDataWeights[I2] = dataAvailWeightAndSign[I2];
  } 

  autoCorrS->loadData(1, &dataAvailCount, measureDataWeights, measureData);
  sCond = autoCorrS->getAverage(1);

  ComplexVector derivatives(5);
  derivatives.setZero();
  derivatives.vectorElements[1].x = 1;
  sCondError = autoCorrS->estimateCombinedError(derivatives);

  //autoCorStime = autoCorrS->estimateAutoCorrelationTime();

  delete[] measureData;
  delete[] measureDataWeights;

  return true;
}


void EvaluateObservableScalarCondensate::generateLatexAndPlotsAndXML() {
  startLatexOutputSummaryTable();

  addXML_And_LatexOutputSummaryTableLine("scond", "Scalar Condensate", "$\\langle \\Phi_x \\Phi_{x+\\mu} \\rangle$", sCond, sCondError, NULL, "%1.3f");  

  endLatexOutputSummaryTable();
}


