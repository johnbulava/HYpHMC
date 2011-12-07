#include "EvaluateObservableWeight.h"

EvaluateObservableWeight::EvaluateObservableWeight(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservable(aIOcon, sdr, "Weight", "weight", relStart, relEnd) { 
  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
  averageWeight = NaN;
  sigmaWeight = NaN;
  avgLogWeight = NaN;
  sigmaLogWeight = NaN;
}


EvaluateObservableWeight::~EvaluateObservableWeight() {
}


int EvaluateObservableWeight::getAnalyzerResultsCount() {
  return 1;
}


void EvaluateObservableWeight::defineObsDependencies() { 
}


bool EvaluateObservableWeight::evaluate() {
  averageWeight = 0;
  sigmaWeight = 0;
  avgLogWeight = 0;
  sigmaLogWeight = 0;
  for (int I2=0; I2<dataAvailCount; I2++) {
    averageWeight += dataAvail[I2][0];
    sigmaWeight += sqr(dataAvail[I2][0]);
    double dummy = log(dataAvail[I2][0]);
    avgLogWeight += dummy;
    sigmaLogWeight += dummy*dummy;
  }
  averageWeight /= dataAvailCount;
  sigmaWeight = sqrt(sigmaWeight/dataAvailCount - sqr(averageWeight));
  avgLogWeight /= dataAvailCount;
  sigmaLogWeight = sqrt(sigmaLogWeight/dataAvailCount - sqr(avgLogWeight));

  return true;
}


LAPsystemPlot* EvaluateObservableWeight::createPlot1(double rescale) {
  LAPsystemPlot* plot = NULL;
  if (rescale!=1) {
    plot = LAPsystem->createNewPlot("WeightRescaled");
  } else {
    plot = LAPsystem->createNewPlot("Weight");
  }

  double** plotData = new double*[dataAvailCount];
  int maxID = 0;
  for (int I=0; I<dataAvailCount; I++) {
    plotData[I] = new double[2];
    plotData[I][0] = dataAvailID[I];
    plotData[I][1] = rescale*dataAvail[I][0];
    if (dataAvailID[I]>maxID) maxID = dataAvailID[I];    
  }
  plot->setPlotData(dataAvailCount, 2, plotData);
  plot->setXLabel("Monte Carlo time");
  plot->setYLabel("Exact weight");
  if (rescale!=1) {
    plot->setCaption("Exact weight rescaled");
  } else {
    plot->setCaption("Exact weight"); 
  }
  plot->setTitle("");
  plot->setYLogScale(false);
  plot->setSize(0.8, 0.8);
  plot->setXRange(0, maxID);
  plot->setYErrorBars(false);  
  plot->setPointType(3);
  plot->setLineType(2);  
  plot->plotData("1:2");
  
  for (int I=0; I<dataAvailCount; I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  
  return plot;
}


void EvaluateObservableWeight::generateLatexAndPlotsAndXML() {
  LAPsystemPlot* plot1 = createPlot1(1.0);
  LAPsystem->addPlot(plot1);
  
  if (averageWeight!=1) {
    LAPsystemPlot* plot2 = createPlot1(1.0/averageWeight);
    LAPsystem->addPlot(plot2);
  }
  
  startLatexOutputSummaryTable();
  
  addXML_And_LatexOutputSummaryTableLine("AvgWeight", "Average Weight", "$\\langle W\\rangle$", averageWeight, sigmaWeight, NULL, "%1.5f");
  addXML_And_LatexOutputSummaryTableLine("AvgLogWeight", "Average Log of Weight", "$\\langle \\log(W)\\rangle$", avgLogWeight, sigmaLogWeight, NULL, "%1.5f");
 
  endLatexOutputSummaryTable();
}


double EvaluateObservableWeight::getAverageWeight() {
  return averageWeight;
}
