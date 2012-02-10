#include "EvaluateObservableFermionMatrixConditionNumberBase.h"

EvaluateObservableFermionMatrixConditionNumberBase::EvaluateObservableFermionMatrixConditionNumberBase(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, const char* oName, const char* nick, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservable(aIOcon, sdr, oName, nick, relStart, relEnd) { 
  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
  drawUpperBound = false;
  drawLowerBound = false;
}


EvaluateObservableFermionMatrixConditionNumberBase::~EvaluateObservableFermionMatrixConditionNumberBase() {
}


int EvaluateObservableFermionMatrixConditionNumberBase::getAnalyzerResultsCount() {
  return 3;
}


void EvaluateObservableFermionMatrixConditionNumberBase::defineObsDependencies() { 
}


bool EvaluateObservableFermionMatrixConditionNumberBase::evaluate() {
  smallestEW = 0;
  largestEW = 0;  
 
  if (dataAvailCount>0) {
    smallestEW = dataAvail[0][0];
    largestEW = 0;  
  }
 
  for (int I=0; I<dataAvailCount; I++) {
    if (dataAvail[I][0]<smallestEW) smallestEW = dataAvail[I][0];
    if (dataAvail[I][1]>largestEW) largestEW = dataAvail[I][1];
  }

  InverseCondNr = smallestEW / largestEW;
  
  return true;
}


LAPsystemPlot* EvaluateObservableFermionMatrixConditionNumberBase::createPlot1(bool low, bool logY) {
  char* name = new char[1000];
  if (low) {
    snprintf(name,1000,"%sLowestEW", getObsName());
  } else {
    snprintf(name,1000,"%sHighestEW", getObsName());  
  }
  LAPsystemPlot* plot = LAPsystem->createNewPlot(name);

  char* plotCmd = new char[1000];
  double** plotData = new double*[dataAvailCount];
  for (int I=0; I<dataAvailCount; I++) {
    plotData[I] = new double[3];
    plotData[I][0] = dataAvailID[I];
    plotData[I][1] = dataAvail[I][0];
    plotData[I][2] = dataAvail[I][1];    
  }
  plot->setPlotData(dataAvailCount, 3, plotData);
  plot->setXLabel("Monte Carlo time");
  plot->setYLabel("Eigenvalues");
  plot->setCaption("");
  plot->setTitle("");
  plot->setPlotTitle(name);
  plot->setYLogScale(logY);
  plot->setSize(1.0, 1.0);
  plot->setXRange(0, plotData[dataAvailCount-1][0]);
  plot->setYErrorBars(false);  
  plot->setPointSize(0.35);  
  plot->setPointType(5);
//  plot->setLineType(2); 
  if (low) { 
    plot->plotData("1:2");
    if (drawLowerBound) {
      snprintf(plotCmd,1000,"replot %1.15e title 'Lower bound'", SDReader->getPolynomLowerBound_P0());    
      plot->plotDirect(plotCmd);
    }
  } else {
    plot->plotData("1:3");  
    if (drawUpperBound) {
      snprintf(plotCmd,1000,"replot %1.15e title 'Upper bound'", SDReader->getPolynomUpperBound());
      plot->plotDirect(plotCmd);
    }
  }
  
  for (int I=0; I<dataAvailCount; I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  delete[] plotCmd;
  delete[] name;
  
  return plot;
}


LAPsystemPlot* EvaluateObservableFermionMatrixConditionNumberBase::createPlot2(bool logY) {
  char* name = new char[1000];
  snprintf(name,1000,"%sCondNr", getObsName());
  LAPsystemPlot* plot = LAPsystem->createNewPlot(name);

  char* plotCmd = new char[1000];
  double** plotData = new double*[dataAvailCount];
  for (int I=0; I<dataAvailCount; I++) {
    plotData[I] = new double[3];
    plotData[I][0] = dataAvailID[I];
    plotData[I][1] = dataAvail[I][1] / dataAvail[I][0];
    plotData[I][2] = 0;    
  }
  plot->setPlotData(dataAvailCount, 3, plotData);
  plot->setXLabel("Monte Carlo time");
  plot->setYLabel("Condition number $\\\\lambda_{max}/\\\\lambda_{min}$");
  plot->setCaption("");
  plot->setTitle("");
  plot->setPlotTitle(name);
  plot->setYLogScale(logY);
  plot->setSize(1.0, 1.0);
  plot->setXRange(0, plotData[dataAvailCount-1][0]);
  plot->setYErrorBars(false);  
  plot->setPointSize(0.35);  
  plot->setPointType(5);
//  plot->setLineType(2); 
  plot->plotData("1:2");
  
  for (int I=0; I<dataAvailCount; I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  delete[] plotCmd;
  delete[] name;
  
  return plot;
}


void EvaluateObservableFermionMatrixConditionNumberBase::generateLatexAndPlotsAndXML() {
  LAPsystemPlot* plot1 = createPlot1(true, false);  
  LAPsystem->addPlot(plot1);
  
  LAPsystemPlot* plot2 = createPlot1(false, false);  
  LAPsystem->addPlot(plot2);

  LAPsystemPlot* plot3 = createPlot2(false);  
  LAPsystem->addPlot(plot3);
  
  startLatexOutputSummaryTable();

  addXML_And_LatexOutputSummaryTableLine("SmallestEW", "Smallest Eigenwert", "$\\lambda_{min}$",smallestEW, 0, NULL, "%1.1e");
  addXML_And_LatexOutputSummaryTableLine("LargestEW", "Largest Eigenwert", "$\\lambda_{max}$", largestEW, 0, NULL, "%1.2f");
  addXML_And_LatexOutputSummaryTableLine("LargestCondNr", "Largest condition number", "$\\gamma_{max}$", 1.0/InverseCondNr, 0, NULL, "%1.1f");

  endLatexOutputSummaryTable();
}
