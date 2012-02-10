#include "EvaluateObservableFermionMatrixSpectrumBase.h"

EvaluateObservableFermionMatrixSpectrumBase::EvaluateObservableFermionMatrixSpectrumBase(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, const char* oName, const char* nick, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservable(aIOcon, sdr, oName, nick, relStart, relEnd) { 
}


EvaluateObservableFermionMatrixSpectrumBase::~EvaluateObservableFermionMatrixSpectrumBase() {
}


void EvaluateObservableFermionMatrixSpectrumBase::defineObsDependencies() { 
}


bool EvaluateObservableFermionMatrixSpectrumBase::evaluate() {
  return true;
}


LAPsystemPlot* EvaluateObservableFermionMatrixSpectrumBase::createPlot1() {
  char* name = new char[1000];
  snprintf(name,1000,"%s", getObsName());
  LAPsystemPlot* plot = LAPsystem->createNewPlot(name);

  char* plotCmd = new char[1000];
  int size = getAnalyzerResultsCount()/2;
  double** plotData = new double*[dataAvailCount*size];
  for (int I=0; I<dataAvailCount; I++) {
    for (int I2=0; I2<size; I2++) {
      plotData[size*I+I2] = new double[2];    
      plotData[size*I+I2][0] = dataAvail[I][2*I2+0];
      plotData[size*I+I2][1] = dataAvail[I][2*I2+1];
    }
  }
  plot->setPlotData(size*dataAvailCount, 2, plotData);
  plot->setXLabel("Real part");
  plot->setYLabel("Imaginary part");
  plot->setCaption("");
  plot->setTitle("");
  plot->setPlotTitle(name);
  plot->setYLogScale(false);
  plot->setSize(1.0, 1.0);
//  plot->setXRange(0, plotData[dataAvailCount-1][0]);
  plot->setYErrorBars(false);  
  plot->setPointSize(0.35);  
  plot->setPointType(5);
//  plot->setLineType(2); 

  plot->plotData("1:2");
  snprintf(plotCmd,1000,"replot 0 notitle\nset yzeroaxis lt 1");    
  plot->plotDirect(plotCmd);
  
  for (int I=0; I<size*dataAvailCount; I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  delete[] plotCmd;
  delete[] name;
  
  return plot;
}



void EvaluateObservableFermionMatrixSpectrumBase::generateLatexAndPlotsAndXML() {
  LAPsystemPlot* plot1 = createPlot1();  
  LAPsystem->addPlot(plot1);
  
  startLatexOutputSummaryTable();
  endLatexOutputSummaryTable();
}
