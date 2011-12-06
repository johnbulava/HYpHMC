EvaluateObservableMultipleTimeScaleIntegration::EvaluateObservableMultipleTimeScaleIntegration(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservable(aIOcon, sdr, "MultipleTimeScaleIntegration", "msint", relStart, relEnd) { 
  MultiplePolynomFlag = false;
  if (SDReader->getSubPolynomCount()>0) MultiplePolynomFlag = true;
  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
  trajectoryLength = NaN;
}


EvaluateObservableMultipleTimeScaleIntegration::~EvaluateObservableMultipleTimeScaleIntegration() {
}


int EvaluateObservableMultipleTimeScaleIntegration::getAnalyzerResultsCount() {
  if (MultiplePolynomFlag) {
    return 96;  
  } else {
    return 96;
  }
}


void EvaluateObservableMultipleTimeScaleIntegration::defineObsDependencies() { 
}


bool EvaluateObservableMultipleTimeScaleIntegration::evaluate() {
  trajectoryLength = dataAvail[0][0]*dataAvail[0][2];
  return true;
}


LAPsystemPlot* EvaluateObservableMultipleTimeScaleIntegration::createPlot1(int startInd, int indCount, char* tag, char* des) {
  LAPsystemPlot* plot = LAPsystem->createNewPlot(tag);

  double** plotData = new double*[indCount];
  double maxeps = 0;
  int offset = 4*startInd;
  for (int I2=0; I2<indCount; I2++) {
    plotData[I2] = new double[3];  
    plotData[I2][1] = 0;
    plotData[I2][2] = 0;
    
    for (int I=0; I<dataAvailCount; I++) {
      plotData[I2][0] = dataAvail[I][offset+I2*4+2];
      double dummy = exp(-dataAvail[I][offset+I2*4+1]);      
      if (dummy>1) dummy = 1;
      dummy /= dataAvail[I][offset+I2*4+3];
      plotData[I2][1] += dummy;
      plotData[I2][2] += sqr(dummy);
      if (dataAvail[I][offset+I2*4+2]>maxeps) maxeps = dataAvail[I][offset+I2*4+2];    
    }
    
    plotData[I2][1] /= dataAvailCount;
    plotData[I2][2] = sqrt(plotData[I2][2]/dataAvailCount - sqr(plotData[I2][1])) / sqrt(dataAvailCount);
  }
  
  plot->setPlotData(indCount, 3, plotData);
  plot->setXLabel("Step size $\\\\epsilon$");
  plot->setYLabel("$\\\\langle p_{acc}\\\\rangle$ / number of matrix applications");
  plot->setCaption(des);
  plot->setTitle("");
  plot->setYLogScale(false);
  plot->setSize(0.8, 0.8);
  plot->setXRange(0, 1.1*maxeps);
  plot->setYErrorBars(true);  
  plot->setPointType(3);
  plot->setLineType(2);  
  plot->plotData("1:2:3");
  
  for (int I=0; I<indCount; I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  
  return plot;
}


void EvaluateObservableMultipleTimeScaleIntegration::generateLatexAndPlotsAndXML() {
  if (!MultiplePolynomFlag) {
    LAPsystemPlot* plot1 = createPlot1(0,16,"SimpleLeapFrogIntegration","Simple Leap Frog Integration");
    LAPsystem->addPlot(plot1);
    LAPsystemPlot* plot2 = createPlot1(16,8,"SimpleOmelyanO4FrogIntegration","Simple Omelyan O(4) Frog Integration");
    LAPsystem->addPlot(plot2);
  } else {
    LAPsystemPlot* plot1 = createPlot1(0,8,"MultipleTSouterLeapFrogIntegration","Multiple Time Scales with Leap Frog as outer Integration");
    LAPsystem->addPlot(plot1);
    LAPsystemPlot* plot2 = createPlot1(8,8,"OnlyOuterOmelyanO4Integration","Only outer Omelyan O4 Integration");
    LAPsystem->addPlot(plot2);
    LAPsystemPlot* plot3 = createPlot1(16,8,"OnlyInnerOmelyanO4Integration","Only inner Omelyan O4 Integration");
    LAPsystem->addPlot(plot3);  
  }
  
  startLatexOutputSummaryTable();
  
  addXML_And_LatexOutputSummaryTableLine("TraLen", "Trajectory Length", "$\\tau$", trajectoryLength, 0, NULL, "%1.5f");
 
  endLatexOutputSummaryTable();
}
