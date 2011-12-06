EvaluateObservableGaussianWeightEstimate::EvaluateObservableGaussianWeightEstimate(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservable(aIOcon, sdr, "GaussianWeightEstimate", "gwest", relStart, relEnd) { 
  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
  averageWeight = NaN;
  averageNCG = NaN;
  averageNCGSigma = NaN;
  averageMMdagApplic = NaN;
  averageMMdagApplicSigma = NaN;
  averageMMdagApplicNewApproach = NaN;
  averageMMdagApplicNewApproachSigma = NaN;
}


EvaluateObservableGaussianWeightEstimate::~EvaluateObservableGaussianWeightEstimate() {
}


int EvaluateObservableGaussianWeightEstimate::getAnalyzerResultsCount() {
  return 6;
}


void EvaluateObservableGaussianWeightEstimate::defineObsDependencies() { 
}


bool EvaluateObservableGaussianWeightEstimate::evaluate() {
  averageWeight = 0;
  averageNCG = 0;
  averageNCGSigma = 0;
  averageMMdagApplic = 0;
  averageMMdagApplicSigma = 0;
  averageMMdagApplicNewApproach = 0;
  averageMMdagApplicNewApproachSigma = 0;

  for (int I2=0; I2<dataAvailCount; I2++) {
    averageWeight += dataAvail[I2][0];
    averageNCG += dataAvail[I2][3];
    averageNCGSigma += sqr(dataAvail[I2][3]);
    averageMMdagApplic += dataAvail[I2][4];
    averageMMdagApplicSigma += sqr(dataAvail[I2][4]);
    averageMMdagApplicNewApproach += dataAvail[I2][5];
    averageMMdagApplicNewApproachSigma += sqr(dataAvail[I2][5]);
  }
  averageWeight /= dataAvailCount;
  averageNCG /= dataAvailCount;
  averageNCGSigma = sqrt(averageNCGSigma/dataAvailCount - sqr(averageNCG));
  averageMMdagApplic /= dataAvailCount;
  averageMMdagApplicSigma = sqrt(averageMMdagApplicSigma/dataAvailCount - sqr(averageMMdagApplic));
  averageMMdagApplicNewApproach /= dataAvailCount;
  averageMMdagApplicNewApproachSigma = sqrt(averageMMdagApplicNewApproachSigma/dataAvailCount - sqr(averageMMdagApplicNewApproach));

  return true;
}


LAPsystemPlot* EvaluateObservableGaussianWeightEstimate::createPlot1(double rescale) {
  LAPsystemPlot* plot = NULL;
  if (rescale!=1) {
    plot = LAPsystem->createNewPlot("GaussianWeightEstimateRescaled");
  } else {
    plot = LAPsystem->createNewPlot("GaussianWeightEstimate");
  }

  double** plotData = new double*[dataAvailCount];
  int maxID = 0;
  for (int I=0; I<dataAvailCount; I++) {
    plotData[I] = new double[3];
    plotData[I][0] = dataAvailID[I];
    plotData[I][1] = rescale*dataAvail[I][0];
    plotData[I][2] = rescale*dataAvail[I][1];
    if (dataAvailID[I]>maxID) maxID = dataAvailID[I];
  }
  plot->setPlotData(dataAvailCount, 3, plotData);
  plot->setXLabel("Monte Carlo time");
  plot->setYLabel("Exact weight");
  if (rescale!=1) {
    plot->setCaption("Gaussian weight rescaled");
  } else {
    plot->setCaption("Gaussian weight"); 
  }
  plot->setTitle("");
  plot->setYLogScale(false);
  plot->setSize(0.8, 0.8);
  plot->setXRange(0, maxID);
  plot->setYErrorBars(true);  
  plot->setPointType(3);
  plot->setLineType(2);  
  plot->plotData("1:2:3");
  
  for (int I=0; I<dataAvailCount; I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  
  return plot;
}


LAPsystemPlot* EvaluateObservableGaussianWeightEstimate::createPlot2(bool gauss) {
  LAPsystemPlot* plot = NULL;
  if (gauss) {
    plot = LAPsystem->createNewPlot("AverageMatrixApplicationsGaussApp");
  } else {
    plot = LAPsystem->createNewPlot("AverageMatrixApplicationsExactApp");
  }

  double** plotData = new double*[dataAvailCount];
  int maxID = 0;
  for (int I=0; I<dataAvailCount; I++) {
    plotData[I] = new double[2];
    plotData[I][0] = dataAvailID[I];
    if (gauss) {
      plotData[I][1] = dataAvail[I][4];
    } else {
      plotData[I][1] = dataAvail[I][5];
    }
    if (dataAvailID[I]>maxID) maxID = dataAvailID[I];
  }
  plot->setPlotData(dataAvailCount, 2, plotData);
  plot->setXLabel("Monte Carlo time");
  plot->setYLabel("Average Matrix Applications");
  if (gauss) {
    plot->setCaption("Gaussian approach");
  } else {
    plot->setCaption("Exact approach");  
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


void EvaluateObservableGaussianWeightEstimate::generateLatexAndPlotsAndXML() {
  LAPsystemPlot* plot1 = createPlot1(1.0);
  LAPsystem->addPlot(plot1);
  
  if (averageWeight!=1) {
    LAPsystemPlot* plot2 = createPlot1(1.0/averageWeight);
    LAPsystem->addPlot(plot2);
  }
 
  LAPsystemPlot* plot3 = createPlot2(true);
  LAPsystem->addPlot(plot3);

  LAPsystemPlot* plot4 = createPlot2(false);
  LAPsystem->addPlot(plot4);
  
  startLatexOutputSummaryTable();
  
  addXML_And_LatexOutputSummaryTableLine("AvgWeight", "Average Weight", "$\\langle W\\rangle$", averageWeight, 0, NULL, "%1.5f");
  addXML_And_LatexOutputSummaryTableLine("AvgNCG", "Average number $N_{CG}$", "$\\langle N_{CG}\\rangle$", averageNCG, averageNCGSigma, NULL, "%1.5f");
  addXML_And_LatexOutputSummaryTableLine("AvgMMdag", "Average Number $N_{MM^\\dagger}$", "$\\langle N_{MM^\\dagger}\\rangle$", averageMMdagApplic, averageMMdagApplicSigma, NULL, "%1.5f");
  addXML_And_LatexOutputSummaryTableLine("AvgMMdagNewApp", "Average Number $N^{new}_{MM^\\dagger}$", "$\\langle N_{MM^\\dagger}\\rangle$", averageMMdagApplicNewApproach, averageMMdagApplicNewApproachSigma, NULL, "%1.5f");
 
  endLatexOutputSummaryTable();
}


double EvaluateObservableGaussianWeightEstimate::getAverageWeight() {
  return averageWeight;
}
