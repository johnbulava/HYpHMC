#include "EvaluateObservableMagnetizations.h"

EvaluateObservableMagnetizations::EvaluateObservableMagnetizations(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservable(aIOcon, sdr, "Magnetizations", "mags", relStart, relEnd) { 
  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
  autoCorM = new AutoCorrelation(5, 100);
  autoCorS = new AutoCorrelation(5, 100);  
  magM = NaN;
  magS = NaN;
  magMsigma = NaN;
  magSsigma = NaN;
  autoCorMtime = NaN;
  autoCorStime = NaN;
  susM = NaN;
  susS = NaN;
  susMsigma = NaN;
  susSsigma = NaN;  
}


EvaluateObservableMagnetizations::~EvaluateObservableMagnetizations() {
  delete autoCorM;
  delete autoCorS;  
}


int EvaluateObservableMagnetizations::getAnalyzerResultsCount() {
  return 8;
}


void EvaluateObservableMagnetizations::defineObsDependencies() { 
}


bool EvaluateObservableMagnetizations::evaluate() {
  int reps = 0;
  if (dataAvailCount <= 10) {
    reps = dataAvailCount;
  } else if (dataAvailCount <= 100) {
    reps = 10;
  } else {
    reps = ((int) sqrt(dataAvailCount));  
  }
  double errorRescaleFactor =  sqrt(reps) * (1.0 - 1.0/reps);
  int blockSize = dataAvailCount / reps;
  if (blockSize<1) blockSize = 1;

  ComplexVector derivatives(5);
  int* RunLengths = &(dataAvailCount);
  double* measureData = new double[dataAvailCount];
  double* measureData2 = new double[dataAvailCount];
  double* measureDataWeights = new double[dataAvailCount];
  
  susM = 0;
  susS = 0;
  susMsigma = 0;
  susSsigma = 0;  
  int L0 = SDReader->getL0();  
  int L1 = SDReader->getL1();  
  int L2 = SDReader->getL2();  
  int L3 = SDReader->getL3();  
  double volume = L0*L1*L2*L3;
  for (int I=0; I<reps; I++) {
    int igStart = I * blockSize;
    int igEnd = (I+1) * blockSize;
    if (igEnd>=dataAvailCount) igEnd = dataAvailCount-1;

    int count = 0;
    for (int I2=0; I2<dataAvailCount; I2++) {
      if ((I2<igStart) || (I2>igEnd)) {
        measureData[count] = dataAvail[I2][0];
        measureData2[count] = dataAvail[I2][1];
        measureDataWeights[count] = dataAvailWeightAndSign[I2];
        count++;
      }
    }

    autoCorM->loadData(1, &count, measureDataWeights, measureData);
    double dummy = volume * (autoCorM->getAverage(2) - sqr(autoCorM->getAverage(1)));
    susM += dummy;
    susMsigma += sqr(dummy);
    
    autoCorS->loadData(1, &count, measureDataWeights, measureData2);
    dummy = volume * (autoCorS->getAverage(2) - sqr(autoCorS->getAverage(1)));
    susS += dummy;
    susSsigma += sqr(dummy);
  }
  susMsigma = errorRescaleFactor * sqrt(susMsigma/reps - sqr(susM/reps));
  susSsigma = errorRescaleFactor * sqrt(susSsigma/reps - sqr(susS/reps));


  for (int I=0; I<dataAvailCount; I++) {
    measureData[I] = dataAvail[I][0];
  }
  autoCorM->loadData(1, RunLengths, dataAvailWeightAndSign, measureData);
  derivatives.setZero();
  derivatives.vectorElements[1].x = 1;
  magM = autoCorM->getAverage(1);
  magMsigma = autoCorM->estimateCombinedError(derivatives);
  autoCorMtime = autoCorM->estimateAutoCorrelationTime();
  susM = volume * (autoCorM->getAverage(2) - sqr(autoCorM->getAverage(1)));
  
  for (int I=0; I<dataAvailCount; I++) {
    measureData[I] = dataAvail[I][1];
  }
  autoCorS->loadData(1, RunLengths, dataAvailWeightAndSign, measureData);
  derivatives.setZero();
  derivatives.vectorElements[1].x = 1;
  magS = autoCorS->getAverage(1);
  magSsigma = autoCorS->estimateCombinedError(derivatives);
  autoCorStime = autoCorS->estimateAutoCorrelationTime();
  susS = volume * (autoCorS->getAverage(2) - sqr(autoCorS->getAverage(1)));

  delete[] measureData;
  delete[] measureData2;
  delete[] measureDataWeights;

  return true;
}


LAPsystemPlot* EvaluateObservableMagnetizations::createPlot1() {
  LAPsystemPlot* plot = LAPsystem->createNewPlot("MagM");

  double** plotData = new double*[dataAvailCount];
  int smallestID=0;
  int largestID=0;
  if (dataAvailCount>0) {
    smallestID=dataAvailID[0];
    largestID=dataAvailID[0];;  
  }
  for (int I=0; I<dataAvailCount; I++) {
    plotData[I] = new double[2];
    plotData[I][0] = dataAvailID[I];
    plotData[I][1] = dataAvail[I][0];
    if (dataAvailID[I]>largestID) largestID = dataAvailID[I];
    if (dataAvailID[I]<smallestID) smallestID = dataAvailID[I];
  }
  plot->setPlotData(dataAvailCount, 2, plotData);
  plot->setXLabel("Monte Carlo time");
  plot->setYLabel("Magnetization m");
  plot->setCaption("Magnetization m");
  plot->setTitle("");
  plot->setYLogScale(false);
  plot->setSize(1.5, 1.5);
  plot->setXRange(smallestID, largestID);
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


LAPsystemPlot* EvaluateObservableMagnetizations::createPlot2() {
  LAPsystemPlot* plot = LAPsystem->createNewPlot("MagS");

  double** plotData = new double*[dataAvailCount];
  int smallestID=0;
  int largestID=0;
  if (dataAvailCount>0) {
    smallestID=dataAvailID[0];
    largestID=dataAvailID[0];;  
  }
  for (int I=0; I<dataAvailCount; I++) {
    plotData[I] = new double[2];
    plotData[I][0] = dataAvailID[I];
    plotData[I][1] = dataAvail[I][1];
    if (dataAvailID[I]>largestID) largestID = dataAvailID[I];
    if (dataAvailID[I]<smallestID) smallestID = dataAvailID[I];
  }
  plot->setPlotData(dataAvailCount, 2, plotData);
  plot->setXLabel("Monte Carlo time");
  plot->setYLabel("Magnetization s");
  plot->setCaption("Magnetization s");
  plot->setTitle("");
  plot->setYLogScale(false);
  plot->setSize(1.0, 1.0);
  plot->setXRange(smallestID, largestID);
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


LAPsystemPlot* EvaluateObservableMagnetizations::createPlot3() {
  LAPsystemPlot* plot = LAPsystem->createNewPlot("AutoCorrM");

  double* redCfun = autoCorM->getCombinedReducedCFunction();
  double* redCfunError = autoCorM->getCombinedReducedCFunctionErrors();
  double** plotData = new double*[autoCorM->getGammaWMax()];
  for (int I=0; I<autoCorM->getGammaWMax(); I++) {
    plotData[I] = new double[3];
    plotData[I][0] = I;
    plotData[I][1] = redCfun[I];
    plotData[I][2] = redCfunError[I];
  }
  plot->setPlotData(autoCorM->getGammaWMax(), 3, plotData);
  plot->setXLabel("$W$");
  plot->setYLabel("Auto-correlation $C(W)/\\\\Gamma(0)$");
  plot->setCaption("Auto correlation from $\\Gamma$-strategy");
  plot->setTitle("");
  plot->setYLogScale(false);
  plot->setSize(0.8, 0.8);
  plot->setXRange(0, autoCorM->getGammaWMax());
  plot->setYErrorBars(true);  
  plot->setPointType(3);
  plot->setLineType(2);  
  plot->plotData("1:2:3");

  char* plotCmd = new char[1000];
  double fitVal = autoCorM->getTotalN()*sqr(magMsigma)/(autoCorM->getCombinedCFunction())[0];
  snprintf(plotCmd,1000,"replot %1.15e notitle",fitVal);
  plot->plotDirect(plotCmd);
  delete[] plotCmd;
  
  for (int I=0; I<autoCorM->getGammaWMax(); I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  
  return plot;
}



void EvaluateObservableMagnetizations::generateLatexAndPlotsAndXML() {
  LAPsystemPlot* plot1 = createPlot1();
  LAPsystemPlot* plot2 = createPlot2();
  LAPsystemPlot* plot3 = createPlot3();
  
  LAPsystem->addPlot(plot1);
  LAPsystem->addPlot(plot2);
  LAPsystem->addPlot(plot3);
  
  startLatexOutputSummaryTable();

  addXML_And_LatexOutputSummaryTableLine("vev", "VEV of Higgs field (bare Nf-Notation)", "$\\langle \\phi \\rangle$", magM, magMsigma, NULL, "%1.3f");  
  addXML_And_LatexOutputSummaryTableLine("autoCorVEV", "Auto-correlation time of vev", "$\\tau_{int}^{vev}$", autoCorMtime, 0, NULL, "%1.1f");  
  addXML_And_LatexOutputSummaryTableLine("stagvev", "Staggered VEV of Higgs field (bare Nf-Notation)", "$\\langle \\phi \\rangle_{stag}$", magS, magSsigma, NULL, "%1.3f");  
  addXML_And_LatexOutputSummaryTableLine("autoCorStagVEV", "Auto-correlation time of staggered vev", "$\\tau_{int}^{stag}$", autoCorStime, 0, NULL, "%1.1f");  
  addXML_And_LatexOutputSummaryTableLine("SuscepM", "Susceptibility of vev", "$\\chi_m$",susM , susMsigma, NULL, "%1.3f");  
  addXML_And_LatexOutputSummaryTableLine("SuscepS", "Susceptibility of staggered vev", "$\\chi_s$",susS , susSsigma, NULL, "%1.3f");  

  endLatexOutputSummaryTable();
}


double EvaluateObservableMagnetizations::getMagnetizationM() {
  return magM;
}


double EvaluateObservableMagnetizations::getMagnetizationMError() {
  return magMsigma;
}


double EvaluateObservableMagnetizations::getMagnetizationS() {
  return magS;
}


double EvaluateObservableMagnetizations::getMagnetizationSError() {
  return magSsigma;
}
