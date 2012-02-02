#include "EvaluateObservableFermionMatrixSingleMFullSpectrumNoPreconditioning.h"

EvaluateObservableFermionMatrixSingleMFullSpectrumNoPreconditioning::EvaluateObservableFermionMatrixSingleMFullSpectrumNoPreconditioning(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservable(aIOcon, sdr, "FermionMatrixSingleMFullSpectrumNoPreconditioning", "specsmfnp", relStart, relEnd) { 
  MaxPhaseExtension = 0;
  MaxPhaseExtensionFac = 1.25;
  BinCount = 0;
  Bins = NULL;
  L0 = SDReader->getL0();
  L1 = SDReader->getL1();
  L2 = SDReader->getL2();
  L3 = SDReader->getL3();
  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
}


EvaluateObservableFermionMatrixSingleMFullSpectrumNoPreconditioning::~EvaluateObservableFermionMatrixSingleMFullSpectrumNoPreconditioning() {
  delete[] Bins;
}


int EvaluateObservableFermionMatrixSingleMFullSpectrumNoPreconditioning::getAnalyzerResultsCount() {
  return 16*(L0*L1*L2*L3);
}


void EvaluateObservableFermionMatrixSingleMFullSpectrumNoPreconditioning::defineObsDependencies() { 
}


bool EvaluateObservableFermionMatrixSingleMFullSpectrumNoPreconditioning::evaluate() {
  double* detPhases = new double[dataAvailCount];
  double avg = 0;
  double sigma = 0;

  for (int I=0; I<dataAvailCount; I++) {
    Complex normedDet(1, 0);
    double logDet = 0;
    for (int I2=0; I2<getAnalyzerResultsCount()/2; I2++) {
      Complex c(dataAvail[I][2*I2+0], dataAvail[I][2*I2+1]);
      double norm = sqrt(sqr(c.x) + sqr(c.y));
      normedDet = (normedDet * c) / norm;
      logDet += log(norm);
    }
    double norm = sqrt(sqr(normedDet.x) + sqr(normedDet.y));
    double phase = asin(normedDet.y/norm);
    if (normedDet.x<0) phase = pi - phase;
    while (phase<-pi) phase += 2*pi;
    while (phase>pi) phase -= 2*pi;
    if (fabs(phase)>MaxPhaseExtension) MaxPhaseExtension = fabs(phase);
    detPhases[I] = phase;
    avg += phase;
    sigma += sqr(phase);
  }
  avg /= dataAvailCount;
  sigma = sqrt(sigma/dataAvailCount - sqr(avg));
  if (sigma <= 1E-10) {
    BinCount = (int) (1+sqrt(dataAvailCount));
  } else {
    if (dataAvailCount>200) {
      BinCount = (int) (1+0.1*sqrt(dataAvailCount)*MaxPhaseExtension*MaxPhaseExtensionFac/sigma);
    } else {
      BinCount = (int) (1+sqrt(dataAvailCount));    
    }
  }

  delete[] Bins;
  Bins = new double[BinCount];
  for (int I=0; I<BinCount; I++) Bins[I] = 0;
  double delta = (2*MaxPhaseExtension*MaxPhaseExtensionFac) / BinCount;
  double min = -MaxPhaseExtension*MaxPhaseExtensionFac;
  
  for (int I=0; I<dataAvailCount; I++) {
    int binNr = (int) ((detPhases[I]-min) / delta);
    if (binNr>=BinCount) binNr=BinCount-1;
    Bins[binNr] += 1;
  }

  delete[] detPhases;
  return true;
}


LAPsystemPlot* EvaluateObservableFermionMatrixSingleMFullSpectrumNoPreconditioning::createPlot1() {
  char* name = new char[1000];
  snprintf(name,1000,"%s", getObsName());
  LAPsystemPlot* plot = LAPsystem->createNewPlot(name);

  char* plotCmd = new char[1000];
  double** plotData = new double*[BinCount];
  double max = 0;
  for (int I=0; I<BinCount; I++) {
    plotData[I] = new double[2];    
    plotData[I][0] = (-MaxPhaseExtension*MaxPhaseExtensionFac) + ((I+0.5) * (2*MaxPhaseExtension*MaxPhaseExtensionFac)) / BinCount;
    plotData[I][1] = Bins[I]*BinCount / (dataAvailCount*(2*MaxPhaseExtension*MaxPhaseExtensionFac));
    if (plotData[I][1] > max) max = plotData[I][1];
  }
  plot->setPlotData(BinCount, 2, plotData);
  plot->setXLabel("Phase of fermion determinant");
  plot->setYLabel("Relative occurrence density");
  plot->setCaption("");
  plot->setTitle("");
  plot->setPlotTitle(name);
  plot->setYLogScale(true);
  plot->setSize(1.0, 1.0);
//  plot->setYRange(, 1.5*max);
  plot->setYErrorBars(false);  
  plot->setPointSize(0.35);  
  plot->setPointType(5);
//  plot->setLineType(2); 

  snprintf(plotCmd,1000,"set xrange [%1.5e:%1.5e]", (-MaxPhaseExtension*MaxPhaseExtensionFac), (+MaxPhaseExtension*MaxPhaseExtensionFac));    
  plot->plotDirect(plotCmd);
  snprintf(plotCmd,1000,"set yrange [:%1.5e]", 1.5*max);    
  plot->plotDirect(plotCmd);
  
  plot->plotData("1:2 with histeps");
  
  for (int I=0; I<BinCount; I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  delete[] plotCmd;
  delete[] name;
  
  return plot;
}



void EvaluateObservableFermionMatrixSingleMFullSpectrumNoPreconditioning::generateLatexAndPlotsAndXML() {
  LAPsystemPlot* plot1 = createPlot1();  
  LAPsystem->addPlot(plot1);
  
  startLatexOutputSummaryTable();
  endLatexOutputSummaryTable();
}
