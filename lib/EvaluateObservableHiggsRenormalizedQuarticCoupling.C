#include "EvaluateObservableHiggsRenormalizedQuarticCoupling.h"

EvaluateObservableHiggsRenormalizedQuarticCoupling::EvaluateObservableHiggsRenormalizedQuarticCoupling(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservable(aIOcon, sdr, "HiggsRenormalizedQuarticCoupling", "hrlam", relStart, relEnd) { 
  int L0 = SDReader->getL0();  
  int L1 = SDReader->getL1();  
  int L2 = SDReader->getL2();  
  int L3 = SDReader->getL3();  
  latticeBins = new LatticeMomentumBins(L0,L1,L2,L3);  
  autoCorr = new AutoCorrelation(3, 100);
  avgLamRen = new double[latticeBins->getMomentumSqrSlotCount()];
  sigmaLamRen = new double[latticeBins->getMomentumSqrSlotCount()];
  sigmaLamRenHelper = new double[latticeBins->getMomentumSqrSlotCount()];
  autoCorrelationTimeHiggs4 = new double[latticeBins->getMomentumSqrSlotCount()];
  pSqr = new double[latticeBins->getMomentumSqrSlotCount()]; 
  considerZeroMomentum = true;
  
  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
}


EvaluateObservableHiggsRenormalizedQuarticCoupling::~EvaluateObservableHiggsRenormalizedQuarticCoupling() {
  delete latticeBins;
  delete autoCorr;
  delete[] avgLamRen;  
  delete[] sigmaLamRen;
  delete[] sigmaLamRenHelper;
  delete[] autoCorrelationTimeHiggs4;
  delete[] pSqr;
}


void EvaluateObservableHiggsRenormalizedQuarticCoupling::defineObsDependencies() { 
  addDependOnObsByName("Magnetizations");
  addDependOnObsByName("HiggsPropagator");
}


int EvaluateObservableHiggsRenormalizedQuarticCoupling::getAnalyzerResultsCount() {
  return 1;
}


void EvaluateObservableHiggsRenormalizedQuarticCoupling::calcLamRen(int ignoreStart, int ignoreEnd) {
  ComplexVector derivatives(5);

  double* data = new double[dataAvailCount];
  double* data2 = new double[dataAvailCount];
  double* dataDet = new double[dataAvailCount];  
  int dataStartOffsetForMag = getAnalyzerResultsCount();  
  int dataStartOffsetForProp = dataStartOffsetForMag + (getDependOnObsByName("Magnetizations"))->getAnalyzerResultsCount();

  int count = 0;
  for (int I2=0; I2<dataAvailCount; I2++) {
    if ((I2<ignoreStart) || (I2>ignoreEnd)) {
      data[count] = dataAvail[I2][dataStartOffsetForMag+0];
      dataDet[count] = dataAvailWeightAndSign[I2];	
      count++;
    }
  }
  autoCorr->loadData(1, &(count), dataDet, data);
  double vevMode = autoCorr->getAverage(1);
  double rescale = SimParaSet->reparametrize_HiggsField(SimulationParameterSet_ContinuumNotation);
  rescale *= sqrt(SDReader->getL0()*SDReader->getL1()*SDReader->getL2()*SDReader->getL3());
  vevMode *= rescale;
  double vevModeSqr = vevMode * vevMode;

  count = 0;
  for (int I2=0; I2<dataAvailCount; I2++) {
    if ((I2<ignoreStart) || (I2>ignoreEnd)) {
      data[count] = (rescale*dataAvail[I2][dataStartOffsetForMag+0]) - vevMode;
      dataDet[count] = dataAvailWeightAndSign[I2];	
      count++;
    }
  }
  autoCorr->loadData(1, &(count), dataDet, data);
  double avgHp0Prop = autoCorr->getAverage(2);
  double avgHp0sqrsqr = autoCorr->getAverage(4);
  avgLamRen[0] = -(avgHp0sqrsqr - 1.5*avgHp0Prop*avgHp0Prop);   //Factor 1.5 = 1.(3-1)*1/4 due to different modes averaged
  avgLamRen[0] /= avgHp0Prop*avgHp0Prop*avgHp0Prop*avgHp0Prop;
  pSqr[0] = latticeBins->getLatMomSqrFromSlotNr(0);
  derivatives.setZero();
  derivatives.vectorElements[4].x = 1;  
  if (ignoreStart<=0) autoCorrelationTimeHiggs4[0] = autoCorr->estimateAutoCorrelationTime(derivatives);

  for (int I=1; I<latticeBins->getMomentumSqrSlotCount(); I++) {
    count = 0;
    for (int I2=0; I2<dataAvailCount; I2++) {
      if ((I2<ignoreStart) || (I2>ignoreEnd)) {
        data[count] = dataAvail[I2][dataStartOffsetForProp+I];
        data2[count] = (dataAvail[I2][dataStartOffsetForProp+0]-vevModeSqr) * dataAvail[I2][dataStartOffsetForProp+I];	
        dataDet[count] = dataAvailWeightAndSign[I2];	
	count++;
      }
    }
    autoCorr->loadData(1, &(count), dataDet, data);
    pSqr[I] = latticeBins->getLatMomSqrFromSlotNr(I);
    double avgHProp = autoCorr->getAverage(1);
    autoCorr->loadData(1, &(count), dataDet, data2);
    double avgHsqrsqr = autoCorr->getAverage(1);
    
    avgLamRen[I] = -(avgHsqrsqr - avgHp0Prop*avgHProp);
    avgLamRen[I] /= avgHp0Prop*avgHProp*avgHp0Prop*avgHProp;
    derivatives.setZero();
    derivatives.vectorElements[1].x = 1;      
    if (ignoreStart<=0) autoCorrelationTimeHiggs4[I] = autoCorr->estimateAutoCorrelationTime(derivatives);
  }
  delete[] data;
  delete[] data2;
  delete[] dataDet;
}


bool EvaluateObservableHiggsRenormalizedQuarticCoupling::evaluate() {
  int reps = 0;
  if (dataAvailCount <= 10) {
    reps = dataAvailCount;
  } else if (dataAvailCount <= 100) {
    reps = 10;
  } else {
    reps = ((int) sqrt(dataAvailCount));  
  }

reps = 10;

  double errorRescaleFactor =  sqrt(reps) * (1.0 - 1.0/reps);
  int blockSize = dataAvailCount / reps;
  if (blockSize<1) blockSize = 1;

  for (int I=0; I<latticeBins->getMomentumSqrSlotCount(); I++) {
    sigmaLamRen[I] = 0;
    sigmaLamRenHelper[I] = 0;    
  }

  for (int I=0; I<reps; I++) {
    int igStart = I * blockSize;
    int igEnd = (I+1) * blockSize;
    if (igEnd>=dataAvailCount) igEnd = dataAvailCount-1;
    calcLamRen(igStart, igEnd);

    for (int I2=0; I2<latticeBins->getMomentumSqrSlotCount(); I2++) {    
      sigmaLamRen[I2] += avgLamRen[I2] * avgLamRen[I2];
      sigmaLamRenHelper[I2] += avgLamRen[I2];
    }
  }
  
  calcLamRen(-1, -1);
  
  for (int I2=0; I2<latticeBins->getMomentumSqrSlotCount(); I2++) {
    sigmaLamRen[I2] = errorRescaleFactor * sqrt((sigmaLamRen[I2]/reps)  - sqr(sigmaLamRenHelper[I2]/reps));
  }

  return true;
}


LAPsystemPlot* EvaluateObservableHiggsRenormalizedQuarticCoupling::createPlot1(double maxP) {
  char* name = new char[1000];
  snprintf(name,1000,"%smaxP%1.2f",getObsName(),maxP);
  LAPsystemPlot* plot = LAPsystem->createNewPlot(name);

  int startP = 1;
  if (considerZeroMomentum) startP = 0;
  char* plotCmd = new char[1000];
  double** plotData = new double*[latticeBins->getMomentumSqrSlotCount()-startP];
  for (int I=startP; I<latticeBins->getMomentumSqrSlotCount(); I++) {
    plotData[I-startP] = new double[3];
    plotData[I-startP][0] = pSqr[I];
    plotData[I-startP][1] = avgLamRen[I];
    plotData[I-startP][2] = sigmaLamRen[I];
  }
  plot->setPlotData(latticeBins->getMomentumSqrSlotCount()-startP, 3, plotData);
  plot->setXLabel("$\\\\hat p^2$");
  plot->setYLabel("Renormalized quartic Higgs coupling");
  plot->setCaption("Caption");
  plot->setTitle("");
  plot->setPlotTitle(name);
  plot->setYLogScale(false);
  plot->setSize(1.0, 1.0);
  plot->setXRange(0, maxP);
  plot->setYErrorBars(true);  
  plot->setPointSize(0.3);
  plot->setPointType(5);
  plot->setLineType(0);  
  plot->plotData("1:2:3");
    
  for (int I=0; I<latticeBins->getMomentumSqrSlotCount()-startP; I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  delete[] plotCmd;
  delete[] name;
  
  return plot;
}


LAPsystemPlot* EvaluateObservableHiggsRenormalizedQuarticCoupling::createPlot2() {
  char* name = new char[1000];
  snprintf(name,1000,"%sAutoCorrelationTimesHiggs4",getObsName());
  LAPsystemPlot* plot = LAPsystem->createNewPlot(name);

  int startP = 0;
  if (considerZeroMomentum) startP = 0;
  char* plotCmd = new char[1000];
  double** plotData = new double*[latticeBins->getMomentumSqrSlotCount()-startP];
  for (int I=startP; I<latticeBins->getMomentumSqrSlotCount(); I++) {
    plotData[I-startP] = new double[3];
    plotData[I-startP][0] = pSqr[I];
    plotData[I-startP][1] = autoCorrelationTimeHiggs4[I];
    plotData[I-startP][2] = 0;
  }
  plot->setPlotData(latticeBins->getMomentumSqrSlotCount()-startP, 3, plotData);
  plot->setXLabel("$\\\\hat p^2$");
  plot->setYLabel("Auto-Correlation Times");
  plot->setCaption("");
  plot->setTitle("");
  plot->setPlotTitle(name);
  plot->setYLogScale(false);
  plot->setSize(1.0, 1.0);
//  plot->setXRange(0, maxP);
  plot->setYErrorBars(false);  
  plot->setPointSize(0.3);
  plot->setPointType(5);
  plot->setLineType(0);  
  plot->plotData("1:2");
  
  for (int I=0; I<latticeBins->getMomentumSqrSlotCount()-startP; I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  delete[] plotCmd;
  delete[] name;
  
  return plot;
}


void EvaluateObservableHiggsRenormalizedQuarticCoupling::generateLatexAndPlotsAndXML() {
  LAPsystemPlot* plot1 = createPlot1(16);  
  LAPsystem->addPlot(plot1);

  LAPsystemPlot* plot2 = createPlot1(1);  
  LAPsystem->addPlot(plot2);
  
  LAPsystemPlot* plot3 = createPlot2();  
  LAPsystem->addPlot(plot3);

  startLatexOutputSummaryTable();
  
  addXML_And_LatexOutputSummaryTableLine("LamRenP1", "Renormalized Higgs quartic coupling at p=(0,0,0,1)", "$\\lambda^{h}_{ren}(0,0,0,1)$", avgLamRen[1], sigmaLamRen[1], NULL, "%1.5f");
  
  endLatexOutputSummaryTable();
}
