#include "EvaluateObservablePsiBarPsiCorrBase.h"

EvaluateObservablePsiBarPsiCorrBase::EvaluateObservablePsiBarPsiCorrBase(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, char* oName, char* nick, double relStart, double relEnd) : EvaluateObservable(aIOcon, sdr, oName, nick, relStart, relEnd) { 
  operatorVEVdataAvail = false;
  massAnalyzer = NULL;
  fitMassCount = 1;
}


EvaluateObservablePsiBarPsiCorrBase::~EvaluateObservablePsiBarPsiCorrBase() {
  delete massAnalyzer;
}


int EvaluateObservablePsiBarPsiCorrBase::getAnalyzerResultsCount() {
  if (operatorVEVdataAvail) {
    return 2 + 6*getLargestL(SDReader);  
  } else {
    return 2*(1 + getLargestL(SDReader));
  }
}


void EvaluateObservablePsiBarPsiCorrBase::defineObsDependencies() { 
}


bool EvaluateObservablePsiBarPsiCorrBase::evaluate() {
  int LargestL = getLargestL(SDReader);

  Complex*** opCorr = new Complex**[LargestL+1];
  for (int I=0; I<LargestL+1; I++) {
    opCorr[I] = new Complex*[1];
    opCorr[I][0] = new Complex[1];
  }
  Complex** avgOpDag = new Complex*[LargestL];
  Complex** avgOp = new Complex*[LargestL];
  for (int I=0; I<LargestL; I++) {
    avgOpDag[I] = new Complex[1];
    avgOp[I] = new Complex[1];    
  }
  
  massAnalyzer->clearData();
  for (int I=0; I<dataAvailCount; I++) {
    for (int t=0; t<LargestL+1; t++) {
      opCorr[t][0][0].x = dataAvail[I][2*t+0];
      opCorr[t][0][0].y = dataAvail[I][2*t+1];      
    }
    if (operatorVEVdataAvail) {
      for (int t=0; t<LargestL; t++) {
        avgOpDag[t][0].x = dataAvail[I][2+2*LargestL+2*t+0];
        avgOpDag[t][0].y = dataAvail[I][2+2*LargestL+2*t+1];      
        avgOp[t][0].x = dataAvail[I][2+4*LargestL+2*t+0];
        avgOp[t][0].y = dataAvail[I][2+4*LargestL+2*t+1];            
      }    
      massAnalyzer->addOperatorCorrelationData(dataAvailID[I], dataAvailWeightAndSign[I], opCorr, avgOpDag, avgOp);
    } else {
      massAnalyzer->addOperatorCorrelationData(dataAvailID[I], dataAvailWeightAndSign[I], opCorr);    
    }
  }
  for (int I=0; I<LargestL+1; I++) {
    delete[] opCorr[I][0];
    delete[] opCorr[I];
  }
  delete[] opCorr;
  for (int I=0; I<LargestL; I++) {
    delete[] avgOpDag[I];
    delete[] avgOp[I];
  }
  delete[] avgOpDag;
  delete[] avgOp;

  massAnalyzer->calcEigenvaluesAndMassesWithJackKnifeAnalysis(fitMassCount);
  
  return true;
}


LAPsystemPlot* EvaluateObservablePsiBarPsiCorrBase::createPlot1(bool logY) {
  int LargestL = getLargestL(SDReader);

  char* name = new char[1000];
  snprintf(name,1000,"%sLogY%d",getObsName(),logY);
  LAPsystemPlot* plot = LAPsystem->createNewPlot(name);
  delete[] name;

  double** plotData = new double*[LargestL+1];
  for (int I=0; I<LargestL+1; I++) {
    plotData[I] = new double[3];
    plotData[I][0] = I;
    plotData[I][1] = massAnalyzer->getMassCorrelationEigenvalue(I, 0);
    if (logY) plotData[I][1] = abs(massAnalyzer->getMassCorrelationEigenvalue(I, 0));
    plotData[I][2] = massAnalyzer->getMassCorrelationEigenvalueError(I, 0);    
  }
  plot->setPlotData(LargestL+1, 3, plotData);
  plot->setXLabel("$\\\\Delta t = t_2 - t_1$");
  plot->setYLabel("$\\\\langle \\\\psi_{t_1}\\\\bar\\\\psi_{t_2}\\\\rangle$");
  plot->setCaption("Caption");
  plot->setTitle("");
  plot->setPlotTitle(getObsName());
  plot->setYLogScale(logY);
  plot->setSize(0.8, 0.8);
  plot->setXRange(0, LargestL);
  plot->setYErrorBars(true);  
  plot->setPointSize(0.35);
  plot->setPointType(5);
//  plot->setLineType(0);  
  plot->plotData("1:2:3");
  char* fitString = new char[1000];
  char* dummyStr = new char[1000];
  snprintf(fitString,1000,"replot %1.15e",massAnalyzer->getFittedVacuumExpectationValue(0));
  for (int I=0; I<fitMassCount; I++) {
    snprintf(dummyStr,1000,"%s + %1.15e*cosh(%1.15e*(x-%f))",fitString, massAnalyzer->getFittedCoefficient(0, I), massAnalyzer->getFittedMass(0, I), 0.5*LargestL);
    snprintf(fitString,1000,"%s", dummyStr);
  }
  snprintf(dummyStr,1000,"%s title '$\\chi^2/dof = %1.2f$'",fitString,massAnalyzer->getFittedChiSquare(0));
  snprintf(fitString,1000,"%s", dummyStr);
  plot->plotDirect(fitString);  
  delete[] fitString;
  delete[] dummyStr;
  
  for (int I=0; I<LargestL+1; I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  
  return plot;
}


LAPsystemPlot* EvaluateObservablePsiBarPsiCorrBase::createPlot2() {
  int LargestL = getLargestL(SDReader);

  char* name = new char[1000];
  snprintf(name,1000,"%sEffectiveMasses",getObsName());
  LAPsystemPlot* plot = LAPsystem->createNewPlot(name);
  delete[] name;

  char* plotCmd = new char[1000];
  double** plotData = new double*[LargestL/2];
  for (int I=0; I<LargestL/2; I++) {
    plotData[I] = new double[3];
    plotData[I][0] = I+0.5;
    plotData[I][1] = massAnalyzer->getFittedEffectiveMasses(I, 0);
    plotData[I][2] = massAnalyzer->getFittedEffectiveMassesError(I, 0);    
  }
  plot->setPlotData(LargestL/2, 3, plotData);
  plot->setXLabel("$\\\\Delta t = t_2 - t_1$");
  plot->setYLabel("Effective masses $m_{eff}(\\\\Delta t)$");
  plot->setCaption("Caption");
  plot->setTitle("");
  plot->setPlotTitle(getObsName());
  plot->setYLogScale(false);
  plot->setSize(0.8, 0.8);
  plot->setXRange(0, LargestL/2+1);
  plot->setYErrorBars(true);  
  plot->setPointSize(0.35);
  plot->setPointType(5);
//  plot->setLineType(0);  
  plot->plotData("1:2:3");
  snprintf(plotCmd,1000,"replot %1.15e notitle", massAnalyzer->getFittedAsymptoticEffectiveMasses( 0));
  plot->plotDirect(plotCmd);  
  
  for (int I=0; I<LargestL/2; I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  delete[] plotCmd;
  
  return plot;
}


void EvaluateObservablePsiBarPsiCorrBase::generateLatexAndPlotsAndXML() {
  LAPsystemPlot* plot1 = createPlot1(false);  
  LAPsystem->addPlot(plot1);
  
  LAPsystemPlot* plot2 = createPlot1(true);  
  LAPsystem->addPlot(plot2); 

  LAPsystemPlot* plot3 = createPlot2();  
  LAPsystem->addPlot(plot3); 
  
  startLatexOutputSummaryTable();

  double latMass = massAnalyzer->getFittedAsymptoticEffectiveMasses(0);
  double latMassError = massAnalyzer->getFittedAsymptoticEffectiveMassesError(0);  
  addXML_And_LatexOutputSummaryTableLine("LatMass", "Mass from effective masses in lattice units", "$m_{lat}$", latMass, latMassError, NULL, "%1.3f");
  addXML_And_LatexOutputSummaryTableLine("PhysMass", "Mass from effective masses in GeV", "$m$",physicalScaleInGEV*latMass, sqrt(sqr(physicalScaleInGEV*latMassError) + sqr(physicalScaleErrorInGEV*latMass)), "GeV", "%1.1f");
  
  double Yr = (latMass * physicalScaleInGEV) / (Physical_VEV_GeV);
  double YrError = Yr * sqrt(sqr(latMassError / latMass) + sqr(physicalScaleErrorInGEV / physicalScaleInGEV));  
  addXML_And_LatexOutputSummaryTableLine("Yr", "Renormalized Yukawa coupling", "$y_r$",Yr, YrError, NULL, "%1.3f");
  
  endLatexOutputSummaryTable();
}
