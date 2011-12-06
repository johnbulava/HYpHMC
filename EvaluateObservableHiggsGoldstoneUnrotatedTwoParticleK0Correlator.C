EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator::EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservable(aIOcon, sdr, "HiggsGoldstoneUnrotatedTwoParticleK0Correlator", "hgu2pk0corr", relStart, relEnd) { 
  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
  
  NrOfNestedVariables = 4;
  NrOfIndependentVariables = 6;
  int LargestL = getLargestL(SDReader);  
  bool genEWmode = false;
  bool iterEValign = true;
  massAnalyzerCombined = new MassCorrelationMatrixAnalyzer(NrOfIndependentVariables, LargestL, false, genEWmode, iterEValign, true, false, 0, 100000,"CombinedCorr");
  massAnalyzerCombinedWithConstFit = new MassCorrelationMatrixAnalyzer(NrOfIndependentVariables, LargestL, false, genEWmode, iterEValign, true, true, 0, 100000, "CombinedCorrWithConstFit");
}


EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator::~EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator() {
  delete massAnalyzerCombined;
  delete massAnalyzerCombinedWithConstFit;  
}


int EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator::getAnalyzerResultsCount() {
  int LargestL = getLargestL(SDReader);  
  int L0 = SDReader->getL0();
  int L1 = SDReader->getL1();
  int L2 = SDReader->getL2();
  int L3 = SDReader->getL3();

  if ((L0==L1) && (L0==L2) && (L0==L3)) {
    return 4 * 4 * LargestL;
  } else {
    return 4 * LargestL;
  }
}


void EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator::defineObsDependencies() { 
}


bool EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator::evaluate() {
  int LargestL = getLargestL(SDReader);
  int L0 = SDReader->getL0();
  int L1 = SDReader->getL1();
  int L2 = SDReader->getL2();
  int L3 = SDReader->getL3();

  Complex** opDataCombined = new Complex*[LargestL];
  for (int I=0; I<LargestL; I++) {
    opDataCombined[I] = new Complex[NrOfIndependentVariables];
  }
  massAnalyzerCombined->clearData();
  massAnalyzerCombinedWithConstFit->clearData();
  
  int dirLoopMax = 1;
  if ((L0==L1) && (L0==L2) && (L0==L3)) {
    dirLoopMax = 4;
  }
  
  for (int I=0; I<dataAvailCount; I++) {
    for (int dirLoopCount=0; dirLoopCount<dirLoopMax; dirLoopCount++) {
      for (int t=0; t<LargestL; t++) {
        opDataCombined[t][0].x = dataAvail[I][dirLoopCount*4*LargestL + t*NrOfNestedVariables+0];
        opDataCombined[t][0].y = 0;
        opDataCombined[t][1].x = dataAvail[I][dirLoopCount*4*LargestL + t*NrOfNestedVariables+1];
        opDataCombined[t][1].y = 0;
        opDataCombined[t][2].x = dataAvail[I][dirLoopCount*4*LargestL + t*NrOfNestedVariables+2];
        opDataCombined[t][2].y = 0;
        opDataCombined[t][3].x = dataAvail[I][dirLoopCount*4*LargestL + t*NrOfNestedVariables+3];
        opDataCombined[t][3].y = 0;
        opDataCombined[t][4].x = sqr(dataAvail[I][dirLoopCount*4*LargestL + t*NrOfNestedVariables+0]);
        opDataCombined[t][4].y = 0;
        opDataCombined[t][5].x = sqr(dataAvail[I][dirLoopCount*4*LargestL + t*NrOfNestedVariables+1]) + sqr(dataAvail[I][dirLoopCount*4*LargestL + t*NrOfNestedVariables+2]) + sqr(dataAvail[I][dirLoopCount*4*LargestL + t*NrOfNestedVariables+3]);
        opDataCombined[t][5].y = 0;
      }
      massAnalyzerCombined->addOperatorData(dirLoopMax*dataAvailID[I]+dirLoopCount, dataAvailWeightAndSign[I], opDataCombined);
      massAnalyzerCombinedWithConstFit->addOperatorData(dirLoopMax*dataAvailID[I]+dirLoopCount, dataAvailWeightAndSign[I], opDataCombined);
    }
  }
  for (int I=0; I<LargestL; I++) {
    delete[] opDataCombined[I];
  }
 
  delete[] opDataCombined;

  massAnalyzerCombined->calcEigenvaluesAndMassesWithJackKnifeAnalysis(1);
  massAnalyzerCombinedWithConstFit->calcEigenvaluesAndMassesWithJackKnifeAnalysis(1);

  return true;
}


LAPsystemPlot* EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator::createPlot1(bool logY) {
  int LargestL = getLargestL(SDReader);

  char* name = new char[1000];
  snprintf(name,1000,"%sCAnaLogY%d", getObsName(), logY);
  LAPsystemPlot* plot = LAPsystem->createNewPlot(name);

  char* plotCmd = new char[1000];
  double** plotData = new double*[LargestL+1];
  for (int I=0; I<LargestL+1; I++) {
    plotData[I] = new double[1+2*NrOfIndependentVariables];
    plotData[I][0] = I;
    for (int I2=0; I2<NrOfIndependentVariables; I2++) {
      plotData[I][2*I2+1] = massAnalyzerCombined->getMassCorrelationEigenvalue(I, I2);
      plotData[I][2*I2+2] = massAnalyzerCombined->getMassCorrelationEigenvalueError(I, I2);    
    }
  }
  plot->setPlotData(LargestL+1, 1+2*NrOfIndependentVariables, plotData);
  plot->setXLabel("$\\\\Delta t = t_2 - t_1$");
  plot->setYLabel("$\\\\langle h_{t_1}h_{t_2}\\\\rangle$");
  plot->setCaption("Caption");
  plot->setTitle("");
  snprintf(name,1000,"Correlator from Mass-Correlation-Matrix-Analysis of observable %s", getObsName());
  plot->setPlotTitle(name);
  plot->setYLogScale(logY);
  plot->setSize(0.8, 0.8);
  plot->setXRange(0, LargestL);
  plot->setYErrorBars(true);  
  plot->setPointSize(0.35);  
  plot->setPointType(5);
//  plot->setLineType(2);  
  for (int I2=0; I2<NrOfIndependentVariables; I2++) {
    snprintf(plotCmd, 1000, "1:%d:%d",2*I2+2,2*I2+3);
    plot->plotData(plotCmd);
   
    snprintf(plotCmd,1000,"replot %1.15e*cosh(%1.15e*(x-%f)) + %1.15e title '$\\chi^2/dof = %1.2f$'", massAnalyzerCombined->getFittedCoefficient(I2, 0), massAnalyzerCombined->getFittedMass(I2, 0), 0.5*LargestL, massAnalyzerCombined->getFittedVacuumExpectationValue(I2), massAnalyzerCombined->getFittedChiSquare(I2));
    plot->plotDirect(plotCmd);
    snprintf(plotCmd,1000,"replot %1.15e*cosh(%1.15e*(x-%f)) + %1.15e title '$\\chi^2/dof = %1.2f$'", massAnalyzerCombinedWithConstFit->getFittedCoefficient(I2, 0), massAnalyzerCombinedWithConstFit->getFittedMass(I2, 0), 0.5*LargestL, massAnalyzerCombinedWithConstFit->getFittedVacuumExpectationValue(I2), massAnalyzerCombinedWithConstFit->getFittedChiSquare(I2));
    plot->plotDirect(plotCmd);    
  }
  
  for (int I=0; I<LargestL+1; I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  delete[] plotCmd;
  delete[] name;
  
  return plot;
}


LAPsystemPlot* EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator::createPlot2(bool fitConst) {
  int LargestL = getLargestL(SDReader);
  MassCorrelationMatrixAnalyzer* massAna = massAnalyzerCombined;
  if (fitConst) {
    massAna = massAnalyzerCombinedWithConstFit;
  }

  LAPsystemPlot* plot;
  char* name = new char[2000];
  if (fitConst)  snprintf(name,2000,"%sCAnaEffMaWFC", getObsName());
  if (!fitConst) snprintf(name,2000,"%sCAnaEffMa", getObsName());
  plot = LAPsystem->createNewPlot(name);  

  char* plotCmd = new char[1000];
  double** plotData = new double*[LargestL];
  for (int I=0; I<LargestL; I++) {
    plotData[I] = new double[1+2*NrOfIndependentVariables];
    plotData[I][0] = I+0.5;
    for (int I2=0; I2<NrOfIndependentVariables; I2++) {
      plotData[I][2*I2+1] = massAna->getFittedEffectiveMasses(I, I2);
      plotData[I][2*I2+2] = massAna->getFittedEffectiveMassesError(I, I2);    
      if ((isNaN(plotData[I][2*I2+1])) || (isNaN(plotData[I][2*I2+2]))) {
        plotData[I][2*I2+1] = 0;
	plotData[I][2*I2+2] = 0;
      }
    }
  }
  plot->setPlotData(LargestL, 1+2*NrOfIndependentVariables, plotData);
  plot->setXLabel("$\\\\Delta t = t_2 - t_1$");
  plot->setYLabel("$m_H^{eff}(\\\\Delta t)$");
  plot->setCaption("Caption");
  plot->setTitle("");
  if (fitConst)  snprintf(name,2000,"Effective masses $m_H^{eff}(\\\\Delta t)$ from combined analysis of %s with fit to const", getObsName());
  if (!fitConst) snprintf(name,2000,"Effective masses $m_H^{eff}(\\\\Delta t)$ from combined analysis of %s", getObsName());
  plot->setPlotTitle(name);
  plot->setYLogScale(false);
  plot->setSize(1.0, 1.0);
  plot->setXRange(0, LargestL);
  plot->setYErrorBars(true);  
  plot->setPointSize(0.35);  
  plot->setPointType(5);
//  plot->setLineType(2);  

  for (int I2=0; I2<NrOfIndependentVariables; I2++) {
    snprintf(plotCmd, 1000, "1:%d:%d",2*I2+2,2*I2+3);
    plot->plotData(plotCmd);
    snprintf(plotCmd,1000,"replot %1.15e notitle", massAna->getFittedMass(I2, 0));
    plot->plotDirect(plotCmd);
  }
  
  for (int I=0; I<LargestL; I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  delete[] plotCmd;
  delete[] name;
  
  return plot;
}


void EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator::generateLatexAndPlotsAndXML() {
  LAPsystemPlot* plot1 = createPlot1(false);  
  LAPsystem->addPlot(plot1);
  
  LAPsystemPlot* plot2 = createPlot1(true);  
  LAPsystem->addPlot(plot2);  
    
  LAPsystemPlot* plot3 = createPlot2(false);  
  LAPsystem->addPlot(plot3);  

  LAPsystemPlot* plot4 = createPlot2(true);  
  LAPsystem->addPlot(plot4);  

  startLatexOutputSummaryTable();
  char* tag = new char[2000];
  char* ext = new char[2000];
    
  for (int I2=0; I2<NrOfIndependentVariables; I2++) {
    snprintf(ext,2000,"ComA%d",I2);
    snprintf(tag,2000,"LatMass%s",ext);
    addXML_And_LatexOutputSummaryTableLine(tag, "Mass from correlation function in lattice units", "$m_{lat}$", massAnalyzerCombined->getFittedMass(I2, 0), massAnalyzerCombined->getFittedMassError(I2, 0), NULL, "%1.3f");
      
    snprintf(tag,2000,"PhysMass%s",ext);
    addXML_And_LatexOutputSummaryTableLine(tag, "Mass from correlation function in GeV", "$m$", physicalScaleInGEV*massAnalyzerCombined->getFittedMass(I2, 0), sqrt(sqr(physicalScaleInGEV*massAnalyzerCombined->getFittedMassError(I2, 0)) + sqr(physicalScaleErrorInGEV*massAnalyzerCombined->getFittedMass(I2, 0))), "GeV", "%1.1f");
    
    snprintf(tag,2000,"LatMassWFC%s",ext);
    addXML_And_LatexOutputSummaryTableLine(tag, "Mass from correlation function in lattice units with fit to const", "$m_{lat}$", massAnalyzerCombinedWithConstFit->getFittedMass(I2, 0), massAnalyzerCombinedWithConstFit->getFittedMassError(I2, 0), NULL, "%1.3f");
    
    snprintf(tag,2000,"FitConstWFC%s",ext);
    addXML_And_LatexOutputSummaryTableLine(tag, "Const from correlation function fit in lattice units", "$C$", massAnalyzerCombinedWithConstFit->getFittedVacuumExpectationValue(I2), massAnalyzerCombinedWithConstFit->getFittedVacuumExpectationValueError(I2), NULL, "%1.1e");
    
    snprintf(tag,2000,"PhysMassWFC%s",ext);
    addXML_And_LatexOutputSummaryTableLine(tag, "Mass from correlation function in GeV with fit to const", "$m$", physicalScaleInGEV*massAnalyzerCombinedWithConstFit->getFittedMass(I2, 0), sqrt(sqr(physicalScaleInGEV*massAnalyzerCombinedWithConstFit->getFittedMassError(I2, 0)) + sqr(physicalScaleErrorInGEV*massAnalyzerCombinedWithConstFit->getFittedMass(I2, 0))), "GeV", "%1.1f");
  }

  delete[] tag;
  delete[] ext;
  endLatexOutputSummaryTable();
}
