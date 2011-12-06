EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator::EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservable(aIOcon, sdr, "HiggsGoldstoneUnrotatedTwoParticleK012Correlator", "hgu2pk012corr", relStart, relEnd) { 
  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
  
  NrOfNestedVariables = 52;
  int NrOfIndependentVariablesK01 = 5;
  int NrOfIndependentVariablesK012 = 7;
  int LargestL = getLargestL(SDReader);  
  bool genEWmode = true;
  bool iterEValign = true;
  massAnalyzerCombinedGoldstone = new MassCorrelationMatrixAnalyzer(3, LargestL, false, false, iterEValign, true, false, 0, 100000, "CombinedCorrGold"); 
  massAnalyzerCombinedK01 = new MassCorrelationMatrixAnalyzer(NrOfIndependentVariablesK01, LargestL, false, genEWmode, iterEValign, true, false, 0, 100000,"CombinedCorrK01");
  massAnalyzerCombinedK012 = new MassCorrelationMatrixAnalyzer(NrOfIndependentVariablesK012, LargestL, false, genEWmode, iterEValign, true, true, 0, 100000, "CombinedCorrK012");
  GoldstoneMass = NaN;
  GoldstoneMassError = NaN;
  L0 = SDReader->getL0();  
  L1 = SDReader->getL1();  
  L2 = SDReader->getL2();  
  L3 = SDReader->getL3();    
}


EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator::~EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator() {
  delete massAnalyzerCombinedGoldstone;
  delete massAnalyzerCombinedK01;
  delete massAnalyzerCombinedK012;  
}


int EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator::getAnalyzerResultsCount() {
  int LargestL = getLargestL(SDReader);  
  int L0 = SDReader->getL0();
  int L1 = SDReader->getL1();
  int L2 = SDReader->getL2();
  int L3 = SDReader->getL3();

  if ((L0==L1) && (L0==L2) && (L0==L3)) {
    return 4 * 52 * LargestL;
  } else {
    return 52 * LargestL;
  }
}


void EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator::defineObsDependencies() { 
}


bool EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator::evaluate() {
  int LargestL = getLargestL(SDReader);
  int L0 = SDReader->getL0();
  int L1 = SDReader->getL1();
  int L2 = SDReader->getL2();
  int L3 = SDReader->getL3();

  Complex** opDataCombined = new Complex*[LargestL];
  Complex** opDataCombinedGoldstone = new Complex*[LargestL];
  for (int I=0; I<LargestL; I++) {
    opDataCombined[I] = new Complex[massAnalyzerCombinedK01->getDimension()+massAnalyzerCombinedK012->getDimension()];
    opDataCombinedGoldstone[I] = new Complex[3];    
  }
  massAnalyzerCombinedGoldstone->clearData();  
  massAnalyzerCombinedK01->clearData();
  massAnalyzerCombinedK012->clearData();
  
  int dirLoopMax = 1;
  if ((L0==L1) && (L0==L2) && (L0==L3)) {
    dirLoopMax = 4;
  }

  for (int I=0; I<dataAvailCount; I++) {
    for (int dirLoopCount=0; dirLoopCount<dirLoopMax; dirLoopCount++) {
      for (int t=0; t<LargestL; t++) {
        //Higgs
        opDataCombined[t][0].x = dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+0];
        opDataCombined[t][0].y = 0;
        //Goldstones
        opDataCombinedGoldstone[t][0].x = dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+1];
        opDataCombinedGoldstone[t][0].y = 0;
        opDataCombinedGoldstone[t][1].x = dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+2];
        opDataCombinedGoldstone[t][1].y = 0;
        opDataCombinedGoldstone[t][2].x = dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+3];
        opDataCombinedGoldstone[t][2].y = 0;      
        //Two Higgs rel mom = 0
        opDataCombined[t][1].x = sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+0]);
        opDataCombined[t][1].y = 0;
        //Two Goldstones rel mom = 0
        opDataCombined[t][2].x = sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+1]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+2]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+3]);
        opDataCombined[t][2].y = 0;
      
        //Two Higgs rel mom = 2
        opDataCombined[t][3].x  = sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+0*8+0]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+0*8+1]);
        opDataCombined[t][3].x += sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+1*8+0]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+1*8+1]);
        opDataCombined[t][3].x += sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+2*8+0]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+2*8+1]);      
        opDataCombined[t][3].y = 0;
      
        //Two Goldstones rel mom = 2
        opDataCombined[t][4].x  = sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+0*8+2]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+0*8+3]);
        opDataCombined[t][4].x += sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+0*8+4]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+0*8+5]);
        opDataCombined[t][4].x += sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+0*8+6]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+0*8+7]);
        opDataCombined[t][4].x += sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+1*8+2]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+1*8+3]);
        opDataCombined[t][4].x += sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+1*8+4]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+1*8+5]);
        opDataCombined[t][4].x += sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+1*8+6]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+1*8+7]);
        opDataCombined[t][4].x += sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+2*8+2]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+2*8+3]);
        opDataCombined[t][4].x += sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+2*8+4]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+2*8+5]);
        opDataCombined[t][4].x += sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+2*8+6]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+0*24+2*8+7]);
        opDataCombined[t][4].y = 0;
      
        //Two Higgs rel mom = 4
        opDataCombined[t][5].x  = sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+0*8+0]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+0*8+1]);
        opDataCombined[t][5].x += sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+1*8+0]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+1*8+1]);
        opDataCombined[t][5].x += sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+2*8+0]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+2*8+1]);      
        opDataCombined[t][5].y = 0;
      
        //Two Goldstones rel mom = 4
        opDataCombined[t][6].x  = sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+0*8+2]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+0*8+3]);
        opDataCombined[t][6].x += sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+0*8+4]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+0*8+5]);
        opDataCombined[t][6].x += sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+0*8+6]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+0*8+7]);
        opDataCombined[t][6].x += sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+1*8+2]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+1*8+3]);
        opDataCombined[t][6].x += sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+1*8+4]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+1*8+5]);
        opDataCombined[t][6].x += sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+1*8+6]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+1*8+7]);
        opDataCombined[t][6].x += sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+2*8+2]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+2*8+3]);
        opDataCombined[t][6].x += sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+2*8+4]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+2*8+5]);
        opDataCombined[t][6].x += sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+2*8+6]) + sqr(dataAvail[I][dirLoopCount*52*LargestL + t*NrOfNestedVariables+4+1*24+2*8+7]);
        opDataCombined[t][6].y = 0;
      }
      massAnalyzerCombinedGoldstone->addOperatorData(dirLoopMax*dataAvailID[I]+dirLoopCount, dataAvailWeightAndSign[I], opDataCombinedGoldstone);      
      massAnalyzerCombinedK01->addOperatorData(dirLoopMax*dataAvailID[I]+dirLoopCount, dataAvailWeightAndSign[I], opDataCombined);
      massAnalyzerCombinedK012->addOperatorData(dirLoopMax*dataAvailID[I]+dirLoopCount, dataAvailWeightAndSign[I], opDataCombined);
    }
  }

  for (int I=0; I<LargestL; I++) {
    delete[] opDataCombined[I];
    delete[] opDataCombinedGoldstone[I];    
  }
 
  delete[] opDataCombined;
  delete[] opDataCombinedGoldstone;
  massAnalyzerCombinedGoldstone->calcEigenvaluesAndMassesWithJackKnifeAnalysis(1);
  massAnalyzerCombinedK01->calcEigenvaluesAndMassesWithJackKnifeAnalysis(1); 
//  massAnalyzerCombinedK012->calcEigenvaluesAndMassesWithJackKnifeAnalysis(1);
  
  int* massIndices = new int[massAnalyzerCombinedK01->getDimension()];
  for (int I=0; I<massAnalyzerCombinedK01->getDimension(); I++) {
    massIndices[I] = I;
  }
  bool change = true;
  while (change) {
    change = false;
    for (int I=0; I<massAnalyzerCombinedK01->getDimension()-1; I++) {
      if (massAnalyzerCombinedK01->getFittedMass(massIndices[I], 0) > massAnalyzerCombinedK01->getFittedMass(massIndices[I+1], 0)) {
        int dummy = massIndices[I];
	massIndices[I] = massIndices[I+1];
	massIndices[I+1] = dummy;
	change = true;
      }
    }    
  }
  
  GoldstoneMass  = massAnalyzerCombinedGoldstone->getFittedMass(0, 0);
  GoldstoneMass += massAnalyzerCombinedGoldstone->getFittedMass(1, 0);
  GoldstoneMass += massAnalyzerCombinedGoldstone->getFittedMass(2, 0);
  GoldstoneMass /= 3.0;
  
  GoldstoneMassError  = sqr(massAnalyzerCombinedGoldstone->getFittedMassError(0, 0));
  GoldstoneMassError += sqr(massAnalyzerCombinedGoldstone->getFittedMassError(1, 0));
  GoldstoneMassError += sqr(massAnalyzerCombinedGoldstone->getFittedMassError(2, 0));
  GoldstoneMassError = sqrt(GoldstoneMassError) / 3.0;
  
  delete[] massIndices;
  
  return true;
}


LAPsystemPlot* EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator::createPlot1(bool logY, bool fullK) {
  int LargestL = getLargestL(SDReader);
  MassCorrelationMatrixAnalyzer* massAna = massAnalyzerCombinedK01;
  if (fullK) massAna = massAnalyzerCombinedK012;

  char* name = new char[1000];
  snprintf(name,1000,"%sCAnaLogY%dFullK%d", getObsName(), logY, fullK);
  LAPsystemPlot* plot = LAPsystem->createNewPlot(name);

  char* plotCmd = new char[1000];
  double** plotData = new double*[LargestL+1];
  for (int I=0; I<LargestL+1; I++) {
    plotData[I] = new double[1+2*massAna->getDimension()];
    plotData[I][0] = I;
    for (int I2=0; I2<massAna->getDimension(); I2++) {
      plotData[I][2*I2+1] = massAna->getMassCorrelationEigenvalue(I, I2);
      plotData[I][2*I2+2] = massAna->getMassCorrelationEigenvalueError(I, I2);    
    }
  }
  plot->setPlotData(LargestL+1, 1+2*massAna->getDimension(), plotData);
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
  for (int I2=0; I2<massAna->getDimension(); I2++) {
    snprintf(plotCmd, 1000, "1:%d:%d",2*I2+2,2*I2+3);
    plot->plotData(plotCmd);
   
    snprintf(plotCmd,1000,"replot %1.15e*cosh(%1.15e*(x-%f)) + %1.15e title '$\\chi^2/dof = %1.2f$'", massAna->getFittedCoefficient(I2, 0), massAna->getFittedMass(I2, 0), 0.5*LargestL, massAna->getFittedVacuumExpectationValue(I2), massAna->getFittedChiSquare(I2));
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


LAPsystemPlot* EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator::createPlot2(bool fullK) {
  int LargestL = getLargestL(SDReader);
  MassCorrelationMatrixAnalyzer* massAna = massAnalyzerCombinedK01;
  if (fullK) massAna = massAnalyzerCombinedK012;

  LAPsystemPlot* plot;
  char* name = new char[2000];
  if (fullK)  snprintf(name,2000,"%sCAnaEffMak012", getObsName());
  if (!fullK) snprintf(name,2000,"%sCAnaEffMak01", getObsName());
  plot = LAPsystem->createNewPlot(name);  

  char* plotCmd = new char[1000];
  double** plotData = new double*[LargestL];
  for (int I=0; I<LargestL; I++) {
    plotData[I] = new double[1+2*massAna->getDimension()];
    plotData[I][0] = I+0.5;
    for (int I2=0; I2<massAna->getDimension(); I2++) {
      plotData[I][2*I2+1] = massAna->getFittedEffectiveMasses(I, I2);
      plotData[I][2*I2+2] = massAna->getFittedEffectiveMassesError(I, I2);    
      if ((isNaN(plotData[I][2*I2+1])) || (isNaN(plotData[I][2*I2+2]))) {
        plotData[I][2*I2+1] = 0;
	plotData[I][2*I2+2] = 0;
      }
    }
  }
  plot->setPlotData(LargestL, 1+2*massAna->getDimension(), plotData);
  plot->setXLabel("$\\\\Delta t = t_2 - t_1$");
  plot->setYLabel("$m_H^{eff}(\\\\Delta t)$");
  plot->setCaption("Caption");
  plot->setTitle("");
  if (fullK)  snprintf(name,2000,"Effective masses $m_H^{eff}(\\\\Delta t)$ from combined analysis of %s for k=0,1,2", getObsName());
  if (!fullK) snprintf(name,2000,"Effective masses $m_H^{eff}(\\\\Delta t)$ from combined analysis of %s for k=0,1", getObsName());
  plot->setPlotTitle(name);
  plot->setYLogScale(false);
  plot->setSize(1.0, 1.0);
  plot->setXRange(0, LargestL);
  plot->setYErrorBars(true);  
  plot->setPointSize(0.35);  
  plot->setPointType(5);
//  plot->setLineType(2);  

  for (int I2=0; I2<massAna->getDimension(); I2++) {
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


void EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator::generateLatexAndPlotsAndXML() {
  LAPsystemPlot* plot1 = createPlot1(false, false);  
  LAPsystem->addPlot(plot1);
  
  LAPsystemPlot* plot2 = createPlot1(true, false);  
  LAPsystem->addPlot(plot2);  

  LAPsystemPlot* plot3 = createPlot2(false);  
  LAPsystem->addPlot(plot3);  

  LAPsystem->addDirectText("\\begin{center}\n");
  LAPsystem->addDirectText("Mass correlation matrix analysis for k=0,1\n");
  LAPsystem->addDirectText("\\end{center}\n");

  startLatexOutputSummaryTable();
  char* tag = new char[2000];
  char* ext = new char[2000];
    
  int kCount = 0;
  for (int I2=0; I2<massAnalyzerCombinedK01->getDimension(); I2++) {
    snprintf(ext,2000,"ComAK01_%d",I2);
    snprintf(tag,2000,"LatMass%s",ext);
    addXML_And_LatexOutputSummaryTableLine(tag, "Mass from correlation function in lattice units", "$m_{lat}$", massAnalyzerCombinedK01->getFittedMass(I2, 0),  massAnalyzerCombinedK01->getFittedMassError(I2, 0), NULL, "%1.3f");
      
    snprintf(tag,2000,"PhysMass%s",ext);
    addXML_And_LatexOutputSummaryTableLine(tag, "Mass from correlation function in GeV", "$m$", physicalScaleInGEV*massAnalyzerCombinedK01->getFittedMass(I2, 0), sqrt(sqr(physicalScaleInGEV*massAnalyzerCombinedK01->getFittedMassError(I2, 0)) + sqr(physicalScaleErrorInGEV*massAnalyzerCombinedK01->getFittedMass(I2, 0))), "GeV", "%1.1f");    
    
    double kSqr = sqr(massAnalyzerCombinedK01->getFittedMass(I2, 0)/2) - sqr(GoldstoneMass);
    if ((kSqr>0) && (kSqr<3*sqr(GoldstoneMass))) {
      double k = sqrt(kSqr);
      double kError = sqrt(sqr(massAnalyzerCombinedK01->getFittedMassError(I2, 0)*massAnalyzerCombinedK01->getFittedMass(I2, 0)/2) + sqr(GoldstoneMassError*GoldstoneMass)) / k;
      
      snprintf(tag,2000,"KValComAK01_%d",kCount);
      addXML_And_LatexOutputSummaryTableLine(tag, "Momentum k in lattice units", "$k_{lat}$", k, kError, NULL, "%1.4f");

      double q = k*L0/(2*pi);
      double qError = kError*L0/(2*pi);
      double phiQ = calcLuescherPhiFunction(q, 1E-5);
      double phiQError = abs(qError * calcLuescherPhiDerivative(q, 1E-5));

      snprintf(tag,2000,"PhiQValComAK01_%d",kCount);
      addXML_And_LatexOutputSummaryTableLine(tag, "Phi(q) in lattice units", "$\\phi(q)$", phiQ, phiQError, NULL, "%1.4f");
      
      kCount++;
    }
  }
  for (kCount=kCount; kCount<massAnalyzerCombinedK01->getDimension(); kCount++) {
    snprintf(tag,2000,"KValComAK01_%d",kCount);
    addXMLonly(tag, "Momentum k in lattice units", -1, 1);
    snprintf(tag,2000,"PhiQValComAK01_%d",kCount);
    addXMLonly(tag, "Phi(q) in lattice units", -1, 1);
  }
  addXML_And_LatexOutputSummaryTableLine("AvgGoldstoneLatMass", "Average Goldstone lattice mass from correlation function", "$m_G$", GoldstoneMass, GoldstoneMassError, NULL, "%1.3f");    

  delete[] tag;
  delete[] ext;
  endLatexOutputSummaryTable();

/*  LAPsystemPlot* plot4 = createPlot1(false, true);  
  LAPsystem->addPlot(plot4);
  
  LAPsystemPlot* plot5 = createPlot1(true, true);  
  LAPsystem->addPlot(plot5);  

  LAPsystemPlot* plot6 = createPlot2(true);  
  LAPsystem->addPlot(plot6);  

  LAPsystem->addDirectText("\\begin{center}\n");
  LAPsystem->addDirectText("Mass correlation matrix analysis for k=0,1,2\n");
  LAPsystem->addDirectText("\\end{center}\n");

  startLatexOutputSummaryTable();
  tag = new char[2000];
  ext = new char[2000];
    
  for (int I2=0; I2<massAnalyzerCombinedK012->getDimension(); I2++) {
    snprintf(ext,2000,"ComAK012_%d",I2);
    snprintf(tag,2000,"LatMass%s",ext);
    addXML_And_LatexOutputSummaryTableLine(tag, "Mass from correlation function in lattice units", "$m_{lat}$", massAnalyzerCombinedK012->getFittedMass(I2, 0), massAnalyzerCombinedK012->getFittedMassError(I2, 0), NULL, "%1.3f");
      
    snprintf(tag,2000,"PhysMass%s",ext);
    addXML_And_LatexOutputSummaryTableLine(tag, "Mass from correlation function in GeV", "$m$", physicalScaleInGEV*massAnalyzerCombinedK012->getFittedMass(I2, 0), sqrt(sqr(physicalScaleInGEV*massAnalyzerCombinedK012->getFittedMassError(I2, 0)) + sqr(physicalScaleErrorInGEV*massAnalyzerCombinedK012->getFittedMass(I2, 0))), "GeV", "%1.1f");    
  }

  delete[] tag;
  delete[] ext;
  endLatexOutputSummaryTable();*/
}
