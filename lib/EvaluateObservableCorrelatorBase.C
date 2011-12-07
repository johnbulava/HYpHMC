#include "EvaluateObservableCorrelatorBase.h"

EvaluateObservableCorrelatorBase::EvaluateObservableCorrelatorBase(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, char* oName, char* nick, double relStart, double relEnd, int nrOfVar, int nrOfIndVars) : EvaluateObservable(aIOcon, sdr, oName, nick, relStart, relEnd) { 
  NrOfNestedVariables = nrOfVar;
  NrOfIndependentVariables = nrOfIndVars;
  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
  doMassCorrMatrixAnalysis = false;
  doSeparateAnalysis = true;
  
  int LargestL = getLargestL(SDReader);  
  massAnalyzerCombined = new MassCorrelationMatrixAnalyzer(NrOfIndependentVariables, LargestL, false, false, true, true, false, 0, 100000,"CombinedCorr");
  massAnalyzerCombinedWithConstFit = new MassCorrelationMatrixAnalyzer(NrOfIndependentVariables, LargestL, false, false, true, true, true, 0, 100000, "CombinedCorrWithConstFit");
  massAnalyzerSeparate = new MassCorrelationMatrixAnalyzer*[NrOfIndependentVariables];
  massAnalyzerSeparateWithConstFit = new MassCorrelationMatrixAnalyzer*[NrOfIndependentVariables];
  for (int I=0; I<NrOfIndependentVariables; I++) {
    massAnalyzerSeparate[I] = new MassCorrelationMatrixAnalyzer(NrOfIndependentVariables, LargestL, false, false, true, true, false, 0, 100000, "SeparateCorr");  
    massAnalyzerSeparateWithConstFit[I] = new MassCorrelationMatrixAnalyzer(NrOfIndependentVariables, LargestL, false, false, true, true, true, 0, 100000, "SeparateCorrWithConstFit");  
  }
}


EvaluateObservableCorrelatorBase::~EvaluateObservableCorrelatorBase() {
  delete massAnalyzerCombined;
  delete massAnalyzerCombinedWithConstFit;  
  for (int I=0; I<NrOfIndependentVariables; I++) {
    delete massAnalyzerSeparate[I];
    delete massAnalyzerSeparateWithConstFit[I];
  }
  delete[] massAnalyzerSeparate;
  delete[] massAnalyzerSeparateWithConstFit;
}


int EvaluateObservableCorrelatorBase::getAnalyzerResultsCount() {
  int LargestL = getLargestL(SDReader);  
  int L0 = SDReader->getL0();
  int L1 = SDReader->getL1();
  int L2 = SDReader->getL2();
  int L3 = SDReader->getL3();

  if ((L0==L1) && (L0==L2) && (L0==L3)) {
    return 4*NrOfNestedVariables * LargestL;
  } else {
    return NrOfNestedVariables * LargestL;
  }
}


double EvaluateObservableCorrelatorBase::transformData(double* data, int index) {
  return data[index];
}


void EvaluateObservableCorrelatorBase::defineObsDependencies() { 
}


bool EvaluateObservableCorrelatorBase::evaluate() {
  int LargestL = getLargestL(SDReader);
  int L0 = SDReader->getL0();
  int L1 = SDReader->getL1();
  int L2 = SDReader->getL2();
  int L3 = SDReader->getL3();

  Complex** opDataSeparate = new Complex*[LargestL];
  Complex** opDataCombined = new Complex*[LargestL];
  for (int I=0; I<LargestL; I++) {
    opDataSeparate[I] = new Complex[1];
    opDataCombined[I] = new Complex[NrOfIndependentVariables];
  }
  massAnalyzerCombined->clearData();
  massAnalyzerCombinedWithConstFit->clearData();
  for (int I=0; I<NrOfIndependentVariables; I++) {
    massAnalyzerSeparate[I]->clearData();
    massAnalyzerSeparateWithConstFit[I]->clearData();
  }
  int dirLoopMax = 1;
  if ((L0==L1) && (L0==L2) && (L0==L3)) {
    dirLoopMax = 4;
  }
  
  for (int I=0; I<dataAvailCount; I++) {
    for (int dirLoopCount=0; dirLoopCount<dirLoopMax; dirLoopCount++) {
      for (int I2=0; I2<NrOfIndependentVariables; I2++) {
        for (int t=0; t<LargestL; t++) {
          opDataSeparate[t][0].x = transformData(&(dataAvail[I][dirLoopCount*NrOfNestedVariables*LargestL + t*NrOfNestedVariables]), I2);
          opDataSeparate[t][0].y = 0;
	  opDataCombined[t][I2].x = transformData(&(dataAvail[I][dirLoopCount*NrOfNestedVariables*LargestL + t*NrOfNestedVariables]), I2);
  	  opDataCombined[t][I2].y = 0;
        }
        massAnalyzerSeparate[I2]->addOperatorData(dirLoopMax*dataAvailID[I]+dirLoopCount, dataAvailWeightAndSign[I], opDataSeparate);
        massAnalyzerSeparateWithConstFit[I2]->addOperatorData(dirLoopMax*dataAvailID[I]+dirLoopCount, dataAvailWeightAndSign[I], opDataSeparate);
      }
      massAnalyzerCombined->addOperatorData(dirLoopMax*dataAvailID[I]+dirLoopCount, dataAvailWeightAndSign[I], opDataCombined);
      massAnalyzerCombinedWithConstFit->addOperatorData(dirLoopMax*dataAvailID[I]+dirLoopCount, dataAvailWeightAndSign[I], opDataCombined);
    }
  }
  for (int I=0; I<LargestL; I++) {
    delete[] opDataSeparate[I];
    delete[] opDataCombined[I];
  }
 
  delete[] opDataSeparate;
  delete[] opDataCombined;

  if (doMassCorrMatrixAnalysis) {
    massAnalyzerCombined->calcEigenvaluesAndMassesWithJackKnifeAnalysis(1);
    massAnalyzerCombinedWithConstFit->calcEigenvaluesAndMassesWithJackKnifeAnalysis(1);
  }
  if (doSeparateAnalysis) {
    for (int I=0; I<NrOfIndependentVariables; I++) {
      massAnalyzerSeparate[I]->calcEigenvaluesAndMassesWithJackKnifeAnalysis(1);
      massAnalyzerSeparateWithConstFit[I]->calcEigenvaluesAndMassesWithJackKnifeAnalysis(1);
    }
  }

  return true;
}


LAPsystemPlot* EvaluateObservableCorrelatorBase::createPlot1(bool logY, bool combinedAna) {
  int LargestL = getLargestL(SDReader);

  char* name = new char[1000];
  if (combinedAna) {
    snprintf(name,1000,"%sCAnaLogY%d", getObsName(), logY);
  } else {
    snprintf(name,1000,"%sSAnaLogY%d", getObsName(), logY);
  }
  LAPsystemPlot* plot = LAPsystem->createNewPlot(name);

  char* plotCmd = new char[1000];
  double** plotData = new double*[LargestL+1];
  for (int I=0; I<LargestL+1; I++) {
    plotData[I] = new double[1+2*NrOfIndependentVariables];
    plotData[I][0] = I;
    for (int I2=0; I2<NrOfIndependentVariables; I2++) {
      if (combinedAna) {
        plotData[I][2*I2+1] = massAnalyzerCombined->getMassCorrelationEigenvalue(I, I2);
        plotData[I][2*I2+2] = massAnalyzerCombined->getMassCorrelationEigenvalueError(I, I2);    
      } else {
        plotData[I][2*I2+1] = massAnalyzerSeparate[I2]->getMassCorrelationEigenvalue(I, 0);
        plotData[I][2*I2+2] = massAnalyzerSeparate[I2]->getMassCorrelationEigenvalueError(I, 0);    
      }
    }
  }
  plot->setPlotData(LargestL+1, 1+2*NrOfIndependentVariables, plotData);
  plot->setXLabel("$\\\\Delta t = t_2 - t_1$");
  plot->setYLabel("$\\\\langle h_{t_1}h_{t_2}\\\\rangle$");
  plot->setCaption("Caption");
  plot->setTitle("");
  if (combinedAna) {
    snprintf(name,1000,"Correlator from Mass-Correlation-Matrix-Analysis of observable %s", getObsName());
  } else {
    snprintf(name,1000,"Correlator from separate analysis of observable %s", getObsName());
  }
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
   
    if (combinedAna) {
      snprintf(plotCmd,1000,"replot %1.15e*cosh(%1.15e*(x-%f)) + %1.15e title '$\\chi^2/dof = %1.2f$'", massAnalyzerCombined->getFittedCoefficient(I2, 0), massAnalyzerCombined->getFittedMass(I2, 0), 0.5*LargestL, massAnalyzerCombined->getFittedVacuumExpectationValue(I2), massAnalyzerCombined->getFittedChiSquare(I2));
      plot->plotDirect(plotCmd);
      snprintf(plotCmd,1000,"replot %1.15e*cosh(%1.15e*(x-%f)) + %1.15e title '$\\chi^2/dof = %1.2f$'", massAnalyzerCombinedWithConstFit->getFittedCoefficient(I2, 0), massAnalyzerCombinedWithConstFit->getFittedMass(I2, 0), 0.5*LargestL, massAnalyzerCombinedWithConstFit->getFittedVacuumExpectationValue(I2), massAnalyzerCombinedWithConstFit->getFittedChiSquare(I2));
      plot->plotDirect(plotCmd);    
    } else {
      snprintf(plotCmd,1000,"replot %1.15e*cosh(%1.15e*(x-%f)) + %1.15e title '$\\chi^2/dof = %1.2f$'", massAnalyzerSeparate[I2]->getFittedCoefficient(0, 0), massAnalyzerSeparate[I2]->getFittedMass(0, 0), 0.5*LargestL,massAnalyzerSeparate[I2]->getFittedVacuumExpectationValue(0), massAnalyzerSeparate[I2]->getFittedChiSquare(0));
      plot->plotDirect(plotCmd);
      snprintf(plotCmd,1000,"replot %1.15e*cosh(%1.15e*(x-%f)) + %1.15e title '$\\chi^2/dof = %1.2f$'", massAnalyzerSeparateWithConstFit[I2]->getFittedCoefficient(0, 0), massAnalyzerSeparateWithConstFit[I2]->getFittedMass(0, 0), 0.5*LargestL,massAnalyzerSeparateWithConstFit[I2]->getFittedVacuumExpectationValue(0), massAnalyzerSeparateWithConstFit[I2]->getFittedChiSquare(0));
      plot->plotDirect(plotCmd);    
    }
  }
  
  for (int I=0; I<LargestL+1; I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  delete[] plotCmd;
  delete[] name;
  
  return plot;
}


LAPsystemPlot* EvaluateObservableCorrelatorBase::createPlot2(bool fitConst, bool combinedAna) {
  int LargestL = getLargestL(SDReader);
  MassCorrelationMatrixAnalyzer* massAna = massAnalyzerCombined;
  MassCorrelationMatrixAnalyzer** massAnaSep = massAnalyzerSeparate;
  if (fitConst) {
    massAna = massAnalyzerCombinedWithConstFit;
    massAnaSep = massAnalyzerSeparateWithConstFit;
  }

  LAPsystemPlot* plot;
  char* name = new char[2000];
  if ((fitConst) && (combinedAna)) snprintf(name,2000,"%sCAnaEffMaWFC", getObsName());
  if ((fitConst) && (!combinedAna)) snprintf(name,2000,"%sSAnaEffMaWFC", getObsName());
  if ((!fitConst) && (combinedAna)) snprintf(name,2000,"%sCAnaEffMa", getObsName());
  if ((!fitConst) && (!combinedAna)) snprintf(name,2000,"%sSAnaEffMa", getObsName());
  plot = LAPsystem->createNewPlot(name);  

  char* plotCmd = new char[1000];
  double** plotData = new double*[LargestL];
  for (int I=0; I<LargestL; I++) {
    plotData[I] = new double[1+2*NrOfIndependentVariables];
    plotData[I][0] = I+0.5;
    for (int I2=0; I2<NrOfIndependentVariables; I2++) {
      if (combinedAna) {
        plotData[I][2*I2+1] = massAna->getFittedEffectiveMasses(I, I2);
        plotData[I][2*I2+2] = massAna->getFittedEffectiveMassesError(I, I2);    
      } else {
        plotData[I][2*I2+1] = massAnaSep[I2]->getFittedEffectiveMasses(I, 0);
        plotData[I][2*I2+2] = massAnaSep[I2]->getFittedEffectiveMassesError(I, 0);    
      }
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
  if ((fitConst) && (combinedAna)) snprintf(name,2000,"Effective masses $m_H^{eff}(\\\\Delta t)$ from combined analysis of %s with fit to const", getObsName());
  if ((fitConst) && (!combinedAna)) snprintf(name,2000,"Effective masses $m_H^{eff}(\\\\Delta t)$ from separate analysis of %s with fit to const", getObsName());
  if ((!fitConst) && (combinedAna)) snprintf(name,2000,"Effective masses $m_H^{eff}(\\\\Delta t)$ from combined analysis of %s", getObsName());
  if ((!fitConst) && (!combinedAna)) snprintf(name,2000,"Effective masses $m_H^{eff}(\\\\Delta t)$ from separate analysis of %s", getObsName());
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
    if (combinedAna) {
      snprintf(plotCmd,1000,"replot %1.15e notitle", massAna->getFittedMass(I2, 0));
      plot->plotDirect(plotCmd);
    } else {
      snprintf(plotCmd,1000,"replot %1.15e notitle", massAnaSep[I2]->getFittedMass(0, 0));
      plot->plotDirect(plotCmd);    
    }
  }
  
  for (int I=0; I<LargestL; I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  delete[] plotCmd;
  delete[] name;
  
  return plot;
}


void EvaluateObservableCorrelatorBase::generateLatexAndPlotsAndXML() {
  if (doSeparateAnalysis) {
  
    LAPsystemPlot* plot1 = createPlot1(false, false);  
    LAPsystem->addPlot(plot1);
  
    LAPsystemPlot* plot2 = createPlot1(true, false);  
    LAPsystem->addPlot(plot2);  
    
    LAPsystemPlot* plot3 = createPlot2(false, false);  
    LAPsystem->addPlot(plot3);  

    LAPsystemPlot* plot4 = createPlot2(true, false);  
    LAPsystem->addPlot(plot4);  
    
    LAPsystem->addDirectText("\\begin{center}\n");
    LAPsystem->addDirectText("Separate mass analysis\n");
    LAPsystem->addDirectText("\\end{center}\n");
    
    startLatexOutputSummaryTable();

    char* tag = new char[2000];
    char* ext = new char[2000];
    
    for (int I2=0; I2<NrOfIndependentVariables; I2++) {
      if (I2>0) {
        snprintf(ext,2000,"SepA%d",I2);
      } else {
        ext[0] = (char) 0;
      }
      snprintf(tag,2000,"LatMass%s",ext);
      addXML_And_LatexOutputSummaryTableLine(tag, "Mass from correlation function in lattice units", "$m_{lat}$", massAnalyzerSeparate[I2]->getFittedMass(0, 0), massAnalyzerSeparate[I2]->getFittedMassError(0, 0), NULL, "%1.3f");
      
      snprintf(tag,2000,"PhysMass%s",ext);
      addXML_And_LatexOutputSummaryTableLine(tag, "Mass from correlation function in GeV", "$m$", physicalScaleInGEV*massAnalyzerSeparate[I2]->getFittedMass(0, 0), sqrt(sqr(physicalScaleInGEV*massAnalyzerSeparate[I2]->getFittedMassError(0, 0)) + sqr(physicalScaleErrorInGEV*massAnalyzerSeparate[I2]->getFittedMass(0, 0))), "GeV", "%1.1f");
    
      snprintf(tag,2000,"LatMassWFC%s",ext);
      addXML_And_LatexOutputSummaryTableLine(tag, "Mass from correlation function in lattice units with fit to const", "$m_{lat}$", massAnalyzerSeparateWithConstFit[I2]->getFittedMass(0, 0), massAnalyzerSeparateWithConstFit[I2]->getFittedMassError(0, 0), NULL, "%1.3f");
    
      snprintf(tag,2000,"FitConstWFC%s",ext);
      addXML_And_LatexOutputSummaryTableLine(tag, "Const from correlation function fit in lattice units", "$C$", massAnalyzerSeparateWithConstFit[I2]->getFittedVacuumExpectationValue(0), massAnalyzerSeparateWithConstFit[I2]->getFittedVacuumExpectationValueError(0), NULL, "%1.1e");
    
      snprintf(tag,2000,"PhysMassWFC%s",ext);
      addXML_And_LatexOutputSummaryTableLine(tag, "Mass from correlation function in GeV with fit to const", "$m$", physicalScaleInGEV*massAnalyzerSeparateWithConstFit[I2]->getFittedMass(0, 0), sqrt(sqr(physicalScaleInGEV*massAnalyzerSeparateWithConstFit[I2]->getFittedMassError(0, 0)) + sqr(physicalScaleErrorInGEV*massAnalyzerSeparateWithConstFit[I2]->getFittedMass(0, 0))), "GeV", "%1.1f");
    }

    delete[] tag;
    delete[] ext;
    endLatexOutputSummaryTable();
  }
  
  if (doMassCorrMatrixAnalysis) {
    LAPsystemPlot* plot1 = createPlot1(false, true);  
    LAPsystem->addPlot(plot1);
  
    LAPsystemPlot* plot2 = createPlot1(true, true);  
    LAPsystem->addPlot(plot2);  
    
    LAPsystemPlot* plot3 = createPlot2(false, true);  
    LAPsystem->addPlot(plot3);  

    LAPsystemPlot* plot4 = createPlot2(true, true);  
    LAPsystem->addPlot(plot4);  

    LAPsystem->addDirectText("\\begin{center}\n");
    LAPsystem->addDirectText("Mass correlation matrix analysis\n");
    LAPsystem->addDirectText("\\end{center}\n");

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
}
