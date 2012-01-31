#include "EvaluateObservableScalarPsiBarPsiVertexBase.h"

EvaluateObservableScalarPsiBarPsiVertexBase::EvaluateObservableScalarPsiBarPsiVertexBase(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, char* oName, char* nick, double relStart, double relEnd) : EvaluateObservable(aIOcon, sdr, oName, nick, relStart, relEnd) { 
  yr = NULL;
  yrError = NULL;
}


EvaluateObservableScalarPsiBarPsiVertexBase::~EvaluateObservableScalarPsiBarPsiVertexBase() {
  delete[] yr;
  delete[] yrError;
}


int EvaluateObservableScalarPsiBarPsiVertexBase::getAnalyzerResultsCount() {
  return 6 * getLargestL(SDReader) * getLargestL(SDReader);  
}


void EvaluateObservableScalarPsiBarPsiVertexBase::defineObsDependencies() { 
  addDependOnObsByName("GoldstonePropagator");
  addDependOnObsByName("TopBarTopPropagatorChiralLeftHandedGauged");
}


bool EvaluateObservableScalarPsiBarPsiVertexBase::evaluate() {
  int L0 = SDReader->getL0();
  int L1 = SDReader->getL1();
  int L2 = SDReader->getL2();
  int L3 = SDReader->getL3();
  double volume = L0*L1*L2*L3;
  int LargestL = getLargestL(SDReader);
  
  delete[] yr;
  delete[] yrError;
  yr = new double[LargestL*LargestL];
  yrError = new double[LargestL*LargestL];  
  for (int I2=0; I2<LargestL*LargestL; I2++) {
    yr[I2] = 0;
    yrError[I2] = 0;    
  }

  //For Jackknife analysis
  int jackIter = (int) (sqrt(dataAvailCount));
  if (dataAvailCount <= 0) jackIter = 0;
  if (dataAvailCount <= 10) jackIter = dataAvailCount;
  if (dataAvailCount <= 100) jackIter = 10;
  if (dataAvailCount >= 10000) jackIter = 100;  
  double errorRescale = sqrt(jackIter) * (1.0 - 1.0/jackIter);
  
  
  for (int I2=0; I2<LargestL*LargestL; I2++) {
    for (int iter=0; iter<jackIter; iter++) {
      int ignoreStart = (iter*dataAvailCount) / jackIter;
      int ignoreEnd = ((iter+1)*dataAvailCount) / jackIter;
      int effectiveDataCount = 0;
      Complex threePointCorr(0,0);
      Complex fermiProp(0,0);
      double hprop = 0;
      Complex hValPExp(0,0);
      
      if (I2==0) {
        for (int I=0; I<dataAvailCount; I++) {
          if ((I<ignoreStart) || (I>ignoreEnd)) {
	    effectiveDataCount++;
	    
            hValPExp = hValPExp + Complex(dataAvail[I][6*I2+4],dataAvail[I][6*I2+5]);
	  }
	}
        hValPExp = hValPExp / effectiveDataCount;
	effectiveDataCount = 0;
      }

      for (int I=0; I<dataAvailCount; I++) {
        if ((I<ignoreStart) || (I>ignoreEnd)) {
	  effectiveDataCount++;
	  
	  Complex hValP = Complex(dataAvail[I][6*I2+4],dataAvail[I][6*I2+5]) - hValPExp;
	  threePointCorr = threePointCorr + hValP * Complex(dataAvail[I][6*I2+0],dataAvail[I][6*I2+1]);
	  fermiProp = fermiProp + Complex(dataAvail[I][6*I2+2],dataAvail[I][6*I2+3]);
	  hprop += hValP.x*hValP.x + hValP.y*hValP.y;
	}
      }
      threePointCorr = threePointCorr / effectiveDataCount;
      fermiProp = fermiProp / effectiveDataCount;
      hprop = hprop / effectiveDataCount;
      Complex valC = -1.0 * sqrt(volume) * threePointCorr / (hprop * fermiProp);
            
      double val = valC.x;
      yr[I2] += val;
      yrError[I2] += val*val;
    }
    yr[I2] /= jackIter;
    yrError[I2] = errorRescale * sqrt(yrError[I2]/jackIter - sqr(yr[I2]));
  }  

  return true;
}


LAPsystemPlot* EvaluateObservableScalarPsiBarPsiVertexBase::createPlot1() {
  int LargestL = getLargestL(SDReader);

  char* name = new char[1000];
  snprintf(name,1000,"%s",getObsName());
  LAPsystemPlot* plot = LAPsystem->createNewPlot(name);
  delete[] name;

  double** plotData = new double*[LargestL*LargestL];
  for (int I=0; I<LargestL*LargestL; I++) {
    plotData[I] = new double[3];
    plotData[I][0] = I;
    plotData[I][1] = yr[I];
    plotData[I][2] = yrError[I];
  }
  plot->setPlotData(LargestL*LargestL, 3, plotData);
  plot->setXLabel("$\\\\Delta p$");
  plot->setYLabel("Vertex-Function without Z-Factors");
  plot->setCaption("Caption");
  plot->setTitle("");
  plot->setPlotTitle(getObsName());
  plot->setYLogScale(false);
  plot->setSize(0.8, 0.8);
  plot->setXRange(-0.5, 0.5*LargestL-0.5/*LargestL*LargestL*/);
  plot->setYErrorBars(true);  
  plot->setPointSize(0.35);
  plot->setPointType(5);
//  plot->setLineType(0);  
  plot->plotData("1:2:3");
  
  if (abs(SimParaSet->getKappaN()) > 1E-10) {
    char* fitString = new char[1000];
    snprintf(fitString,1000,"replot %1.3e", SimParaSet->getY0());
    plot->plotDirect(fitString);  
    delete[] fitString;
  }
  
  for (int I=0; I<LargestL*LargestL; I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  
  return plot;
}



void EvaluateObservableScalarPsiBarPsiVertexBase::generateLatexAndPlotsAndXML() {
  LAPsystemPlot* plot1 = createPlot1();  
  LAPsystem->addPlot(plot1);

  double ZG = ((EvaluateObservablePropagatorBase*)(getDependOnObsByName("GoldstonePropagator")))->getPropagatorZFactor();
  double ZG_Error = ((EvaluateObservablePropagatorBase*)(getDependOnObsByName("GoldstonePropagator")))->getPropagatorZFactorError();
  double ZPsi = ((EvaluateObservablePsiBarPsiPropagatorBase*)(getDependOnObsByName("TopBarTopPropagatorChiralLeftHandedGauged")))->getPropagatorZFactor();
  double ZPsi_Error = ((EvaluateObservablePsiBarPsiPropagatorBase*)(getDependOnObsByName("TopBarTopPropagatorChiralLeftHandedGauged")))->getPropagatorZFactorError();


  startLatexOutputSummaryTable();

  double yrp0 = sqrt(ZG)*ZPsi*yr[0];
  double yrp0_Error = abs(yrp0*sqrt(sqr(yrError[0]/yr[0]) + sqr(0.5*ZG_Error/ZG) + sqr(ZPsi_Error/ZPsi)));
  addXML_And_LatexOutputSummaryTableLine("Yrp0", "Renormalized Yukawa coupling at $\\Delta p=0$", "$y_r(\\Delta p=0)$",yrp0, yrp0_Error, NULL, "%1.3f");
  double yrp1 = sqrt(ZG)*ZPsi*yr[1];
  double yrp1_Error = abs(yrp1*sqrt(sqr(yrError[1]/yr[1]) + sqr(0.5*ZG_Error/ZG) + sqr(ZPsi_Error/ZPsi)));
  addXML_And_LatexOutputSummaryTableLine("Yrp1", "Renormalized Yukawa coupling at $\\Delta p=1$", "$y_r(\\Delta p=1)$",yrp1, yrp1_Error, NULL, "%1.3f");
  
  endLatexOutputSummaryTable();
  
  printf("%f %f\n",ZG , ZG_Error);
  printf("%f %f\n",ZPsi , ZPsi_Error);
  
}
