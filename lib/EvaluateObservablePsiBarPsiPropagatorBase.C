#include "EvaluateObservablePsiBarPsiPropagatorBase.h"
#include "EvaluateObservableMagnetizations.h"
#include "AutoCorrelation.h"

EvaluateObservablePsiBarPsiPropagatorBase::EvaluateObservablePsiBarPsiPropagatorBase(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, const char* oName, const char* nick, double relStart, double relEnd) : EvaluateObservable(aIOcon, sdr, oName, nick, relStart, relEnd) { 
  PropagatorValues = NULL;
  PropagatorErrors = NULL;
  FreePropagatorValues = NULL;
  pSqr = NULL;
  pHatSqr = NULL;
  
  redChiSqr_Lin = NaN;
  fittedMass_Lin = NaN;
  fittedZ_Lin = NaN;
  fitThreshold_Lin = 0.5;
  kickOutPHatThreshold_Lin = 1.0;
  fittedMass_Lin_Error = NaN;
  fittedZ_Lin_Error = NaN;

  redChiSqr_FreeAna = NaN;
  fittedMass_FreeAna = NaN;
  fittedZ_FreeAna = NaN;
  fitThreshold_FreeAna = 0.5;
  kickOutPHatThreshold_FreeAna = 1.0;
  fittedMass_FreeAna_Error = NaN;
  fittedZ_FreeAna_Error = NaN;
  
  RescaleFactor = 1.0;
}


EvaluateObservablePsiBarPsiPropagatorBase::~EvaluateObservablePsiBarPsiPropagatorBase() {
  delete[] PropagatorValues;
  delete[] PropagatorErrors;
  delete[] FreePropagatorValues;
  delete[] pSqr;
  delete[] pHatSqr;
}


int EvaluateObservablePsiBarPsiPropagatorBase::getAnalyzerResultsCount() {
  return 2*getLargestL(SDReader);  
}


void EvaluateObservablePsiBarPsiPropagatorBase::defineObsDependencies() { 
  addDependOnObsByName("Magnetizations");
}


double EvaluateObservablePsiBarPsiPropagatorBase::chiSqrForFreeAnaFit(double mass, double Zpsi) {
  int LargestL = getLargestL(SDReader);
  double chiSqr = 0;
  int countP = 0;

  for (int p=0; p<LargestL; p++) {
    if (pHatSqr[2*p + 0]<=fitThreshold_FreeAna) {
      vector4D pvec;
      pvec[0] = 0;
      pvec[1] = 0;
      pvec[2] = 0;
      pvec[3] = 2*pi*p/LargestL;
      
      double val = (freeInvPropAnalytical(pvec, mass, Zpsi)).x;
      chiSqr += sqr((val-PropagatorValues[2*p+0]) / PropagatorErrors[2*p+0]);
      countP++;
    }
  }
  
  return chiSqr / countP;
}


EvaluateObservablePsiBarPsiPropagatorBase* EvaluateObservablePsiBarPsiPropagatorBase_FreeAnaFitClassObjectPointer = NULL;
double EvaluateObservablePsiBarPsiPropagatorBase_FreeAnaFitMinFunction(double* x) {
  return EvaluateObservablePsiBarPsiPropagatorBase_FreeAnaFitClassObjectPointer->chiSqrForFreeAnaFit(x[0], x[1]);
}


bool EvaluateObservablePsiBarPsiPropagatorBase::evaluate() {
  int LargestL = getLargestL(SDReader);
  PropagatorValues = new double[2*LargestL];
  PropagatorErrors = new double[2*LargestL];
  FreePropagatorValues = new double[2*LargestL];
  pSqr = new double[2*LargestL];
  pHatSqr = new double[2*LargestL];

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

  double* measureData = new double[dataAvailCount];
  double* detData = new double[dataAvailCount];  
  for (int p=0; p<2*LargestL; p++) {
    PropagatorValues[p] = 0;
    PropagatorErrors[p] = 0;
  
    for (int I=0; I<reps; I++) {
      int igStart = I * blockSize;
      int igEnd = (I+1) * blockSize;
      if (igEnd>=dataAvailCount) igEnd = dataAvailCount-1;

      int count = 0;
      for (int I2=0; I2<dataAvailCount; I2++) {
        if ((I2<igStart) || (I2>igEnd)) {
          measureData[count] = RescaleFactor*dataAvail[I2][p];
          detData[count] = dataAvailWeightAndSign[I2];	
	  count++;
        }
      }
      AutoCorrelation* autoCorr = new AutoCorrelation(5,100);
      autoCorr->loadData(1, &count, detData, measureData);
      double dummy = 1.0 / autoCorr->getAverage(1);
      PropagatorValues[p] += dummy;
      PropagatorErrors[p] += dummy*dummy;
      delete autoCorr;
    }
  }

  for (int p=0; p<2*LargestL; p++) {
    PropagatorValues[p] /= reps;
    PropagatorErrors[p] = errorRescaleFactor * sqrt((PropagatorErrors[p]/reps) - sqr(PropagatorValues[p]));            

    int ip = p/2;
//    pSqr[p] = 4*sqr(sin(ip*pi/(LargestL)));  //\hat p squared
    pSqr[p] = sqr(sin(ip*2*pi/(LargestL)));  //Not \hat p but \tilde p squared
//    pSqr[p] = sqr((ip*2*pi/(LargestL)));  //p square
    pHatSqr[p] = 4*sqr(sin(ip*pi/LargestL));  //\hat p squared
  }
  
  //Fit with analytical free inverse Propagator
  double mag = ((EvaluateObservableMagnetizations*)(getDependOnObsByName("Magnetizations")))->getMagnetizationM();
  double* fitParameter = new double[2];
  fitParameter[0] = SDReader->getYN() * mag;
  fitParameter[1] = 1.0;  
  redChiSqr_FreeAna = 0;
  fittedMass_FreeAna = 0;
  fittedZ_FreeAna = 0;
  fittedMass_FreeAna_Error = 0;
  fittedZ_FreeAna_Error = 0;
  bool b0 = true;
  
  for (int I=0; I<reps; I++) {
    int igStart = I * blockSize;
    int igEnd = (I+1) * blockSize;
    if (igEnd>=dataAvailCount) igEnd = dataAvailCount-1;

    for (int p=0; p<2*LargestL; p++) {
      int count = 0;
      for (int I2=0; I2<dataAvailCount; I2++) {
        if ((I2<igStart) || (I2>igEnd)) {
          measureData[count] = RescaleFactor*dataAvail[I2][p];
          detData[count] = dataAvailWeightAndSign[I2];	
          count++;
        }
      }
      AutoCorrelation* autoCorr = new AutoCorrelation(5,100);
      autoCorr->loadData(1, &count, detData, measureData);
      double dummy = 1.0 / autoCorr->getAverage(1);
      PropagatorValues[p] = dummy;
      delete autoCorr;
    }
  
    fitParameter[0] = SDReader->getYN() * mag;
    fitParameter[1] = 1.0;  
    EvaluateObservablePsiBarPsiPropagatorBase_FreeAnaFitClassObjectPointer = this;
    b0 = b0 & GradientMinimization(EvaluateObservablePsiBarPsiPropagatorBase_FreeAnaFitMinFunction, 2, 1E-3, 1E-7, 1E-9, fitParameter, NULL, NULL, NULL, 3, 10000);
    
    redChiSqr_FreeAna += chiSqrForFreeAnaFit(fitParameter[0], fitParameter[1]);
    fittedMass_FreeAna += fitParameter[0];
    fittedZ_FreeAna += fitParameter[1];
    fittedMass_FreeAna_Error += sqr(fitParameter[0]);
    fittedZ_FreeAna_Error += sqr(fitParameter[1]);
  }
  delete[] fitParameter;
  
  redChiSqr_FreeAna /= reps;
  fittedMass_FreeAna /= reps;
  fittedZ_FreeAna /= reps;
  fittedMass_FreeAna_Error = errorRescaleFactor * sqrt(fittedMass_FreeAna_Error/reps - sqr(fittedMass_FreeAna));
  fittedZ_FreeAna_Error = errorRescaleFactor * sqrt(fittedZ_FreeAna_Error/reps - sqr(fittedZ_FreeAna));
  
  //Compute full averages
  for (int p=0; p<2*LargestL; p++) {
    int count = 0;
    for (int I2=0; I2<dataAvailCount; I2++) {
      measureData[count] = RescaleFactor*dataAvail[I2][p];
      detData[count] = dataAvailWeightAndSign[I2];	
      count++;
    }
    AutoCorrelation* autoCorr = new AutoCorrelation(5,100);
    autoCorr->loadData(1, &count, detData, measureData);
    double dummy = 1.0 / autoCorr->getAverage(1);
    PropagatorValues[p] = dummy;
    delete autoCorr;
  }
  

  //Linear Fit  
  char* functionBody = new char[1000];
  double* fitRes = new double[2];
  double* fitErr = new double[2];

  int incFit = 0;
  int pointsConsidered = 0;
  double* fitData_X = new double[LargestL];
  double* fitData_Y = new double[LargestL];
  double* fitData_E = new double[LargestL];

  snprintf(functionBody,1000,"(x+A2*A2)/(A1*A2)");
  
  bool b = true;
  fittedMass_Lin = 0;
  fittedMass_Lin_Error = 0;
  fittedZ_Lin = 0;
  fittedZ_Lin_Error = 0;
  redChiSqr_Lin = 0;
  for (int I=0; I<reps; I++) {
    int igStart = I * blockSize;
    int igEnd = (I+1) * blockSize;
    if (igEnd>=dataAvailCount) igEnd = dataAvailCount-1;

    incFit = 0;
    pointsConsidered = 0;
    for (int p=0; p<LargestL; p++) {
      if (pHatSqr[2*p + 0]<=kickOutPHatThreshold_Lin) {
        fitData_X[pointsConsidered] = pSqr[2*p + 0];
        fitData_Y[pointsConsidered] = 0;
        fitData_E[pointsConsidered] = PropagatorErrors[2*p + 0];
        if (fitData_X[pointsConsidered]<=fitThreshold_Lin) incFit++;

        int count = 0;
        for (int I2=0; I2<dataAvailCount; I2++) {
          if ((I2<igStart) || (I2>igEnd)) {
            measureData[count] = RescaleFactor*dataAvail[I2][2*p];
            detData[count] = dataAvailWeightAndSign[I2];	
	    count++;
          }
        }
        AutoCorrelation* autoCorr = new AutoCorrelation(5,100);
        autoCorr->loadData(1, &count, detData, measureData);
        fitData_Y[pointsConsidered] = 1.0 / autoCorr->getAverage(1);
        delete autoCorr;
	pointsConsidered++;
      }
    }
    
    fitRes[0] = 1.0;
    fitRes[1] = 1.0;
    double redChi;
    
    //Sorting data
    bool change = true;
    while (change) {
      change = false;
      for (int p=0; p<pointsConsidered-1; p++) {
        if (fitData_X[p]>fitData_X[p+1]) {
	  double dummy_X = fitData_X[p];
	  double dummy_Y = fitData_Y[p];
	  double dummy_E = fitData_E[p];
	  fitData_X[p] = fitData_X[p+1];
	  fitData_Y[p] = fitData_Y[p+1];
	  fitData_E[p] = fitData_E[p+1];
	  fitData_X[p+1] = dummy_X;
	  fitData_Y[p+1] = dummy_Y;
	  fitData_E[p+1] = dummy_E;
	  change = true;
	}
      }
    }
           
    
    b = b & performGnuplotFit(functionBody, fitData_X, fitData_Y, fitData_E, incFit, 2, fitRes, fitErr, redChi);
    fittedMass_Lin += fitRes[1];
    fittedMass_Lin_Error += sqr(fitRes[1]);    
    fittedZ_Lin += fitRes[0];
    fittedZ_Lin_Error += sqr(fitRes[0]);
    redChiSqr_Lin += redChi;
  }
    
  fittedMass_Lin /= reps;
  fittedMass_Lin_Error = errorRescaleFactor * sqrt((fittedMass_Lin_Error/reps) - sqr(fittedMass_Lin));
  fittedZ_Lin /= reps;
  fittedZ_Lin_Error = errorRescaleFactor * sqrt((fittedZ_Lin_Error/reps) - sqr(fittedZ_Lin));
  redChiSqr_Lin /= reps;

  delete[] fitData_X;
  delete[] fitData_Y;
  delete[] fitData_E;
  delete[] functionBody;
  delete[] fitRes;
  delete[] fitErr;
  delete[] measureData;
  delete[] detData;
    
  return b & b0;
}


LAPsystemPlot* EvaluateObservablePsiBarPsiPropagatorBase::createPlot1(int dataShift, double xRange) {
  int LargestL = getLargestL(SDReader);

  char* name = new char[1000];
  if (dataShift==0) {
    snprintf(name,1000,"%sRealMaxP%1.2f",getObsName(),xRange);
  } else {
    snprintf(name,1000,"%sImgMaxP%1.2f",getObsName(),xRange);
  }
  LAPsystemPlot* plot = LAPsystem->createNewPlot(name);
  delete[] name;

  double** plotData = new double*[LargestL*5];
  for (int I=0; I<LargestL*5; I++) {
    plotData[I] = new double[5];
    plotData[I][0] = -1;
    plotData[I][1] = -1;
    plotData[I][2] = -1;
    vector4D pvec;
    pvec[0] = 0;
    pvec[1] = 0;
    pvec[2] = 0;
    pvec[3] = 2*pi*I/(LargestL*20);
    plotData[I][3] = sqr(sin(pvec[3]));
    plotData[I][4] = (freeInvPropAnalytical(pvec, fittedMass_FreeAna, fittedZ_FreeAna)).x;
  }
  
  int pointsConsidered = 0;
  for (int I=0; I<LargestL; I++) {
    if (pHatSqr[2*I + dataShift] <= kickOutPHatThreshold_Lin) {
      plotData[pointsConsidered][0] = pSqr[2*I + dataShift];
      plotData[pointsConsidered][1] = PropagatorValues[2*I + dataShift];
      plotData[pointsConsidered][2] = PropagatorErrors[2*I + dataShift];
      pointsConsidered++;
    }
  }

  plot->setPlotData(LargestL*5, 5, plotData);
  plot->setXLabel("$\\\\tilde p^2$");
  plot->setYLabel("$\\\\langle \\\\psi_{p}\\\\bar\\\\psi_{p}\\\\rangle$");
  plot->setCaption("Caption");
  plot->setTitle("");
  plot->setPlotTitle(getObsName());
  plot->setYLogScale(false);
  plot->setSize(0.8, 0.8);
  plot->setXRange(0, xRange);
  plot->setYErrorBars(true);  
  plot->setPointSize(0.35);
  plot->setPointType(5);
//  plot->setLineType(0);  
  plot->plotData("1:2:3");

  plot->setPointType(1);
  plot->setYErrorBars(false);  
  char* title = new char[1000];
  snprintf(title, 1000, "$\\\\chi^2/dof = %1.2f$", redChiSqr_FreeAna);
  plot->setTitle(title);
  plot->plotData("4:5 with lines");
  delete[] title;


  if (dataShift==0) {
    char* fitString = new char[1000];
    snprintf(fitString,1000,"replot (x+%1.15e)/%1.15e  title '$\\chi^2/dof = %1.2f$'",sqr(fittedMass_Lin), fittedMass_Lin*fittedZ_Lin, redChiSqr_Lin);
    plot->plotDirect(fitString);  
    delete[] fitString;
  }
  
  for (int I=0; I<LargestL*5; I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  
  return plot;
}




void EvaluateObservablePsiBarPsiPropagatorBase::generateLatexAndPlotsAndXML() {
  LAPsystemPlot* plot1 = createPlot1(0, 0.80);  
  LAPsystem->addPlot(plot1);

  startLatexOutputSummaryTable();

  addXML_And_LatexOutputSummaryTableLine("Z", "Renormalization constant from propagator", "$Z$", fittedZ_FreeAna, fittedZ_FreeAna_Error, NULL, "%1.3f");  
  addXML_And_LatexOutputSummaryTableLine("LatMass", "Mass from propagator in lattice units", "$m_{lat}$", fittedMass_FreeAna, fittedMass_FreeAna_Error, NULL, "%1.3f");

  addXML_And_LatexOutputSummaryTableLine("ZLin", "Renormalization constant from propagator (Linear Fit)", "$Z$", fittedZ_Lin, fittedZ_Lin_Error, NULL, "%1.3f");  
  addXML_And_LatexOutputSummaryTableLine("LatMassLin", "Mass from propagator in lattice units (Linear Fit)", "$m_{lat}$", fittedMass_Lin, fittedMass_Lin_Error, NULL, "%1.3f");

  addXML_And_LatexOutputSummaryTableLine("ZFreeAna", "Renormalization constant from propagator (Fit with free prop)", "$Z$", fittedZ_FreeAna, fittedZ_FreeAna_Error, NULL, "%1.3f");  
  addXML_And_LatexOutputSummaryTableLine("LatMassFreeAna", "Mass from propagator in lattice units (Fit with free prop)", "$m_{lat}$", fittedMass_FreeAna, fittedMass_FreeAna_Error, NULL, "%1.3f");
  
  endLatexOutputSummaryTable();
}


double EvaluateObservablePsiBarPsiPropagatorBase::getPropagatorZFactor() {
  return fittedZ_FreeAna;
}


double EvaluateObservablePsiBarPsiPropagatorBase::getPropagatorZFactorError() {
  return fittedZ_FreeAna_Error;
}


double EvaluateObservablePsiBarPsiPropagatorBase::getPropagatorMass() {
  return fittedMass_FreeAna;
}


double EvaluateObservablePsiBarPsiPropagatorBase::getPropagatorMassError() {
  return fittedMass_FreeAna_Error;
}


Complex EvaluateObservablePsiBarPsiPropagatorBase::freeInvPropAnalytical(vector4D pvec, double mass, double Zpsi) {
  NeubergerMatrix* ovDiracOp = new NeubergerMatrix(SDReader->getRho(), SDReader->getR(), 1,1,1,1,1);

  ComplexVector ev[4];
  ev[0].resize(4);
  ev[1].resize(4);
  ev[2].resize(4);
  ev[3].resize(4);
  ovDiracOp->analyticalEigenvectors(pvec, ev);
  Complex nup = ovDiracOp->analyticalEigenvalue(pvec);
  Complex num = adj(nup);    
  Complex ewp = Complex(1,0) / (nup + mass*(Complex(1,0) - (0.5/SDReader->getRho())*nup));
  Complex ewm = adj(ewp);
  ComplexMatrix matMinv = ewp*ComplexMatrix(ev[0]);
  matMinv = (ewp*ComplexMatrix(ev[1])) + matMinv;
  matMinv = (ewm*ComplexMatrix(ev[2])) + matMinv;
  matMinv = (ewm*ComplexMatrix(ev[3])) + matMinv;

  ComplexMatrix Dop = nup*ComplexMatrix(ev[0]);
  Dop = (nup*ComplexMatrix(ev[1])) + Dop;
  Dop = (num*ComplexMatrix(ev[2])) + Dop;
  Dop = (num*ComplexMatrix(ev[3])) + Dop;
    
  Complex Tres = Complex(0,0);
  for (int I=0; I<4; I++) {
    Complex nu = nup;
    if (I>=2) nu = num;
    
    ComplexVector dummy = ev[I]; 
    dummy.vectorElements[0] = Complex(0,0);
    dummy.vectorElements[1] = Complex(0,0);
      
    ComplexVector dummy2 = matMinv * dummy;
      
    ComplexVector dummy3 = 0.5 * (1.0/SDReader->getRho()) * (Dop * dummy2);
    dummy3.vectorElements[2] = dummy2.vectorElements[2] - dummy3.vectorElements[2];
    dummy3.vectorElements[3] = dummy2.vectorElements[3] - dummy3.vectorElements[3];

    for (int I2=0; I2<4; I2++) {
      Tres = Tres + adj(ev[I].vectorElements[I2]) * dummy3.vectorElements[I2];
    }

  }
  Complex invProp = Complex(1.0/(RescaleFactor*Zpsi),0) / Tres;
  delete ovDiracOp;
  
  return invProp;
}
