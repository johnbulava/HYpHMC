#include "EvaluateObservablePsiBarPsiCondensateBase.h"

EvaluateObservablePsiBarPsiCondensateBase::EvaluateObservablePsiBarPsiCondensateBase(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, char* oName, char* nick, double relStart, double relEnd) : EvaluateObservable(aIOcon, sdr, oName, nick, relStart, relEnd) {
  condensate.x = NaN;
  condensate.y = NaN; 
  condensateSigma.x = NaN;
  condensateSigma.y = NaN; 
  autoCorR = new AutoCorrelation(5, 100);
  autoCorI = new AutoCorrelation(5, 100);  
  numberOfRepeatedMeasurements = 5;
  displayedFormula = new char[1000];
  snprintf(displayedFormula, 1000, "\\bar\\psi\\psi");
}


EvaluateObservablePsiBarPsiCondensateBase::~EvaluateObservablePsiBarPsiCondensateBase() {
  delete autoCorR;
  delete autoCorI;  
  delete[] displayedFormula;
}


int EvaluateObservablePsiBarPsiCondensateBase::getAnalyzerResultsCount() {
  return 2*numberOfRepeatedMeasurements;  
}


void EvaluateObservablePsiBarPsiCondensateBase::defineObsDependencies() { 
}


bool EvaluateObservablePsiBarPsiCondensateBase::evaluate() {
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
  double* measureData = new double[numberOfRepeatedMeasurements*dataAvailCount];
  double* measureData2 = new double[numberOfRepeatedMeasurements*dataAvailCount];
  double* measureDataWeights = new double[numberOfRepeatedMeasurements*dataAvailCount];
  
  condensate.x = 0;
  condensate.y = 0;
  condensateSigma.x = 0;
  condensateSigma.y = 0;

  for (int I=0; I<reps; I++) {
    int igStart = I * blockSize;
    int igEnd = (I+1) * blockSize;
    if (igEnd>=dataAvailCount) igEnd = dataAvailCount-1;

    int count = 0;
    for (int I2=0; I2<dataAvailCount; I2++) {
      if ((I2<igStart) || (I2>igEnd)) {
        for (int I3=0; I3<numberOfRepeatedMeasurements; I3++) {
          measureData[count] = dataAvail[I2][2*I3+0];
          measureData2[count] = dataAvail[I2][2*I3+1];
          measureDataWeights[count] = dataAvailWeightAndSign[I2];
          count++;
        }
      }
    }

    autoCorR->loadData(1, &count, measureDataWeights, measureData);
    double dummy = autoCorR->getAverage(1);
    condensate.x += dummy;
    condensateSigma.x += sqr(dummy);
    
    autoCorI->loadData(1, &count, measureDataWeights, measureData2);
    dummy = autoCorI->getAverage(1);
    condensate.y += dummy;
    condensateSigma.y += sqr(dummy);
  }
  condensateSigma.x = errorRescaleFactor * sqrt(condensateSigma.x/reps - sqr(condensate.x/reps));
  condensateSigma.y = errorRescaleFactor * sqrt(condensateSigma.y/reps - sqr(condensate.y/reps));
  condensate.x /= reps;
  condensate.y /= reps;

  delete[] measureData;
  delete[] measureData2;
  delete[] measureDataWeights;

  return true;
}


LAPsystemPlot* EvaluateObservablePsiBarPsiCondensateBase::createPlot1(bool realPart) {
  LAPsystemPlot* plot = NULL;
  if (realPart) {
    plot = LAPsystem->createNewPlot("CondensateRealPart");
  } else {
    plot = LAPsystem->createNewPlot("CondensateImaginaryPart");
  }

  double** plotData = new double*[numberOfRepeatedMeasurements*dataAvailCount];
  int smallestID=0;
  int largestID=0;
  if (dataAvailCount>0) {
    smallestID=dataAvailID[0];
    largestID=dataAvailID[0];;  
  }
  int count = 0;
  int extraAdd = 0;
  if (!realPart) extraAdd = 1;
  for (int I=0; I<dataAvailCount; I++) {
    for (int I2=0; I2<numberOfRepeatedMeasurements; I2++) {
      plotData[count] = new double[2];
      plotData[count][0] = dataAvailID[I];
      plotData[count][1] = dataAvail[I][2*I2+extraAdd];
      if (dataAvailID[I]>largestID) largestID = dataAvailID[I];
      if (dataAvailID[I]<smallestID) smallestID = dataAvailID[I];
      count++;
    }
  }
  plot->setPlotData(numberOfRepeatedMeasurements*dataAvailCount, 2, plotData);
  plot->setXLabel("Monte Carlo time");
  if (realPart) {
    plot->setYLabel("Real part of fermion condensate");
  } else {
    plot->setYLabel("Imaginary part of fermion condensate");
  }
  plot->setCaption("Fermion Condensate");
  plot->setTitle("");
  plot->setYLogScale(false);
  plot->setSize(1.5, 1.5);
  plot->setXRange(smallestID, largestID);
  plot->setYErrorBars(false);  
  plot->setPointType(3);
  plot->setLineType(2);  
  plot->plotData("1:2");
  
  for (int I=0; I<numberOfRepeatedMeasurements*dataAvailCount; I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  
  return plot;
}


void EvaluateObservablePsiBarPsiCondensateBase::generateLatexAndPlotsAndXML() {
  LAPsystemPlot* plot1 = createPlot1(true);  
  LAPsystem->addPlot(plot1);
  
  LAPsystemPlot* plot2 = createPlot1(false);  
  LAPsystem->addPlot(plot2); 

  startLatexOutputSummaryTable();

  char* strtext = new char[1000];
  snprintf(strtext, 1000, "$Re(\\langle%s\\rangle)$", displayedFormula);

  addXML_And_LatexOutputSummaryTableLine("FermCondReal", "Real part of fermion condensate", strtext, condensate.x, condensateSigma.x, NULL, "%1.6f");

  snprintf(strtext, 1000, "$Im(\\langle%s\\rangle)$", displayedFormula);

  addXML_And_LatexOutputSummaryTableLine("FermCondImag", "Imaginary part of fermion condensate", strtext, condensate.y, condensateSigma.y, NULL, "%1.6f");

  delete[] strtext;

  endLatexOutputSummaryTable();
}
