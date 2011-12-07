#include "EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling.h"

EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling::EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservable(aIOcon, sdr, "PhiRenormalizedUnsymmetricQuarticCoupling", "phirulam", relStart, relEnd) { 
  int L0 = sdr->getL0();  
  int L1 = sdr->getL1();  
  int L2 = sdr->getL2();  
  int L3 = sdr->getL3();  
  smallL0 = L0/8;
  smallL1 = L1/8;
  smallL2 = L2/8;
  smallL3 = L3/8;
  if (smallL0<=1) smallL0 = 2;
  if (smallL1<=1) smallL1 = 2;
  if (smallL2<=1) smallL2 = 2;
  if (smallL3<=1) smallL3 = 2;

  if ((L0!=L1) || (L1!=L2)) {
    printf("ERROR in EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling: Lattice geometry not supported!\n");
    exit(0);
  }

  latticeBins = new LatticeMomentumBins(L0,L1,L2,L3);  
  autoCorr = new AutoCorrelation(3, 100);
  avgLamRen = new double[smallL3*(smallL0-1)];
  sigmaLamRen = new double[smallL3*(smallL0-1)];
  sigmaLamRenHelper = new double[smallL3*(smallL0-1)];
  propTemp = new double[getAnalyzerResultsCount() / 2];  
  autoCorrelationTime = new double[smallL3*(smallL0-1)];
  pSqr = new double[smallL3*(smallL0-1)]; 
  
  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
}


EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling::~EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling() {
  delete latticeBins;
  delete autoCorr;
  delete[] propTemp;
  delete[] avgLamRen;  
  delete[] sigmaLamRen;
  delete[] sigmaLamRenHelper;
  delete[] autoCorrelationTime;
  delete[] pSqr;
}


int EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling::getIndex(int i0, int i1, int i2, int i3, int locInd) {
  if ((i0<1-smallL0) || (i0>smallL0-1)) {
    printf("ERROR in EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling::getIndex: invalid index: %d %d %d %d %d\n", i0, i1, i2, i3, locInd);
  }
  if ((i1<1-smallL1) || (i1>smallL1-1)) {
    printf("ERROR in EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling::getIndex: invalid index: %d %d %d %d %d\n", i0, i1, i2, i3, locInd);
  }
  if ((i2<1-smallL2) || (i2>smallL2-1)) {
    printf("ERROR in EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling::getIndex: invalid index: %d %d %d %d %d\n", i0, i1, i2, i3, locInd);
  }
  if ((i3<1-smallL3) || (i3>smallL3-1)) {
    printf("ERROR in EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling::getIndex: invalid index: %d %d %d %d %d\n", i0, i1, i2, i3, locInd);
  }
  if ((locInd<0) || (locInd>3)) {
    printf("ERROR in EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling::getIndex: invalid index: %d %d %d %d %d\n", i0, i1, i2, i3, locInd);
  }
  
  int p0 = i0+smallL0-1;  
  int p1 = i1+smallL1-1;    
  int p2 = i2+smallL2-1;      
  int p3 = i3+smallL3-1;
  int index = p3 + p2*(2*smallL3-1) + p1*(2*smallL2-1)*(2*smallL3-1) + p0*(2*smallL1-1)*(2*smallL2-1)*(2*smallL3-1);
  index *= 8;
  index += 2*locInd;
  return index;
}


void EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling::defineObsDependencies() { 
}


int EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling::getAnalyzerResultsCount() {
  return 8*(2*smallL0-1)*(2*smallL1-1)*(2*smallL2-1)*(2*smallL3-1);
}


void EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling::calcLamRen(int ignoreStart, int ignoreEnd) {
  ComplexVector derivatives(5);
  
  int vecs[4][4];
  int perm[6][3];
  
  perm[0][0] = 0;
  perm[0][1] = 1;
  perm[0][2] = 2;
  
  perm[1][0] = 0;
  perm[1][1] = 2;
  perm[1][2] = 1;

  perm[2][0] = 1;
  perm[2][1] = 0;
  perm[2][2] = 2;
  
  perm[3][0] = 1;
  perm[3][1] = 2;
  perm[3][2] = 0;
  
  perm[4][0] = 2;
  perm[4][1] = 0;
  perm[4][2] = 1;
  
  perm[5][0] = 2;
  perm[5][1] = 1;
  perm[5][2] = 0;


  double* data = new double[dataAvailCount];
  double* dataDet = new double[dataAvailCount];  
  int L0 = SDReader->getL0();  
  int L1 = SDReader->getL1();  
  int L2 = SDReader->getL2();  
  int L3 = SDReader->getL3();  
  double volume = L0*L1*L2*L3;


  for (int I=0; I<getAnalyzerResultsCount(); I+=2) {
    int count = 0;
    for (int I2=0; I2<dataAvailCount; I2++) {
      if ((I2<ignoreStart) || (I2>ignoreEnd)) {
        data[count] = dataAvail[I2][I+0]*dataAvail[I2][I+0] + dataAvail[I2][I+1]*dataAvail[I2][I+1];
        dataDet[count] = dataAvailWeightAndSign[I2];	
        count++;
      }
    }
    autoCorr->loadData(1, &(count), dataDet, data);  
    propTemp[I/2] = autoCorr->getAverage(1);
  }

  int avgCounter = 0;
  for (int nu=0; nu<smallL3; nu++) {
    for (int p=1; p<smallL0; p++) {
      avgLamRen[nu*(smallL0-1)+p-1] = 0;  
      autoCorrelationTime[nu*(smallL0-1)+p-1] = 0;
      pSqr[nu*(smallL0-1)+p-1] = 0;
      
      for (int a0=-1; a0<2; a0+=2) {
      for (int a1=-1; a1<1; a1+=2) {
      for (int a2=-1; a2<2; a2+=2) {
      for (int a3=-1; a3<1; a3+=2) {
      for (int a4=-1; a4<1; a4+=2) {
        vecs[0][0] = 0;
        vecs[0][1] = 0;
        vecs[0][2] = 0;
        vecs[0][3] = a0*nu;
		
        vecs[1][0] = 0;
        vecs[1][1] = a2*p;
        vecs[1][2] = -a3*p;
        vecs[1][3] = a1*nu;
	
        vecs[2][0] = a4*p;
        vecs[2][1] = -a2*p;
        vecs[2][2] = 0;
        vecs[2][3] = -a0*nu;
	
	for (int i=0; i<3; i++) vecs[3][i] = -(vecs[0][i] + vecs[1][i] + vecs[2][i]);
	
        for (int pI=0; pI<1; pI++) {
	  for (int locInd0=0; locInd0<4; locInd0++) {
	  for (int locInd1=0; locInd1<4; locInd1++) {
	  for (int locInd2=0; locInd2<4; locInd2++) {
	  for (int locInd3=0; locInd3<4; locInd3++) {
	    if (((locInd0==locInd1) && (locInd2==locInd3)) || ((locInd0==locInd2) && (locInd1==locInd3)) || ((locInd0==locInd3) && (locInd1==locInd2))) {
   	      int index0 = getIndex(vecs[0][perm[pI][0]],vecs[0][perm[pI][1]], vecs[0][perm[pI][2]],vecs[0][3], locInd0);
	      int index1 = getIndex(vecs[1][perm[pI][0]],vecs[1][perm[pI][1]], vecs[1][perm[pI][2]],vecs[1][3], locInd1);
	      int index2 = getIndex(vecs[2][perm[pI][0]],vecs[2][perm[pI][1]], vecs[2][perm[pI][2]],vecs[2][3], locInd2);
	      int index3 = getIndex(vecs[3][perm[pI][0]],vecs[3][perm[pI][1]], vecs[3][perm[pI][2]],vecs[3][3], locInd3);
	      
              int count = 0;
              for (int I2=0; I2<dataAvailCount; I2++) {
                if ((I2<ignoreStart) || (I2>ignoreEnd)) {
	          Complex c0 = Complex(dataAvail[I2][index0+0], dataAvail[I2][index0+1]);
	          Complex c1 = Complex(dataAvail[I2][index1+0], dataAvail[I2][index1+1]);
	          Complex c2 = Complex(dataAvail[I2][index2+0], dataAvail[I2][index2+1]);
	          Complex c3 = Complex(dataAvail[I2][index3+0], dataAvail[I2][index3+1]);
                  data[count] = (c0*c1*c2*c3).x;
                  dataDet[count] = dataAvailWeightAndSign[I2];	
                  count++;
                }
              }
              autoCorr->loadData(1, &(count), dataDet, data);  
              avgLamRen[nu*(smallL0-1)+p-1] += -(volume/24.0)*autoCorr->getAverage(1) / (propTemp[index0/2]*propTemp[index1/2]*propTemp[index2/2]*propTemp[index3/2]);
	      avgCounter++;
              autoCorrelationTime[nu*(smallL0-1)+p-1] = autoCorr->estimateAutoCorrelationTime();
              pSqr[nu*(smallL0-1)+p-1] = nu*(smallL0-1)+p;
	    }
	  }}}}
	}
      }}}}}
      avgLamRen[nu*(smallL0-1)+p-1] /= avgCounter;  
    }
  }

  delete[] data;
  delete[] dataDet;
}


bool EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling::evaluate() {
  int reps = 0;
  if (dataAvailCount <= 10) {
    reps = dataAvailCount;
  } else if (dataAvailCount <= 100) {
    reps = 10;
  } else {
    reps = ((int) sqrt(dataAvailCount));  
  }

reps = 3;

  double errorRescaleFactor =  sqrt(reps) * (1.0 - 1.0/reps);
  int blockSize = dataAvailCount / reps;
  if (blockSize<1) blockSize = 1;

  for (int I=0; I<smallL3*(smallL0-1); I++) {
    sigmaLamRen[I] = 0;
    sigmaLamRenHelper[I] = 0;    
  }

  for (int I=0; I<reps; I++) {
    int igStart = I * blockSize;
    int igEnd = (I+1) * blockSize;
    if (igEnd>=dataAvailCount) igEnd = dataAvailCount-1;
    calcLamRen(igStart, igEnd);

    for (int I2=0; I2<smallL3*(smallL0-1); I2++) {    
      sigmaLamRen[I2] += avgLamRen[I2] * avgLamRen[I2];
      sigmaLamRenHelper[I2] += avgLamRen[I2];
    }
  }
  
  calcLamRen(-1, -1);
  
  for (int I2=0; I2<smallL3*(smallL0-1); I2++) {
    sigmaLamRen[I2] = errorRescaleFactor * sqrt((sigmaLamRen[I2]/reps)  - sqr(sigmaLamRenHelper[I2]/reps));
  }

  return true;
}


LAPsystemPlot* EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling::createPlot1(double maxP) {
  char* name = new char[1000];
  snprintf(name,1000,"%smaxP%1.2f",getObsName(),maxP);
  LAPsystemPlot* plot = LAPsystem->createNewPlot(name);

  char* plotCmd = new char[1000];
  double** plotData = new double*[smallL3*(smallL0-1)];
  for (int I=0; I<smallL3*(smallL0-1); I++) {
    plotData[I] = new double[3];
    plotData[I][0] = pSqr[I];
    plotData[I][1] = avgLamRen[I];
    plotData[I][2] = sigmaLamRen[I];
  }
  plot->setPlotData(smallL3*(smallL0-1), 3, plotData);
  plot->setXLabel("$\\\\hat p^2$");
  plot->setYLabel("Renormalized quartic Phi coupling");
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
    
  for (int I=0; I<smallL3*(smallL0-1); I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  delete[] plotCmd;
  delete[] name;
  
  return plot;
}




void EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling::generateLatexAndPlotsAndXML() {
  LAPsystemPlot* plot1 = createPlot1(16);  
  LAPsystem->addPlot(plot1);

  LAPsystemPlot* plot2 = createPlot1(1);  
  LAPsystem->addPlot(plot2);
  

  startLatexOutputSummaryTable();
  
//  addXML_And_LatexOutputSummaryTableLine("LamRenP1", "Renormalized Phi quartic coupling at p=(0,0,0,1)", "$\\lambda^{\\Phi}_{ren}(0,0,0,1)$", avgLamRen[1], sigmaLamRen[1], NULL, "%1.5f");
  
  endLatexOutputSummaryTable();
}
