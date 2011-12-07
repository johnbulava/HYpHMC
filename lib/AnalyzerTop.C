#include "AnalyzerTop.h"

void AnalyzerTop::ini(int l0, int l1, int l2, int l3, double y, ControlLogger* log, bool FLAG_GFit, char* fileNameExt, double tol) {
  totalN = 0;
  TOL = tol;
  TopTimeSliceCorrelatorDataCount = 0;
  Parameter_Y = y;
  Parameter_RHO = 1.0;
  Parameter_R = 0.5;
  FLAG_GnuplotFit = FLAG_GFit;
  L0 = l0;
  L1 = l1;
  L2 = l2;
  L3 = l3;  
  LatticeResult_PhysicalTopMass = NaN;
  LatticeResult_PhysicalTopMassError = NaN;
  TopMass = new AutoCorrelation(5,100);
  timeDirection = 0;
  confFileNameExtension = new char[1000];
  snprintf(confFileNameExtension,1000,"%s",fileNameExt);
  LargestL = L0;
  if (L1>LargestL) {
    LargestL = L1;
    timeDirection = 1;
  }
  if (L2>LargestL) {
    LargestL = L2;
    timeDirection = 2;
  }
  if (L3>LargestL) {
    LargestL = L3;
    timeDirection = 3;
  }

  MassControlLog = log;
  fermiOps = new FermionMatrixOperations(L0, L1, L2, L3, Parameter_RHO, Parameter_R, Parameter_Y);
  LeftVector  = fermiOps->createFermionVector();
  RightVector = fermiOps->createFermionVector();
  SolutionVector = fermiOps->createFermionVector();
  
  weightData = new double[AnalyzerTopDataMAX];
  sinPSqr = new double[L0*L1*L2*L3];
  TopTimeSliceCorrelatorData = new Complex*[AnalyzerTopDataMAX];
  TopTimeSliceCorrelatorConfNr = new int[AnalyzerTopDataMAX];
  confPoolIndOrder = new int[AnalyzerTopDataMAX];
  TopTimeSliceCorrelator = new AutoCorrelation*[1+LargestL];
  TopMassData = new double[AnalyzerTopDataMAX];

  TopMassAnalyzer = new MassCorrelationMatrixAnalyzer(1, LargestL, true, false, true, false, false, 1, 100000,"Top");
  
  int I;
  for (I=0; I<1+LargestL; I++) {
    TopTimeSliceCorrelator[I] = new AutoCorrelation(5,100);
  }
  
  int I0,I1,I2,I3;
  vector4D p;
  int count = 0;
  for (I0=0; I0<L0; I0++) {
    p[0] = 2*I0*pi/L0;
    for (I1=0; I1<L1; I1++) {
      p[1] = 2*I1*pi/L1;
      for (I2=0; I2<L2; I2++) {
        p[2] = 2*I2*pi/L2;
        for (I3=0; I3<L3; I3++) {
          p[3] = 2*I3*pi/L3;
	  
	  sinPSqr[count] = 4.0 * (sqr(sin(0.5*p[0])) + sqr(sin(0.5*p[1])) + sqr(sin(0.5*p[2])) + sqr(sin(0.5*p[3])));
	  count++;
	}
      }
    }
  }
}


void AnalyzerTop::desini() {
  delete[] sinPSqr;
  int I;
  for (I=0; I<1+LargestL; I++) {
    delete TopTimeSliceCorrelator[I];
  }
  for (I=0; I<TopTimeSliceCorrelatorDataCount; I++) {
    delete[] TopTimeSliceCorrelatorData[I];
  }  
  delete[] TopTimeSliceCorrelator;
  delete[] TopTimeSliceCorrelatorData;
  delete[] TopTimeSliceCorrelatorConfNr;
  delete[] confPoolIndOrder;
  delete[] weightData;
  delete TopMassAnalyzer;
  delete TopMass;
  delete[] TopMassData;
  fermiOps->destroyFermionVector(LeftVector);
  fermiOps->destroyFermionVector(RightVector);
  fermiOps->destroyFermionVector(SolutionVector);
  delete fermiOps;
  delete[] confFileNameExtension;
}


AnalyzerTop::AnalyzerTop(int l0, int l1, int l2, int l3, double y, ControlLogger* log, bool FLAG_GFit, char* fileNameExt, double tol) {
  ini(l0,l1,l2,l3,y, log, FLAG_GFit, fileNameExt,tol);
}


AnalyzerTop::~AnalyzerTop() {
  desini();
}


int AnalyzerTop::getTotalN() {
  return totalN;
}


void AnalyzerTop::sampleFermioncScanVector(Complex* v, int t, int FermionIndex) {
  int I;
  int VLxtr = fermiOps->getVectorLengthXtrSize();
  
  for (I=0; I<VLxtr; I++) {
    v[I].x = 0;
    v[I].y = 0;   
  }
  
  int ind[4];
  for (ind[0]=0; ind[0]<L0; ind[0]++) {
    for (ind[1]=0; ind[1]<L1; ind[1]++) {
      for (ind[2]=0; ind[2]<L2; ind[2]++) {
        for (ind[3]=0; ind[3]<L3; ind[3]++) {
	  if (ind[timeDirection] == t) {
	    int index = FermionIndex;
	    index += 8 * ind[3];
	    index += 8 * (L3+xtraSize3) * ind[2];
	    index += 8 * (L3+xtraSize3) * (L2+xtraSize2) * ind[1];
	    index += 8 * (L3+xtraSize3) * (L2+xtraSize2) * (L1+xtraSize1) * ind[0];
            v[index].x = 1.0;
	  }
	}
      }
    }
  }
}


void AnalyzerTop::calcPsiPsiBarMatrixForPhiField(vector4D* phiField, Complex**** PsiPsiBarMatrix) {
  bool Quse;
  double mu, beta;
  fermiOps->getQPreconditionerParameter(Quse, mu, beta);

  for (int t1=0; t1<LargestL; t1++) {
    for (int fInd1=0; fInd1<8; fInd1++) {
      if (LogLevel>2) printf("Psi-PsiBar-Matrix for t1 = %d, fInd1 = %d\n", t1, fInd1);          
      sampleFermioncScanVector(LeftVector, t1, fInd1);
      fermiOps->executeFermionMatrixMultiplication(LeftVector, LeftVector, (double*) phiField, false, NULL, NULL, 2, 1);
      int neededIter = 0;
      fermiOps->solveFermionMatrixLGS(LeftVector, SolutionVector, (double*) phiField, TOL, true, false, -1, neededIter);    
      if (Quse) {
        fermiOps->executeQPreconditionerMatrixMultiplication(SolutionVector, SolutionVector, false, false);
      }  

      for (int t2=0; t2<LargestL; t2++) {
        for (int fInd2=0; fInd2<8; fInd2++) {
          sampleFermioncScanVector(RightVector, t2, fInd2);
          Complex scalar(0,0);     
          SSE_ComplexScalarProduct(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, RightVector, SolutionVector, scalar);
          PsiPsiBarMatrix[t1][t2][fInd1][fInd2] = scalar;
        }   
      }
    }   
  }
}


void AnalyzerTop::analyzeHiggsField(vector4D* phiField, double weight, int confNr) {


printf("Spectrum Analysis...\n");
ComplexMatrix FermionMatrix(1);

fermiOps->constructNeubergerWithXiFermionMatrix(FermionMatrix, phiField, 0.5);
printf("Spectrum Analysis step 2...\n");

bool b = FermionMatrix.calcEigenvalues();
printf("eigenvalues are ok=%d\n",b);

FILE* file = fopen("eigenVofSplitYukawa.dat","w");
for (int i=0; i<FermionMatrix.matrixSize; i++) {
  fprintf(file,"%1.15e %1.15e\n",FermionMatrix.eigenvalues[i].x,FermionMatrix.eigenvalues[i].y);

}
fclose(file);
exit(0);


  if (LogLevel>2) printf("Top-Analyzing Higgs Field nr. %d with weight=%f and time direction = %d...\n",confNr,weight,timeDirection);

  int confPoolInd = confInPoolIndex(confNr);
  weightData[totalN] = weight;
  confPoolIndOrder[totalN] = confPoolInd;

  if (confPoolInd < 0) {
    //Top-Quark: Time-Slice Data
    TopTimeSliceCorrelatorData[TopTimeSliceCorrelatorDataCount] = new Complex[1+LargestL];
    TopTimeSliceCorrelatorConfNr[TopTimeSliceCorrelatorDataCount] = confNr;        
    confPoolIndOrder[totalN] = TopTimeSliceCorrelatorDataCount;
    
    int I;
    for (I=0; I<1+LargestL; I++) {
      TopTimeSliceCorrelatorData[TopTimeSliceCorrelatorDataCount][I].x = 0;
      TopTimeSliceCorrelatorData[TopTimeSliceCorrelatorDataCount][I].y = 0;    
    }



    Complex**** PsiPsiBarMatrix = new Complex***[LargestL];
    for (int t1=0; t1<LargestL; t1++) {
      PsiPsiBarMatrix[t1] = new Complex**[LargestL];
      for (int t2=0; t2<LargestL; t2++) {
        PsiPsiBarMatrix[t1][t2] = new Complex*[8];
        for (int fInd1=0; fInd1<8; fInd1++) {
          PsiPsiBarMatrix[t1][t2][fInd1] = new Complex[8];
          for (int fInd2=0; fInd2<8; fInd2++) {
            PsiPsiBarMatrix[t1][t2][fInd1][fInd2].x = NaN;
            PsiPsiBarMatrix[t1][t2][fInd1][fInd2].y = NaN;	  
	  }	  
	}
      }
    }

    fermiOps->setPreconditioner(true, 1.1, 0);
    fermiOps->setQPreconditioner(true, 0.25, 0.25);    
    calcPsiPsiBarMatrixForPhiField(phiField, PsiPsiBarMatrix);


    for (int t1=0; t1<LargestL; t1++) {
      for (int t2=0; t2<LargestL; t2++) {
        int deltaT = t1-t2;
        if (deltaT<0) deltaT = -deltaT;
        int deltaT2 = LargestL-deltaT; 
      
        Complex scalar(0,0); 
	for (I=0; I<8; I++) {
	  scalar = scalar + PsiPsiBarMatrix[t1][t2][I][I];
	}

        TopTimeSliceCorrelatorData[TopTimeSliceCorrelatorDataCount][deltaT].x += scalar.x;
        TopTimeSliceCorrelatorData[TopTimeSliceCorrelatorDataCount][deltaT].y += scalar.y;
        TopTimeSliceCorrelatorData[TopTimeSliceCorrelatorDataCount][deltaT2].x += scalar.x;
        TopTimeSliceCorrelatorData[TopTimeSliceCorrelatorDataCount][deltaT2].y += scalar.y;
      }   
    }


    for (int t1=0; t1<LargestL; t1++) {
      for (int t2=0; t2<LargestL; t2++) {
        for (int fInd1=0; fInd1<8; fInd1++) {
          delete[] PsiPsiBarMatrix[t1][t2][fInd1];
	}
        delete[] PsiPsiBarMatrix[t1][t2];	
      }
      delete[] PsiPsiBarMatrix[t1];
    }
    delete[] PsiPsiBarMatrix;


    
/*    int I2;
    int neededIter;
    fermiOps->setPreconditioner(true, 1.1, 0);
    for (I=0; I<LargestL; I++) {
      sampleFermioncScanVector(LeftVector, I, 1);
      fermiOps->executeFermionMatrixMultiplication(LeftVector, LeftVector, (double*) phiField, false, NULL, NULL, 2, 0);
      fermiOps->solveFermionMatrixLGS(LeftVector, SolutionVector, (double*) phiField, TOL, true, false, -1, neededIter);
      
      for (I2=0; I2<LargestL; I2++) {
        int deltaT = I-I2;
        if (deltaT<0) deltaT = -deltaT;
        int deltaT2 = LargestL-deltaT; 
      
        sampleFermioncScanVector(RightVector, I2, 1);
        Complex scalar;     
        SSE_ComplexScalarProduct(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, SolutionVector, RightVector, scalar);

        TopTimeSliceCorrelatorData[TopTimeSliceCorrelatorDataCount][deltaT].x += scalar.x;
        TopTimeSliceCorrelatorData[TopTimeSliceCorrelatorDataCount][deltaT].y += scalar.y;
        TopTimeSliceCorrelatorData[TopTimeSliceCorrelatorDataCount][deltaT2].x += scalar.x;
        TopTimeSliceCorrelatorData[TopTimeSliceCorrelatorDataCount][deltaT2].y += scalar.y;
      }   
    }*/


/*    for (I=0; I<LargestL; I++) {
      sampleFermioncScanVector(RightVector, I, 1);
      fermiOps->solveFermionMatrixLGS(RightVector, SolutionVector, (double*) phiField, TOL, false, false, -1, neededIter);
      
      for (I2=0; I2<LargestL; I2++) {
        int deltaT = I2-I;
        if (deltaT<0) deltaT = -deltaT;
        int deltaT2 = LargestL-deltaT; 
      
        sampleFermioncScanVector(LeftVector, I2, 1);
        Complex scalar;     
        SSE_ComplexScalarProduct(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, LeftVector, SolutionVector, scalar);

        TopTimeSliceCorrelatorData[TopTimeSliceCorrelatorDataCount][deltaT].x += scalar.x;
        TopTimeSliceCorrelatorData[TopTimeSliceCorrelatorDataCount][deltaT].y += scalar.y;
        TopTimeSliceCorrelatorData[TopTimeSliceCorrelatorDataCount][deltaT2].x += scalar.x;
        TopTimeSliceCorrelatorData[TopTimeSliceCorrelatorDataCount][deltaT2].y += scalar.y;
      }   
    }*/



    for (I=0; I<1+LargestL; I++) {
      int normFac = 2*LargestL;
      if ((I==0) || (I==LargestL)) normFac = LargestL;
      normFac *= L0*L1*L2*L3;
      normFac /= LargestL;    
      TopTimeSliceCorrelatorData[TopTimeSliceCorrelatorDataCount][I].x /= normFac;
      TopTimeSliceCorrelatorData[TopTimeSliceCorrelatorDataCount][I].y /= normFac;
    }

    TopTimeSliceCorrelatorDataCount++;
    saveTopTimeSliceCorrelator();
  }

  totalN++;  
}


void AnalyzerTop::plotTopTimeSliceCorrelator() {
  if (LogLevel>2) printf("Plotting Top-Time-Slice Correlator...\n");
  FILE* file = fopen("data/TimeSliceCorrTop.dat","w");
  int I;
  ComplexVector derivatives(5);
  derivatives.setZero();
  derivatives.vectorElements[1].x = 1;
  double* x = new double[1+LargestL];
  double* y = new double[1+LargestL];
  double* err = new double[1+LargestL];
  
  for (I=0; I<1+LargestL; I++) {
    x[I] = I;
    y[I] = TopTimeSliceCorrelator[I]->getAverage(1);
    err[I] = TopTimeSliceCorrelator[I]->estimateCombinedError(derivatives);
    double autocorr = TopTimeSliceCorrelator[I]->estimateAutoCorrelationTime();
    printf("Auto Correlation Time: %f\n",autocorr);
    fprintf(file,"%d %1.15f  %1.15f\n",I, y[I], err[I]);
  }

  char* fitCommand = new char[1000];
  if (FLAG_GnuplotFit) {
    char* functionBody = new char[1000];
    double* fitRes = new double[2];
    double* fitErr = new double[2];
    double redChiSqr = NaN;
  
    snprintf(functionBody,1000,"A1*cosh(A2*(x-%d))",LargestL/2);
    fitRes[0] = 1.0;
    fitRes[1] = 1.0;
    performGnuplotFit(functionBody, &(x[1]), &(y[1]), &(err[1]), LargestL-2, 2, fitRes, fitErr, redChiSqr);

    snprintf(fitCommand,1000,"replot %1.15f*cosh(%1.15f*(x-%d)) notitle",fitRes[0],fitRes[1],LargestL/2);
    printf("Determined Top Mass from Final Fit: %f +- %f\n",fitRes[1],fitErr[1]);
    
    delete[] functionBody;
    delete[] fitRes;
    delete[] fitErr;    
  } else {
    snprintf(fitCommand,1000,"\n");  
  }
  MassControlLog->addSection("Physical Top Mass");
  MassControlLog->addPlot("", NULL, "Time slice correlator", "$\\Delta t=|t_2-t_1|$", "$\\langle\\psi_{t_1}\\bar\\psi_{t_2}\\rangle$", x, y, err, NULL, NULL, 1+LargestL, fitCommand);
 
  
  fclose(file);
  delete[] fitCommand;
  delete[] x;
  delete[] y;
  delete[] err;
}


int TOPFIT_DataSet = -1;
int TOPFIT_DataCount = -1;
Complex** TOPFIT_Data = NULL;
double TOPFIT_fitQuality(double* p) {
  int I;
  double res = 0;
    
  for (I=1; I<TOPFIT_DataCount-1; I++) {
    double v = sqrt(sqr(TOPFIT_Data[TOPFIT_DataSet][I].x) + sqr(TOPFIT_Data[TOPFIT_DataSet][I].y));
    double f = p[0]*cosh(p[1]*(I-TOPFIT_DataCount/2));
    res += sqr(v-f);
  }
  
  return res;
}


void AnalyzerTop::calcTopTimeSliceCorrelator() {
  if (LogLevel>2) printf("Calculating Top-Time-Slice Correlators...\n");
  int I,I2;
  int* RunLengths = new int[1];
  RunLengths[0] = totalN;
  double* dummyData = new double[totalN];
  
  for (I=0; I<1+LargestL; I++) {
    for (I2=0; I2<totalN; I2++) {
      dummyData[I2] = sqrt(sqr(TopTimeSliceCorrelatorData[confPoolIndOrder[I2]][I].x) + sqr(TopTimeSliceCorrelatorData[confPoolIndOrder[I2]][I].y));
    }
    TopTimeSliceCorrelator[I]->loadData(1, RunLengths, weightData, dummyData);
  }
  
  TOPFIT_DataCount = LargestL;;
  TOPFIT_Data = TopTimeSliceCorrelatorData;
  for (I=0; I<totalN; I++) {
    TOPFIT_DataSet = confPoolIndOrder[I];
    double* pos = new double[2];
    double* bounds1 = new double[2];
    double* bounds2 = new double[2];
    bounds1[0] = 0;
    bounds1[1] = 0;
    bounds2[0] = 10;
    bounds2[1] = 10;
    double v1 = sqrt(sqr(TopTimeSliceCorrelatorData[TOPFIT_DataSet][1].x) + sqr(TopTimeSliceCorrelatorData[TOPFIT_DataSet][1].y));
    double v2 = sqrt(sqr(TopTimeSliceCorrelatorData[TOPFIT_DataSet][LargestL/2].x) + sqr(TopTimeSliceCorrelatorData[TOPFIT_DataSet][LargestL/2].y));

    pos[0] = v2 / cosh(0);
    pos[1] = log(2*v1/v2)/(LargestL/2);
    
    printf("Starting Fit Search from (%f,%f)\n",pos[0],pos[1]);
    bool b = GradientMinimization(&(TOPFIT_fitQuality), 2, 1E-3, 1E-6, 1E-8, pos, bounds1, bounds2, NULL, 3, 1000);
    if (b) {
      printf("Fit found at (%f,%f)\n",pos[0],pos[1]);    
      printf("Top mass from configuration %d in pool pos. %d: %f\n",TopTimeSliceCorrelatorConfNr[TOPFIT_DataSet],TOPFIT_DataSet,pos[1]);
      TopMassData[I] = pos[1];
    }   
    delete[] pos;
    delete[] bounds1;
    delete[] bounds2;     
  }
  
  TopMass->loadData(1, RunLengths, weightData, TopMassData);
  
  ComplexVector derivatives(5);
  derivatives.setZero();
  derivatives.vectorElements[1].x = 1;
  LatticeResult_PhysicalTopMass = TopMass->getAverage(1);
  LatticeResult_PhysicalTopMassError = TopMass->estimateCombinedError(derivatives);
  printf("Determined Top Mass: %f +- %f\n",LatticeResult_PhysicalTopMass,LatticeResult_PhysicalTopMassError);
  double autocorr = TopMass->estimateAutoCorrelationTime();
  printf("Auto Correlation Time: %f\n",autocorr);
  
  delete[] RunLengths;
  delete[] dummyData;  
}


int AnalyzerTop::confInPoolIndex(int confNr) {
  int I;
  for (I=0; I<TopTimeSliceCorrelatorDataCount; I++) {
    if (TopTimeSliceCorrelatorConfNr[I] == confNr) {
      return I;
    }
  }
  return -1;
}


void AnalyzerTop::loadTopTimeSliceCorrelator() {
  char* fileName = new char[1000];
  snprintf(fileName,1000,"%s/data/results/pHMC/analysis/TopTimeSliceCorrelatorTOL%1.2e%s.dat",DataBaseDirectory, TOL,confFileNameExtension);
  if (LogLevel>0) printf("Loding Top Slice Correlators from file: %s...\n",fileName);

  
  FILE* file = fopen(fileName,"r");
  if (file==NULL) {
    if (LogLevel>0) printf("Could not open Data File! ==> Starting Top Evaluation from scratch\n");
    return;
  }
  int I;
  bool endOfFile = false;
  while (!endOfFile) {
    TopTimeSliceCorrelatorData[TopTimeSliceCorrelatorDataCount] = new Complex[1+LargestL];
    if (fscanf(file,"%d ",&(TopTimeSliceCorrelatorConfNr[TopTimeSliceCorrelatorDataCount]))!=1) {
      endOfFile = true;
      break;        
    }
//TopTimeSliceCorrelatorConfNr[TopTimeSliceCorrelatorDataCount] = TopTimeSliceCorrelatorDataCount+1;
    
    for (I=0; I<1+LargestL; I++) {
      if (fscanf(file,"%lf %lf ", &(TopTimeSliceCorrelatorData[TopTimeSliceCorrelatorDataCount][I].x), &(TopTimeSliceCorrelatorData[TopTimeSliceCorrelatorDataCount][I].y))!=2) {
        endOfFile = true;
	break;
      }
    } 
    if (!endOfFile) TopTimeSliceCorrelatorDataCount++;    
    if (endOfFile) {
      delete[] TopTimeSliceCorrelatorData[TopTimeSliceCorrelatorDataCount];
    }
  }
  if (LogLevel>0) printf("Ready: %d data sets were read.\n",TopTimeSliceCorrelatorDataCount);

  fclose(file);
  delete[] fileName;
}


void AnalyzerTop::saveTopTimeSliceCorrelator() {
  char* fileName = new char[1000];
  snprintf(fileName,1000,"%s/data/results/pHMC/analysis/TopTimeSliceCorrelatorTOL%1.2e%s.dat",DataBaseDirectory, TOL,confFileNameExtension);
  if (LogLevel>0) printf("Saving Top Slice Correlators to file: %s...\n",fileName);
  FILE* file = fopen(fileName,"w");
  int I,I2;  
  for (I=0; I<TopTimeSliceCorrelatorDataCount; I++) {
    fprintf(file,"%d ", TopTimeSliceCorrelatorConfNr[I]);
    for (I2=0; I2<1+LargestL; I2++) {
      fprintf(file,"%1.15f %1.15f ",TopTimeSliceCorrelatorData[I][I2].x,TopTimeSliceCorrelatorData[I][I2].y);
    } 
    fprintf(file,"\n");  
  }
  fclose(file);
  delete[] fileName;

  if (LogLevel>0) printf("ready.\n");
}


void AnalyzerTop::plotTopMasses() {
  Complex*** opCorr = new Complex**[LargestL+1];
  for (int I=0; I<LargestL+1; I++) {
    opCorr[I] = new Complex*[1];
    opCorr[I][0] = new Complex[1];
  }
  for (int I=0; I<totalN; I++) {
    int ind = confPoolIndOrder[I];
    for (int t=0; t<LargestL+1; t++) {
      opCorr[t][0][0] =TopTimeSliceCorrelatorData[ind][t];
    }
    TopMassAnalyzer->addOperatorCorrelationData(TopTimeSliceCorrelatorConfNr[ind], weightData[ind], opCorr);
  }
  for (int I=0; I<LargestL+1; I++) {
    delete[] opCorr[I][0];
    delete[] opCorr[I];
  }
  delete[] opCorr;

  TopMassAnalyzer->plotEigenvalues();
  TopMassAnalyzer->plotEigenvalues(true);
  TopMassAnalyzer->plotEffectiveMasses();
  TopMassAnalyzer->plotEigenvalues(MassControlLog, true);
  TopMassAnalyzer->plotEffectiveMasses(MassControlLog);
  
  LatticeResult_PhysicalTopMass = TopMassAnalyzer->getFittedMass(0,0);
  LatticeResult_PhysicalTopMassError = TopMassAnalyzer->getFittedMassError(0,0);
  
}
