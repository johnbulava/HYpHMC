#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include <pthread.h>

#include "Tools.h"


#define GammaHFac 20

char* fileName = NULL;
int FitIterations = 0;
int dataCount = 0;
double** data;
double* dataX = NULL;
double* dataY = NULL;
double* dataXErr = NULL;
double* dataYErr = NULL;
double xmin = 0;
double ymin = 0;
double xmax = 0;
double ymax = 0;
double thresholdForShiftX = 0;
double thresholdForShiftY = 0;
double mG = 0;
int ScanPoints=10000;
int FitMode = 1;


void invalidParameters() {
  printf("Parameters must be: filename \n");
  exit(0);
}


void loadCommandLineParameters(int argc,char **argv) {
  bool error = false;
  if (argc<10) invalidParameters();
  fileName = new char[2000];
  snprintf(fileName, 2000, "%s", argv[1]);
  if (sscanf(argv[2],"%lf",&mG)!=1) error = true;
  if (sscanf(argv[3],"%lf",&xmin)!=1) error = true;
  if (sscanf(argv[4],"%lf",&ymin)!=1) error = true;
  if (sscanf(argv[5],"%lf",&xmax)!=1) error = true;
  if (sscanf(argv[6],"%lf",&ymax)!=1) error = true;
  if (sscanf(argv[7],"%lf",&thresholdForShiftX)!=1) error = true;
  if (sscanf(argv[8],"%lf",&thresholdForShiftY)!=1) error = true;
  if (sscanf(argv[9],"%d",&FitIterations)!=1) error = true;

  if (error) invalidParameters();
  
  printf("Read Parameters: \n");
  printf("  -> Filename           : %s\n", fileName);  
  printf("  -> mG                 : %f\n", mG);
  printf("  -> xmin               : %f\n", xmin);
  printf("  -> ymin               : %f\n", ymin);
  printf("  -> xmax               : %f\n", xmax);
  printf("  -> ymax               : %f\n", ymax);
  printf("  -> thresholdForShiftX : %f\n", thresholdForShiftX);
  printf("  -> thresholdForShiftY : %f\n", thresholdForShiftY);
  printf("  -> FitIterations      : %d\n", FitIterations);
  printf("  -> xmin               : %f\n", xmin);
}


void loadData() {
  printf("Loading file %s\n", fileName);
  int colMax = 4;
  
  FILE* file = fopen(fileName, "r");
  double dummy = 0;
  dataCount = 0;
  bool error = false;
  char* restLine = new char[3000];
  while (!error) {
    for (int I=0; I<colMax; I++) {
      if (fscanf(file, "%lf", &dummy) != 1) {
        error = true;
        break;
      }
    }
    dataCount++;
    fgets(restLine, 2500, file);
  }
  dataCount--;  
  fclose(file);

  data = new double*[dataCount];
  file = fopen(fileName, "r");
  dataCount = 0;
  error = false;
  while (!error) {
    data[dataCount] = new double[colMax];
    for (int I=0; I<colMax; I++) {
      if (fscanf(file, "%lf", &(data[dataCount][I])) != 1) {
        error = true;
        break;
      }
    }
    dataCount++;
    fgets(restLine, 2500, file);
  }
  dataCount--;
  fclose(file);
  
  delete[] restLine;
  
  printf("Read %d data lines\n", dataCount);
  dataX = new double[dataCount];
  dataY = new double[dataCount];
  dataXErr = new double[dataCount];
  dataYErr = new double[dataCount];
  for (int I=0; I<dataCount; I++) {
    dataX[I] = data[I][0];
    dataXErr[I] = data[I][1];
    dataY[I] = data[I][2];
    dataYErr[I] = data[I][3];

    printf(" **> %d: %1.3e %1.3e %1.3e %1.3e\n",I,dataX[I],dataXErr[I],dataY[I],dataYErr[I]);
    if ((dataX[I] > thresholdForShiftX) && (dataY[I] < thresholdForShiftY)) {
      dataY[I] += pi;
      printf(" --- SHIFTED DATUM ---\n");
      printf(" **> %d: %1.3e %1.3e %1.3e %1.3e\n",I,dataX[I],dataXErr[I],dataY[I],dataYErr[I]);
    }    
  }
}


void writeResultToDisk(double lamRen, double lamRenError, double mH, double mHError, double chiSqr) {
  FILE* file = fopen("ScatteringPhaseFitParameters.gnu", "a");
  fprintf(file, "mH%d = %1.15e\n", FitMode, mH);
  fprintf(file, "mHError%d = %1.15e\n", FitMode, mHError);
  if (FitMode==1) {
    fprintf(file, "lambdaRen%d = %1.15e\n", FitMode, lamRen);
    fprintf(file, "lambdaRenError%d = %1.15e\n", FitMode, lamRenError);
  } else {
    fprintf(file, "GammaH%d = %1.15e\n", FitMode, lamRen/GammaHFac);
    fprintf(file, "GammaHError%d = %1.15e\n", FitMode, lamRenError/GammaHFac);  
  }
  fprintf(file, "mG%d = %1.15e\n", FitMode, mG);
  fprintf(file, "thresholdForShiftX = %1.15e\n", thresholdForShiftX);
  fprintf(file, "thresholdForShiftY = %1.15e\n", thresholdForShiftY);
  fprintf(file, "ChirSqr%d = %1.15e\n", FitMode, chiSqr);
  
  fclose(file);
}


double fitFunction(double k, double lambdaRenGammaH, double mH) {
  double W = 2*sqrt(mG*mG + k*k);
  double WSqr = W*W;
  double mHSqr = mH*mH;
  double mGSqr = mG*mG;
  double kSqr = k*k;
  double delta = NaN;  
  
  if (FitMode==1) {
    double deltas = atan(lambdaRenGammaH*3*k*(mHSqr-mGSqr) / (48*pi*(mHSqr-WSqr)*W));
    double deltar = lambdaRenGammaH*(mHSqr-mGSqr)*log((4*kSqr+mHSqr)/mHSqr)/(96*pi*k*W) - lambdaRenGammaH*5*k/(48*pi*W);

    delta = deltas + deltar;
    if (delta<0) delta += pi;
  } else if (FitMode==2) {
    lambdaRenGammaH /= GammaHFac;
    double A = mHSqr*sqr(lambdaRenGammaH) / (mHSqr-2*mGSqr);
    delta = asin(sqrt(A*(WSqr-2*mGSqr) / (sqr(WSqr-mHSqr) + mHSqr*sqr(lambdaRenGammaH)) ));
    if (isNaN(delta)) delta = 0.5*pi;
    if (W>mH) delta = pi-delta;
  }
  return delta;
}


double calcChiSqr(double lambdaRen, double mH) {
  double chiSqr = 0;
  int count = 0;
  for (int I=0; I<dataCount; I++) {
    if ((dataX[I]>xmin) && (dataY[I]>ymin) && (dataX[I]<xmax) && (dataY[I]<ymax)) {
      count++;
      double smallest = -1;
      for (int I2=0; I2<=ScanPoints; I2++) {
        double k = (I2*(xmax-xmin))/ScanPoints + xmin;
	double d = fitFunction(k, lambdaRen, mH);
	double dummy = sqr((k-dataX[I])/dataXErr[I]) + sqr((d-dataY[I])/dataYErr[I]);
	if ((smallest<0) || (dummy<smallest)) smallest = dummy;
      }
      chiSqr += smallest;
    } 
  }
  chiSqr /= count;
  return chiSqr;
}


void plotFit(double lambdaRen, double mH, char* fileName) {
  FILE* file = NULL;
  file = fopen(fileName, "w");
  for (int I2=0; I2<=ScanPoints; I2++) {
    double k = (I2*(xmax-xmin))/ScanPoints + xmin;
    double d = fitFunction(k, lambdaRen, mH);
    fprintf(file, "%1.15e %1.15e \n", k, d);
  }
  
  fclose(file);
}


double performFitHelper(double* para) {
  return calcChiSqr(para[0], para[1]);
}


double performFit(double& lambdaRen, double& mH) {
  double pos[2];
  pos[0] = lambdaRen;
  pos[1] = mH;
  
printf("start: %f %f \n", lambdaRen, mH);  
  GradientMinimization(&performFitHelper, 2, 1E-6, 1E-6, 1E-5, pos, NULL, NULL, NULL, 3, 1000);
  
  lambdaRen = pos[0];
  mH = pos[1];
printf("end: %f %f \n", lambdaRen, mH);  
  
  return performFitHelper(pos);
}


void findErrorEstimatesOfFitParameters(double& lamRen, double& lamRenError, double& mH, double& mHError) {
  double* saveDataX = new double[dataCount];
  double* saveDataY = new double[dataCount];
  
  for (int I=0; I<dataCount; I++) {
    saveDataX[I] = dataX[I];
    saveDataY[I] = dataY[I];    
  }

  performFit(lamRen, mH);
  double startValLamRen = lamRen;
  double startValMH = mH;


  lamRen = 0;
  lamRenError = 0;
  mH = 0;
  mHError = 0;
  
  for (int I2=0; I2<FitIterations; I2++) {
    double z1, z2;
    for (int I=0; I<dataCount; I++) {
      AdvancedGaussZufall(AdvancedSeed, z1, z2);
      dataX[I] += z1 * dataXErr[I];
      dataY[I] += z2 * dataYErr[I];
    }
    
    double dummyLam, dummyMH;
    dummyLam = startValLamRen*exp(1.0*(AdvancedZufall(AdvancedSeed)-0.5));
    dummyMH = startValMH*exp(1.0*(AdvancedZufall(AdvancedSeed)-0.5));
    double chiSqr = performFit(dummyLam, dummyMH);
    if (FitMode==1) {
      printf("lambdaRen: %f , mH: %f, chiSqr: %f\n", dummyLam, dummyMH, chiSqr);    
    } else {
      printf("GammaH: %f , mH: %f, chiSqr: %f\n", dummyLam/GammaHFac, dummyMH, chiSqr);        
    }
    
    if (chiSqr<30) {
      lamRen += dummyLam;
      lamRenError += dummyLam*dummyLam;
      mH += abs(dummyMH);
      mHError += dummyMH*dummyMH;
    } else {
      printf("Fit-Result rejected due to chiScr>30!!!\n");
      I2--;
    }
  
    for (int I=0; I<dataCount; I++) {
      dataX[I] = saveDataX[I];
      dataY[I] = saveDataY[I];    
    }
  }
  
  lamRen /= FitIterations;
  lamRenError = sqrt(lamRenError/FitIterations - sqr(lamRen));
  mH /= FitIterations;
  mHError = sqrt(mHError/FitIterations - sqr(mH));
  
  delete[] saveDataX;
  delete[] saveDataY;
}



int main(int argc, char **argv) {
  loadCommandLineParameters(argc, argv);
  iniTools(5517);
  LogLevel = 0;
  loadData();  
  
  system("rm -f ScatteringPhaseFitParameters.gnu");
  
  FitMode = 1;
  double mH = 0.2;
  double mHError = 0;

  double lambdaRen = 0.2;
  double lambdaRenError = 0;

  double chiSqr = performFit(lambdaRen, mH);
  
  double lambdaRen2 = lambdaRen;
  double mH2 = mH;
  findErrorEstimatesOfFitParameters(lambdaRen2, lambdaRenError, mH2, mHError);
  
  writeResultToDisk(lambdaRen2, lambdaRenError, mH2, mHError, chiSqr);
  plotFit(lambdaRen, mH, "FitScatteringPhases1.dat");



  FitMode = 2;
  double GammaH = 0.015*GammaHFac;
  double GammaHError = 0;
  chiSqr = performFit(GammaH, mH);

  findErrorEstimatesOfFitParameters(GammaH, GammaHError, mH, mHError);
  chiSqr = performFit(GammaH, mH);

  findErrorEstimatesOfFitParameters(GammaH, GammaHError, mH, mHError);
  chiSqr = performFit(GammaH, mH);

  findErrorEstimatesOfFitParameters(GammaH, GammaHError, mH, mHError);
  chiSqr = performFit(GammaH, mH);

  
  writeResultToDisk(GammaH, GammaHError, mH, mHError, chiSqr);
  plotFit(GammaH, mH, "FitScatteringPhases2.dat");  
  plotFit(GammaH/2.0, mH, "FitScatteringPhases2Div2.dat");
  plotFit(GammaH/3.0, mH, "FitScatteringPhases2Div3.dat");
  plotFit(GammaH/4.0, mH, "FitScatteringPhases2Div4.dat");
  plotFit(2.0*GammaH, mH, "FitScatteringPhases2Mul2.dat");
  plotFit(3.0*GammaH, mH, "FitScatteringPhases2Mul3.dat");
  plotFit(4.0*GammaH, mH, "FitScatteringPhases2Mul4.dat");
  
  delete[] fileName;
  delete[] dataX;
  delete[] dataY;
  delete[] dataXErr;
  delete[] dataYErr;
}
