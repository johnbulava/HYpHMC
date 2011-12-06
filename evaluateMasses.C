#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <stddef.h>
#include <sys/types.h>
#include <dirent.h>

#include "Complex.h"
#include "Quat.h"
#include "Global.h"
#include "BootStrapClass.h"
#include "ControlLogger.h"
#include "AutoCorrelation.h"
#include "AnalyzerHiggs.h"
#include "AnalyzerTop.h"
#include "LatticeMomentumBins.h"

#define Critical_Nu 0.5
#define Critical_Gamma 1.2  
#define Physical_VEV_GeV 246
#define Physical_TopMass_GeV 175

  
//Variables
double* phiField;
double exactWeight;
AnalyzerTop* analyzerTop;
AnalyzerHiggs* analyzerHiggs;
char* confFileNameExtension;
ControlLogger MassControlLog;
bool FLAG_TexLog = true;
bool FLAG_GnuplotFit = true;
double HiggsFieldScale_C = NaN;
int Parameter_L0;
int Parameter_L1;
int Parameter_L2;
int Parameter_L3;
int Parameter_Nf;
double Parameter_KappaN;
double Parameter_LambdaN;
double Parameter_YN;
double Parameter_M0Sqr;
double Parameter_Lambda0;
double Parameter_Y0;
double Parameter_RHO;
double Parameter_R;
double Parameter_TopCorrTOL;
int Parameter_PolyDegree;
double Parameter_PolyAlpha;
char* Parameter_filenameSuffix;	 
double Parameter_MaxRunTime;
double startTime;



void startTimer() {
  startTime = clock();
}


double timePassed() {
  double time = (clock()-startTime)/CLOCKS_PER_SEC;
  return time;
}


bool timeOver() {
  if (timePassed()>Parameter_MaxRunTime*3600) return true;
  return false;
}



void buildOutputFileNameExtension() {
  confFileNameExtension = new char[300];
  
  snprintf(confFileNameExtension,300,"L%dx%dx%dx%dNf%dKap%1.5fLam%1.5fY%1.5fRho%1.3fR%1.3fPolDeg%dPolAl%1.3f_%s",
   Parameter_L0,Parameter_L1,Parameter_L2,Parameter_L3,Parameter_Nf,Parameter_KappaN,Parameter_LambdaN,Parameter_YN,Parameter_RHO,Parameter_R,Parameter_PolyDegree,Parameter_PolyAlpha,Parameter_filenameSuffix);	 
}  



void ini() {
  iniTools(5517);
//  optimizeFermionVectorEmbedding(Parameter_L0, Parameter_L1, Parameter_L2, Parameter_L3,0,1,1,1,1);
  buildOutputFileNameExtension(); 
  MassControlLog.setLogging(FLAG_TexLog);
  char* LogFileName = new char[1000];
  snprintf(LogFileName,1000,"MassControlLog");
  MassControlLog.setLatexFileName(LogFileName);
  delete[] LogFileName;
  MassControlLog.setPLOT_PointSize(0.2);
  MassControlLog.setPLOT_FrameSize(1.0, 1.0);
  exactWeight = NaN;
  phiField = new double[4*Parameter_L0*Parameter_L1*Parameter_L2*Parameter_L3];

  HiggsFieldScale_C = 1.0 / sqrt(2*Parameter_KappaN);
  Parameter_M0Sqr = (1.0 - 2.0*Parameter_Nf*Parameter_LambdaN - 8.0*Parameter_KappaN) / Parameter_KappaN;
  Parameter_Lambda0 = sqr(sqr(HiggsFieldScale_C)) * Parameter_LambdaN;
  Parameter_Y0 = HiggsFieldScale_C * Parameter_YN;

  analyzerHiggs = new AnalyzerHiggs(Parameter_L0, Parameter_L1, Parameter_L2, Parameter_L3, &MassControlLog, FLAG_GnuplotFit);
  analyzerTop = new AnalyzerTop(Parameter_L0, Parameter_L1, Parameter_L2, Parameter_L3, Parameter_YN, &MassControlLog, FLAG_GnuplotFit, confFileNameExtension, Parameter_TopCorrTOL);
  startTimer();
}


bool loadPhiconfiguration(char* fileName, bool rescale) {
  if (rescale) {
    if (LogLevel>2) printf("Loading rescaled Configuration: %s...",fileName);
  } else {
    if (LogLevel>2) printf("Loading Configuration: %s...",fileName);  
  }
  std::fstream confFile;
  confFile.open(fileName, std::ios::in);
  if (!confFile.good()) {
    if (LogLevel>2) printf("ERROR !!!\n");
    return false;
  }

  confFile.read((char*)phiField, 32*Parameter_L0*Parameter_L1*Parameter_L2*Parameter_L3);
  if (confFile.eof()) {
    if (LogLevel>2) printf("ERROR !!!\n");
    return false;
  }

  confFile.read((char*)(&exactWeight),8);
  if (!confFile.eof()) {
    if (LogLevel>2) printf(" with weight: %f ",exactWeight);
  } else {
    if (LogLevel>2) printf(" without weight ");
    exactWeight = 1.0;
  }

  confFile.close();

  if (rescale) {
    int I;
    for (I=0; I<4*Parameter_L0*Parameter_L1*Parameter_L2*Parameter_L3; I++) {
      phiField[I] /= HiggsFieldScale_C;
    }
  }

  if (LogLevel>2) printf("successfully.\n");
  return true;
}  


bool loadPhiconfiguration(int fileNr, bool rescale) {
  char* fileName = new char[1000];
  snprintf(fileName,1000,"%s/data/results/pHMC/configurations/subFolder%s/PhiConf%s_%d.dat",DataBaseDirectory,confFileNameExtension,confFileNameExtension,fileNr);
  
  bool b = loadPhiconfiguration(fileName, rescale);
  
  delete[] fileName;
  return b;
}


void multiplyHiggsFieldWithConst(double fac, double* phiField) {
  for (int I=0; I<4*Parameter_L0*Parameter_L1*Parameter_L2*Parameter_L3; I++) {
    phiField[I] *= fac;  
  }
}


void getHiggsFieldDirection(vector4D* phiField, vector4D dir) {
  dir[0] = 0;
  dir[1] = 0;
  dir[2] = 0;
  dir[3] = 0;
  
  for (int I=0; I<Parameter_L0*Parameter_L1*Parameter_L2*Parameter_L3; I++) {
    dir[0] += phiField[I][0];
    dir[1] += phiField[I][1];
    dir[2] += phiField[I][2];
    dir[3] += phiField[I][3];
  }
  
  dir[0] /= (Parameter_L0*Parameter_L1*Parameter_L2*Parameter_L3);
  dir[1] /= (Parameter_L0*Parameter_L1*Parameter_L2*Parameter_L3);
  dir[2] /= (Parameter_L0*Parameter_L1*Parameter_L2*Parameter_L3);
  dir[3] /= (Parameter_L0*Parameter_L1*Parameter_L2*Parameter_L3);
  if (LogLevel>2) {
    printf("Higgs-Field-Direction: (%f,%f,%f,%f)\n",dir[0],dir[1],dir[2],dir[3]);
  }
}


void rotateHiggsField(vector4D* phiField, int ind1, int ind2, double w) {
  double c = cos(w);
  double s = sin(w);
  for (int I=0; I<Parameter_L0*Parameter_L1*Parameter_L2*Parameter_L3; I++) {
    double x = phiField[I][ind1];
    double y = phiField[I][ind2];
    phiField[I][ind1] = c*x+s*y;
    phiField[I][ind2] = -s*x+c*y;    
  }
}


void alignHiggsFieldDirection(double* phiField) {
  vector4D dir;
  double w;
  getHiggsFieldDirection((vector4D*) phiField, dir);
  w = getAngle(dir[1], dir[0]);
  rotateHiggsField((vector4D*) phiField, 1, 0, -w);

  getHiggsFieldDirection((vector4D*) phiField, dir);
  w = getAngle(dir[2], dir[0]);
  rotateHiggsField((vector4D*) phiField, 2, 0, -w);

  getHiggsFieldDirection((vector4D*) phiField, dir);
  w = getAngle(dir[3], dir[0]);
  rotateHiggsField((vector4D*) phiField, 3, 0, -w);

  
  getHiggsFieldDirection((vector4D*) phiField, dir);
  
  
  
}



void printPhysicalAnalysis() {
  double Physical_Lambda_GeV = Physical_VEV_GeV*sqrt(analyzerHiggs->LatticeResult_GoldstoneZFactor) / analyzerHiggs->LatticeResult_VEV;
  double Physical_LambdaError_GeV = Physical_Lambda_GeV*sqrt(sqr(0.5*analyzerHiggs->LatticeResult_GoldstoneZFactorSigma/analyzerHiggs->LatticeResult_GoldstoneZFactor) + sqr(analyzerHiggs->LatticeResult_VEVsigma/analyzerHiggs->LatticeResult_VEV)); 
  double Physical_PhysicalTopMass_Gev = Physical_Lambda_GeV * analyzerTop->LatticeResult_PhysicalTopMass;
  double Physical_PhysicalTopMassError_Gev = sqrt(sqr(Physical_Lambda_GeV * analyzerTop->LatticeResult_PhysicalTopMassError) + sqr(Physical_LambdaError_GeV*analyzerTop->LatticeResult_PhysicalTopMass));
  double Physical_HiggsPropagatorMass_Gev = Physical_Lambda_GeV * analyzerHiggs->LatticeResult_HiggsPropagatorMass;
  double Physical_HiggsPropagatorMassError_Gev = sqrt(sqr(Physical_Lambda_GeV * analyzerHiggs->LatticeResult_HiggsPropagatorMassSigma) + sqr(Physical_LambdaError_GeV*analyzerHiggs->LatticeResult_HiggsPropagatorMass));
  double Physical_PhysicalHiggsMass_GeV = Physical_Lambda_GeV * analyzerHiggs->LatticeResult_PhysicalHiggsMass;
  double Physical_PhysicalHiggsMassError_GeV = sqrt(sqr(Physical_Lambda_GeV * analyzerHiggs->LatticeResult_PhysicalHiggsMassSigma) + sqr(Physical_LambdaError_GeV*analyzerHiggs->LatticeResult_PhysicalHiggsMass));
  double Physical_RenormalizedLambda = 0.5*sqr(analyzerHiggs->LatticeResult_PhysicalHiggsMass*sqrt(analyzerHiggs->LatticeResult_GoldstoneZFactor)/analyzerHiggs->LatticeResult_VEV);
  double Physical_RenormalizedLambdaError = Physical_RenormalizedLambda*sqrt(sqr(analyzerHiggs->LatticeResult_PhysicalHiggsMassSigma/analyzerHiggs->LatticeResult_PhysicalHiggsMass) + sqr(analyzerHiggs->LatticeResult_GoldstoneZFactorSigma/analyzerHiggs->LatticeResult_GoldstoneZFactor) + sqr(analyzerHiggs->LatticeResult_VEVsigma/analyzerHiggs->LatticeResult_VEV));
  double Physical_RenormalizedY = sqrt(analyzerHiggs->LatticeResult_GoldstoneZFactor)*analyzerTop->LatticeResult_PhysicalTopMass/analyzerHiggs->LatticeResult_VEV;
  double Physical_RenormalizedYError = Physical_RenormalizedY*sqrt(sqr(analyzerHiggs->LatticeResult_GoldstoneZFactorSigma/analyzerHiggs->LatticeResult_GoldstoneZFactor) + sqr(analyzerTop->LatticeResult_PhysicalTopMassError/analyzerTop->LatticeResult_PhysicalTopMass) + sqr(analyzerHiggs->LatticeResult_VEVsigma/analyzerHiggs->LatticeResult_VEV));

  printf("Lattice vev (rescaled):                     %1.5f +- %1.5f\n",analyzerHiggs->LatticeResult_VEV,analyzerHiggs->LatticeResult_VEVsigma);
  printf("Lattice Goldstone Z-Faktor (rescaled):      %1.5f +- %1.5f\n",analyzerHiggs->LatticeResult_GoldstoneZFactor,analyzerHiggs->LatticeResult_GoldstoneZFactorSigma);
  printf("Lattice Top-Mass:                           %1.5f +- %1.5f\n",analyzerTop->LatticeResult_PhysicalTopMass, analyzerTop->LatticeResult_PhysicalTopMassError);
  printf("Lattice Higgs Propagator Mass:              %1.5f +- %1.5f\n",analyzerHiggs->LatticeResult_HiggsPropagatorMass,analyzerHiggs->LatticeResult_HiggsPropagatorMassSigma);
  printf("Lattice Physical HiggsMass:                 %1.5f +- %1.5f\n",analyzerHiggs->LatticeResult_PhysicalHiggsMass,analyzerHiggs->LatticeResult_PhysicalHiggsMassSigma);
  printf("Lattice Size:                        %dx%dx%dx%d\n",Parameter_L0,Parameter_L1,Parameter_L2,Parameter_L3);
  printf("Lattice KappaN:                      %1.5f\n",Parameter_KappaN);
  printf("Lattice LambdaN:                     %1.5f\n",Parameter_LambdaN);
  printf("Lattice YN:                          %1.5f\n",Parameter_YN);



  printf("Rescaling Factor C:                      %1.5f\n",HiggsFieldScale_C);
  printf("Continuum Notation m0^2:                 %1.5f\n",Parameter_M0Sqr);
  printf("Continuum Notation y0:                   %1.5f\n",Parameter_Y0);
  printf("Continuum Notation lambda0:              %1.5f\n",Parameter_Lambda0);



  printf("Physical Cut-off:             (%1.3f +- %1.3f) GeV\n",Physical_Lambda_GeV, Physical_LambdaError_GeV);
  printf("Physical Top-Mass:            (%1.3f +- %1.3f) GeV\n",Physical_PhysicalTopMass_Gev, Physical_PhysicalTopMassError_Gev);
  printf("Physical Prop. Higgs-Mass:    (%1.3f +- %1.3f) GeV\n",Physical_HiggsPropagatorMass_Gev, Physical_HiggsPropagatorMassError_Gev);
  printf("Physical Physical Higgs-Mass: (%1.3f +- %1.3f) GeV\n",Physical_PhysicalHiggsMass_GeV, Physical_PhysicalHiggsMassError_GeV);
  
  printf("Renormalized lambda: %1.3f +- %1.3f\n",Physical_RenormalizedLambda, Physical_RenormalizedLambdaError);
  printf("Renormalized Y: %1.3f +- %1.3f\n",Physical_RenormalizedY, Physical_RenormalizedYError);
  
  
  char* text = new char[10000];
  char* text2 = new char[10000];
  snprintf(text,10000,"Physical Cut-off:    (%1.3f $\\pm$ %1.3f) GeV\n\n",Physical_Lambda_GeV, Physical_LambdaError_GeV);
  snprintf(text2,10000,"%s Physical Top-Mass:         (%1.3f $\\pm$ %1.3f) GeV\n\n",text, Physical_PhysicalTopMass_Gev, Physical_PhysicalTopMassError_Gev);
  snprintf(text,10000,"%s Physical Prop. Higgs-Mass: (%1.3f $\\pm$ %1.3f) GeV\n\n",text2, Physical_HiggsPropagatorMass_Gev, Physical_HiggsPropagatorMassError_Gev);
  snprintf(text2,10000,"%s Physical Physical Higgs-Mass (%1.3f $\\pm$ %1.3f) GeV\n\n",text, Physical_PhysicalHiggsMass_GeV, Physical_PhysicalHiggsMassError_GeV);
  snprintf(text,10000,"%s Lattice vev (rescaled): %1.5f$\\pm$ %1.5f\n\n",text2, analyzerHiggs->LatticeResult_VEV, analyzerHiggs->LatticeResult_VEVsigma);
  snprintf(text2,10000,"%s Lattice Goldstone Z-Faktor (rescaled):     %1.5f $\\pm$ %1.5f\n\n",text, analyzerHiggs->LatticeResult_GoldstoneZFactor, analyzerHiggs->LatticeResult_GoldstoneZFactorSigma);
  snprintf(text,10000,"%s Lattice Top-Mass:               %1.5f $\\pm$ %1.5f\n\n",text2, analyzerTop->LatticeResult_PhysicalTopMass, analyzerTop->LatticeResult_PhysicalTopMassError);
  snprintf(text2,10000,"%s Lattice Higgs Propagator Mass:  %1.5f $\\pm$ %1.5f\n\n",text, analyzerHiggs->LatticeResult_HiggsPropagatorMass, analyzerHiggs->LatticeResult_HiggsPropagatorMassSigma);
  snprintf(text,10000,"%s Lattice Physical Higgs-Mass:          %1.5f $\\pm$ %1.5f\n\n",text2, analyzerHiggs->LatticeResult_PhysicalHiggsMass, analyzerHiggs->LatticeResult_PhysicalHiggsMassSigma);
  snprintf(text2,10000,"%s Lattice Size:  %dx%dx%dx%d\n\n",text, Parameter_L0,Parameter_L1,Parameter_L2,Parameter_L3);
  snprintf(text,10000,"%s Lattice $\\kappa_N$:          %1.5f\n\n",text2, Parameter_KappaN);
  snprintf(text2,10000,"%s Lattice $\\lambda_N$:  %1.5f\n\n",text, Parameter_LambdaN);
  snprintf(text,10000,"%s Lattice $y_N$:          %1.5f\n\n",text2, Parameter_YN);
  
  
  snprintf(text2,10000,"%s Rescaling Factor C:  %1.5f\n\n",text, HiggsFieldScale_C);
  snprintf(text,10000,"%s Continuum Notation $m_0^2$:          %1.5f\n\n",text2, Parameter_M0Sqr);
  snprintf(text2,10000,"%s Continuum Notation $y_0$:  %1.5f\n\n",text, Parameter_Y0);
  snprintf(text,10000,"%s Continuum Notation $\\lambda_0$:          %1.5f\n\n",text2, Parameter_Lambda0);
  
  snprintf(text2,10000,"%s Renormalized Lambda $\\lambda_R$:  %1.3f $\\pm$ %1.3f\n\n",text, Physical_RenormalizedLambda, Physical_RenormalizedLambdaError);
  snprintf(text,10000,"%s Renormalized Y $y_R$:          %1.3f $\\pm$ %1.3f\n\n",text2, Physical_RenormalizedY, Physical_RenormalizedYError);
  snprintf(text2,10000,"%s Evaluated configurations for Higgs-analysis:  %d\n\n",text, analyzerHiggs->getTotalN());
  snprintf(text,10000,"%s Evaluated configurations for Top-analysis:              %d\n\n",text2, analyzerTop->getTotalN());


  MassControlLog.addSection("Summary:");
  MassControlLog.addDirectText(text);  
  delete[] text;
  delete[] text2;
}


int main(int argc,char **argv) {
  LogLevel = 5;

  int ParaSelect = -1;
  if (argc>=2) {
    if (sscanf(argv[1],"%d",&ParaSelect)!=1)  ParaSelect = -1;
  }
  printf("Parameter-Selector is %d\n",ParaSelect);
  if ((abs(ParaSelect) >4) || (ParaSelect==0)) {
    printf("Invalid Parameter Selector!!!\n");
    exit(0);
  }
  
  
  Parameter_MaxRunTime = 16.0;

//  if (ParaSelect==-1) {
    FLAG_TexLog = true;
    FLAG_GnuplotFit = true;
//  } else {
//    FLAG_TexLog = false;
//    FLAG_GnuplotFit = false;  
//  }
  

  if (abs(ParaSelect) == 1) {
    printf("Chosen Parameter SET: 1\n");
    Parameter_L0 = 4;
    Parameter_L1 = 4;
    Parameter_L2 = 4;
    Parameter_L3 = 8;
    Parameter_Nf = 1;
//    Parameter_KappaN = 0.23754;
    Parameter_KappaN = 0.24001;
    Parameter_LambdaN = 1.0;  
    Parameter_YN = 0.711;
    Parameter_RHO = 1.0;
    Parameter_R = 0.5;
//    Parameter_PolyDegree = 16;
    Parameter_PolyDegree = 32;
    Parameter_PolyAlpha = 0.5;
    Parameter_filenameSuffix = "level8";	 
//    Parameter_filenameSuffix = "level6";	 
    Parameter_TopCorrTOL = 1E-7;
  }
  
  if (abs(ParaSelect) == 2) {
    printf("Chosen Parameter SET: 2\n");
    Parameter_L0 = 8;
    Parameter_L1 = 8;
    Parameter_L2 = 8;
    Parameter_L3 = 16;
    Parameter_Nf = 1;
//    Parameter_KappaN = 0.23760;
    Parameter_KappaN = 0.23894;
    Parameter_LambdaN = 1.0;  
    Parameter_YN = 0.711;
    Parameter_RHO = 1.0;
    Parameter_R = 0.5;
//    Parameter_PolyDegree = 24;
    Parameter_PolyDegree = 32;
    Parameter_PolyAlpha = 0.5;
//    Parameter_filenameSuffix = "level8";	 
    Parameter_filenameSuffix = "level8";	 
    Parameter_TopCorrTOL = 1E-7;
  }
  
  if (abs(ParaSelect) == 3) {
    printf("Chosen Parameter SET: 3\n");
    Parameter_L0 = 16;
    Parameter_L1 = 16;
    Parameter_L2 = 16;
    Parameter_L3 = 32;
    Parameter_Nf = 1;
    Parameter_KappaN = 0.24003;
    Parameter_LambdaN = 1.0;  
    Parameter_YN = 0.711;
    Parameter_RHO = 1.0;
    Parameter_R = 0.5;
    Parameter_PolyDegree = 32;
    Parameter_PolyAlpha = 0.5;
    Parameter_filenameSuffix = "level8";	 
    Parameter_TopCorrTOL = 1E-7;
  }

  if (abs(ParaSelect) == 4) {
    printf("Chosen Parameter SET: 4\n");
    Parameter_L0 = 20;
    Parameter_L1 = 20;
    Parameter_L2 = 20;
    Parameter_L3 = 32;
    Parameter_Nf = 1;
    Parameter_KappaN = 0.24006;
    Parameter_LambdaN = 1.0;  
    Parameter_YN = 0.711;
    Parameter_RHO = 1.0;
    Parameter_R = 0.5;
    Parameter_PolyDegree = 32;
    Parameter_PolyAlpha = 0.5;
    Parameter_filenameSuffix = "level8";	 
    Parameter_TopCorrTOL = 1E-7;
  }

  ini();

  int act;
/*  for (act=1; act<220; act+=1) {
    loadPhiconfiguration(act, true);
    
//multiplyHiggsFieldWithConst(-1.0, phiField);
    
//exactWeight = 1.0;
    analyzerHiggs->analyzeHiggsField((vector4D*) phiField, exactWeight);
  }
  printf("Total number of read configurations for Higgs-Analysis: %d \n",analyzerHiggs->getTotalN());
*/
  analyzerTop->loadTopTimeSliceCorrelator();
  for (act=1; act<460; act+=10) {
    loadPhiconfiguration(act, false);

    analyzerTop->analyzeHiggsField((vector4D*) phiField, exactWeight, act);
    
    if (timeOver()) {
      printf("Time Limit reached. ==> EXITING!!!\n");
      break;
    }
  }
  printf("Total number of read configurations for Top-Analysis: %d \n",analyzerTop->getTotalN());

  analyzerHiggs->plotSlottedGoldstonePropagator();
  analyzerHiggs->plotSlottedHiggsPropagator();
  analyzerHiggs->calcHiggsVEV();
  analyzerHiggs->calcHiggsTimeSliceCorrelator();
  analyzerHiggs->plotHiggsTimeSliceCorrelator();
  analyzerHiggs->calcGoldstoneTimeSliceCorrelator();
  analyzerHiggs->plotGoldstoneTimeSliceCorrelator();

  
  MassControlLog.clearPage();
  analyzerTop->calcTopTimeSliceCorrelator();
  analyzerTop->plotTopTimeSliceCorrelator(); 
  analyzerTop->saveTopTimeSliceCorrelator();


  analyzerHiggs->plotHiggsMasses();
  analyzerHiggs->plotGoldstone2ParticleMasses();
  analyzerHiggs->plotHiggsGoldstoneMasses();
  analyzerTop->plotTopMasses();


  printPhysicalAnalysis();


  MassControlLog.generateLatex();
}
