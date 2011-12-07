#include "AnalyzerObservableMultipleTimeScaleIntegration4.h"

AnalyzerObservableMultipleTimeScaleIntegration4::AnalyzerObservableMultipleTimeScaleIntegration4(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "MultipleTimeScaleIntegration4", "msint4") { 
  if (fermiOps->get1DSizeLargest() >= 32) {
    analyzeEveryXXXconf = 8;
  }
  MultiplePolynomFlag = false;
  if (SDReader->getSubPolynomCount()>0) MultiplePolynomFlag = true;
  ini(getAnalyzerResultsCount());
  pHMCProp = NULL;  
}


AnalyzerObservableMultipleTimeScaleIntegration4::~AnalyzerObservableMultipleTimeScaleIntegration4() {
  delete pHMCProp;
}


void AnalyzerObservableMultipleTimeScaleIntegration4::readFourierAccelerationData() {
  char* fileName1 = new char[600]; 
  char* fileName2 = new char[600]; 
  char* fileName3 = new char[600]; 
  char* fileName4 = new char[600]; 
  char* fileName5 = new char[600]; 
  char* fileName6 = new char[600]; 
  char* fileName7 = new char[600]; 
  char* fileName8 = new char[600]; 
  char* fileName9 = new char[600]; 
  char* fileName10 = new char[600]; 
  char* outputFileNameExtension = SDReader->getFileNameExtension();
  snprintf(fileName1,600,"%s/data/results/pHMC/FACC/PhiTotalForcePREC%s.dat",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName2,600,"%s/data/results/pHMC/FACC/PhiTotalForceGLOBAL%s.dat",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName3,600,"%s/data/results/pHMC/FACC/PhiHiggsForcePREC%s.dat",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName4,600,"%s/data/results/pHMC/FACC/PhiHiggsForceGLOBAL%s.dat",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName5,600,"%s/data/results/pHMC/FACC/PhiFermionForcePREC%s",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName6,600,"%s/data/results/pHMC/FACC/PhiFermionForceGLOBAL%s",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName7,600,"%s/data/results/pHMC/FACC/PhiChangeNoFACC%s.dat",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName8,600,"%s/data/results/pHMC/FACC/PhiChangeFACC%s.dat",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName9,600,"%s/data/results/pHMC/FACC/PhiTotalSphericalProjectedForcePREC%s.dat",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName10,600,"%s/data/results/pHMC/FACC/PhiTotalSphericalProjectedForceGLOBAL%s.dat",DataBaseDirectory,outputFileNameExtension);
  
  pHMCProp->loadPhiTotalForceFourierComponentsPREC(fileName1);
  pHMCProp->loadPhiTotalForceFourierComponentsGLOBAL(fileName2);
  pHMCProp->loadPhiHiggsForceFourierComponentsPREC(fileName3);
  pHMCProp->loadPhiHiggsForceFourierComponentsGLOBAL(fileName4);
  pHMCProp->loadPhiFermionForceFourierComponentsPREC(fileName5);
  pHMCProp->loadPhiFermionForceFourierComponentsGLOBAL(fileName6);
  pHMCProp->loadPhiChangeFourierComponentsNoFACC(fileName7);   
  pHMCProp->loadPhiChangeFourierComponentsFACC(fileName8);   
  pHMCProp->loadPhiTotalSphericalProjectedForceFourierComponentsPREC(fileName9);   
  pHMCProp->loadPhiTotalSphericalProjectedForceFourierComponentsGLOBAL(fileName10);     
  
  delete[] fileName1;
  delete[] fileName2;
  delete[] fileName3;
  delete[] fileName4;
  delete[] fileName5;
  delete[] fileName6;
  delete[] fileName7;
  delete[] fileName8;
  delete[] fileName9;
  delete[] fileName10;
  delete[] outputFileNameExtension;
}


int AnalyzerObservableMultipleTimeScaleIntegration4::howManyMatrixApplicationsForIntegrator(int* iter, int* Parameter_PolyDegree, int* Parameter_IntegratorType, int polyMode) {
  int howMany = 0;
  int startI = 4;
  int endI = 0;
  if (polyMode==1) startI = 0;
  if (polyMode==2) endI = 1;
  
  for (int I=startI; I>=endI; I--) {
    double add = Parameter_PolyDegree[I];
    
    for (int I2=I; I2>=endI; I2--) {
      add *= iter[I2];
      if (Parameter_IntegratorType[I2]==1) add *= 2;
      if (Parameter_IntegratorType[I2]==2) add *= 5;
    }

    howMany += roundToInt(add);
  }
  
  return howMany;
}


double AnalyzerObservableMultipleTimeScaleIntegration4::PerformIntegration(int outerIntSteps, int outerIntType, int HiggsIntSteps, int HiggsIntType, int &Nmmdag, double &epsi, int polyMode) {
  int I;
  int Parameter_SubPolyCount = SDReader->getSubPolynomCount();
  if ((polyMode>=2) && (!MultiplePolynomFlag)) {
    printf("ERROR in AnalyzerObservableMultipleTimeScaleIntegration4::PerformIntegration: invalid value for polyMode\n");
    exit(0);
  }
  int LevelCount = 2 + Parameter_SubPolyCount;
  int Parameter_Iterations[5];
  Parameter_Iterations[0] = SDReader->getPolynomIterations_P0();
  Parameter_Iterations[1] = SDReader->getPolynomIterations_P1();
  Parameter_Iterations[2] = SDReader->getPolynomIterations_P2();
  Parameter_Iterations[3] = SDReader->getPolynomIterations_P3();
  Parameter_Iterations[4] = SDReader->getPolynomIterations_P4();
  int Parameter_IntegratorType[5];
  Parameter_IntegratorType[0] = SDReader->getPolynomIntegrationType_P0();
  Parameter_IntegratorType[1] = SDReader->getPolynomIntegrationType_P1();
  Parameter_IntegratorType[2] = SDReader->getPolynomIntegrationType_P2();
  Parameter_IntegratorType[3] = SDReader->getPolynomIntegrationType_P3();
  Parameter_IntegratorType[4] = SDReader->getPolynomIntegrationType_P4();
  double Parameter_Epsilon = SDReader->getMolecularDynamicsEpsilon();
  int Parameter_PolyDegree[5];
  Parameter_PolyDegree[0] = SDReader->getPolynomDegree_P0();
  Parameter_PolyDegree[1] = SDReader->getPolynomDegree_P1();
  Parameter_PolyDegree[2] = SDReader->getPolynomDegree_P2();
  Parameter_PolyDegree[3] = SDReader->getPolynomDegree_P3();
  Parameter_PolyDegree[4] = SDReader->getPolynomDegree_P4();
  
  
/*printf("subpolcnt: %d\n", Parameter_SubPolyCount);
printf("iter0: %d\n", Parameter_Iterations[0]);
printf("iter1: %d\n", Parameter_Iterations[1]);
printf("iter2: %d\n", Parameter_Iterations[2]);
printf("iter3: %d\n", Parameter_Iterations[3]);
printf("iter4: %d\n", Parameter_Iterations[4]);
printf("intType0: %d\n", Parameter_IntegratorType[0]);
printf("intType1: %d\n", Parameter_IntegratorType[1]);
printf("intType2: %d\n", Parameter_IntegratorType[2]);
printf("intType3: %d\n", Parameter_IntegratorType[3]);
printf("intType4: %d\n", Parameter_IntegratorType[4]);
printf("eps: %f\n", Parameter_Epsilon);*/




  
  int* iterations = new int[LevelCount];
  for (I=0; I<1+Parameter_SubPolyCount; I++) {
    iterations[1+Parameter_SubPolyCount-I] = Parameter_Iterations[I];
  }
  iterations[0] = SDReader->getHiggsDynamicsIterations();
  
//printf("iterations Higgs: %d\n", iterations[0]);
  
  
  int* integrators = new int[LevelCount];
   for (I=0; I<1+Parameter_SubPolyCount; I++) {
    integrators[1+Parameter_SubPolyCount-I] = Parameter_IntegratorType[I];
  }
  integrators[0] = SDReader->getHiggsDynamicsIntegrationType();
//printf("intTypeHiggs: %d\n", integrators[0]);

  //Set explicit value for Higgs-Integration
  integrators[0] = HiggsIntType;
  iterations[0] = HiggsIntSteps;

  //Set explicit outer integrator
  if (polyMode==0) {
    Parameter_Epsilon *= Parameter_Iterations[0];
    Parameter_Epsilon /= outerIntSteps;
    epsi = Parameter_Epsilon;
    integrators[1+Parameter_SubPolyCount] = outerIntType;
    iterations[1+Parameter_SubPolyCount] = outerIntSteps;
    Parameter_IntegratorType[0] = outerIntType;
    Parameter_Iterations[0] = outerIntSteps;
  } else if (polyMode==1) {
    Parameter_Epsilon *= Parameter_Iterations[0];
    Parameter_Epsilon /= outerIntSteps;
    epsi = Parameter_Epsilon;
    integrators[1] = outerIntType;
    iterations[1] = outerIntSteps;
    Parameter_IntegratorType[0] = outerIntType;
    Parameter_Iterations[0] = outerIntSteps;
  } else {  
    Parameter_Epsilon *= Parameter_Iterations[0];
    Parameter_Epsilon /= outerIntSteps;
    epsi = Parameter_Epsilon;
    integrators[2+Parameter_SubPolyCount] = outerIntType;
    iterations[2+Parameter_SubPolyCount] = outerIntSteps;
    Parameter_IntegratorType[1] = outerIntType;
    Parameter_Iterations[1] = outerIntSteps;  
  }

  Nmmdag = howManyMatrixApplicationsForIntegrator(&(Parameter_Iterations[0]), &(Parameter_PolyDegree[0]), &(Parameter_IntegratorType[0]), polyMode);
  
  int* subPolNr = new int[LevelCount];
  if (polyMode==0) {
    for (I=2; I<LevelCount; I++) subPolNr[I] = 2*Parameter_SubPolyCount-2*I+3;
    subPolNr[0] = 0;
    subPolNr[1] = 2*Parameter_SubPolyCount;  
  } 
  if (polyMode==1) {
    subPolNr[0] = 0;
    subPolNr[1] = 0;  
    LevelCount = 2;
  } 
  if (polyMode==2) {
    for (I=2; I<LevelCount; I++) subPolNr[I] = 2*Parameter_SubPolyCount-2*I+3;
    subPolNr[0] = 0;
    subPolNr[1] = 2*Parameter_SubPolyCount;  
    LevelCount--;
  } 
    
  double* propTOL = NULL;
  double finalTOL = 0;

  pHMCProp->saveALLfields();
  pHMCProp->synchronizedChangeOfTuneMode(true);
  pHMCProp->multiTimeScalePropagation(LevelCount, iterations, Parameter_Epsilon, integrators, subPolNr, propTOL, finalTOL);
  double deltaS = pHMCProp->SafterProp - pHMCProp->SbeforeProp;
  pHMCProp->restoreALLfields(true);
  
  delete[] iterations;
  delete[] integrators;
  delete[] subPolNr;

//printf("deltaS: %1.15f\n",deltaS);

  return deltaS;
}


bool AnalyzerObservableMultipleTimeScaleIntegration4::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  if (pHMCProp==NULL) {
    double* polEps = new double[5];
    double* polLam = new double[5];
    int* polDeg = new int[5];
    
    polEps[0] = SDReader->getPolynomLowerBound_P0();
    polEps[1] = SDReader->getPolynomLowerBound_P1();
    polEps[2] = SDReader->getPolynomLowerBound_P2();
    polEps[3] = SDReader->getPolynomLowerBound_P3();
    polEps[4] = SDReader->getPolynomLowerBound_P4();

    polLam[0] = SDReader->getPolynomUpperBound();
    polLam[1] = SDReader->getPolynomUpperBound();
    polLam[2] = SDReader->getPolynomUpperBound();
    polLam[3] = SDReader->getPolynomUpperBound();
    polLam[4] = SDReader->getPolynomUpperBound();
    
    polDeg[0] = SDReader->getPolynomDegree_P0();
    polDeg[1] = SDReader->getPolynomDegree_P1();
    polDeg[2] = SDReader->getPolynomDegree_P2();
    polDeg[3] = SDReader->getPolynomDegree_P3();
    polDeg[4] = SDReader->getPolynomDegree_P4();
    
    pHMCProp = new pHMCPropagator(fermiOps, SDReader->getLambda(), SDReader->getKappa(), SDReader->getExternalCurrent(), 
                                  SDReader->getModelParameterC6(), SDReader->getModelParameterC8(), SDReader->getModelParameterC10(),
                                  SDReader->getModelParameterLambda6(), SDReader->getModelParameterLambda8(), SDReader->getModelParameterLambda10(),
                                  SDReader->getNf(), 1.0, 
                                  SDReader->getSphericalHiggsIntegrationMode(), SDReader->getZetaForHiggsIntegrationMode(),
                                  SDReader->getTheta(),
                                  SDReader->getSubPolynomCount(), polEps, polLam, polDeg, 0, NULL,
                                  SDReader->getPolynomDigits(), SDReader->getAlpha(), SDReader->getMaximumPolynomDegreePerNode(),10);

    bool quasiHermiteanMode = SDReader->getUseQHM();				  
    pHMCProp->setPhiForceFourierType(SDReader->getFourierAccelerationType(), SDReader->getFourieAccelerationParameter());  
    pHMCProp->setOmegaMassAdaptionMode(0);
    pHMCProp->resetExactMMdagInverseSQRTOmegaAction();
    pHMCProp->synchronizedChangeOfQuasiHermiteanMode(quasiHermiteanMode);
    pHMCProp->synchronizedChangeOfModelSelection(1);
    pHMCProp->getNodesReady();
    pHMCProp->activateForceStoring(true);
    readFourierAccelerationData();
    pHMCProp->calcPhiMomentumMasses(SDReader->getPolynomIterations_P0() * SDReader->getMolecularDynamicsEpsilon());    
  }
  int L0 = SDReader->getL0();
  int L1 = SDReader->getL1();
  int L2 = SDReader->getL2();
  int L3 = SDReader->getL3();


  double* phiField = phiFieldConf->getPhiFieldCopy(); 
  for (int I=0; I<L0*L1*L2*L3; I++) {
    pHMCProp->phiField[I][0] = phiField[4*I+0];
    pHMCProp->phiField[I][1] = phiField[4*I+1];
    pHMCProp->phiField[I][2] = phiField[4*I+2];
    pHMCProp->phiField[I][3] = phiField[4*I+3];    
  }

  for (int I=0; I<getAnalyzerResultsCount(); I++) {
    analyzerResults[I] = 0;
  }
  
  bool useP = SDReader->getUseP();
  bool useQ = SDReader->getUseQ();
  bool useR = SDReader->getUseR();
  double Pm = SDReader->getPPrecondParameterM();
  double Ps = SDReader->getPPrecondParameterS();
  double Qmu = SDReader->getQPrecondParameterMu();
  double Qbeta = SDReader->getQPrecondParameterBeta();
  double Rm = SDReader->getRPrecondParameterM();
  double Rf = SDReader->getRPrecondParameterF();

  fermiOps->setPreconditioner(useP, Pm, Ps);
  fermiOps->setQPreconditioner(useQ, Qmu, Qbeta);
  fermiOps->setRPreconditioner(useR, Rm, Rf);

  pHMCProp->sampleALLMomenta();
  pHMCProp->sampleOmegaFields();

  if (MultiplePolynomFlag) {
    int count = 0;
    for (int I=1; I<=40; I++) {
      if (I>10) I+=4;
      double epsi = NaN;
      int Nmmdag = 0;
      double deltaS = PerformIntegration(I, 0, 5, 2, Nmmdag, epsi, 1);
      analyzerResults[count*4 + 0] = I;
      analyzerResults[count*4 + 1] = deltaS;
      analyzerResults[count*4 + 2] = epsi;
      analyzerResults[count*4 + 3] = Nmmdag;
      count++;
    }
    
/*    for (int I=1; I<=40; I++) {
      if (I>10) I+=4;
      double epsi = NaN;
      int Nmmdag = 0;
      double deltaS = PerformIntegration(I, 1, 5, 2, Nmmdag, epsi, 1);
      analyzerResults[count*4 + 0] = I;
      analyzerResults[count*4 + 1] = deltaS;
      analyzerResults[count*4 + 2] = epsi;
      analyzerResults[count*4 + 3] = Nmmdag;
      count++;
    }
    
    for (int I=1; I<=8; I++) {
      double epsi = NaN;
      int Nmmdag = 0;
      double deltaS = PerformIntegration(I, 2, 5, 2, Nmmdag, epsi, 1);
      analyzerResults[count*4 + 0] = I;
      analyzerResults[count*4 + 1] = deltaS;
      analyzerResults[count*4 + 2] = epsi;
      analyzerResults[count*4 + 3] = Nmmdag;
      count++;
    }*/
    
  } else {
    int count = 0;
    for (int I=1; I<=40; I++) {
      if (I>10) I+=4;
      double epsi = NaN;
      int Nmmdag = 0;
      double deltaS = PerformIntegration(I, 0, 1, 0, Nmmdag, epsi, 0);
      analyzerResults[count*4 + 0] = I;
      analyzerResults[count*4 + 1] = deltaS;
      analyzerResults[count*4 + 2] = epsi;
      analyzerResults[count*4 + 3] = Nmmdag;
      count++;
    }

    for (int I=1; I<=40; I++) {
      if (I>10) I+=4;
      double epsi = NaN;
      int Nmmdag = 0;
      double deltaS = PerformIntegration(I, 1, 1, 0, Nmmdag, epsi, 0);
      analyzerResults[count*4 + 0] = I;
      analyzerResults[count*4 + 1] = deltaS;
      analyzerResults[count*4 + 2] = epsi;
      analyzerResults[count*4 + 3] = Nmmdag;
      count++;
    }

    for (int I=1; I<=8; I++) {
      double epsi = NaN;
      int Nmmdag = 0;
      double deltaS = PerformIntegration(I, 2, 1, 0, Nmmdag, epsi, 0);
      analyzerResults[count*4 + 0] = I;
      analyzerResults[count*4 + 1] = deltaS;
      analyzerResults[count*4 + 2] = epsi;
      analyzerResults[count*4 + 3] = Nmmdag;
      count++;
    }
    

    for (int I=1; I<=40; I++) {
      if (I>10) I+=4;
      double epsi = NaN;
      int Nmmdag = 0;
      double deltaS = PerformIntegration(I, 0, 5, 2, Nmmdag, epsi, 0);
      analyzerResults[count*4 + 0] = I;
      analyzerResults[count*4 + 1] = deltaS;
      analyzerResults[count*4 + 2] = epsi;
      analyzerResults[count*4 + 3] = Nmmdag;
      count++;
    }

    for (int I=1; I<=40; I++) {
      if (I>10) I+=4;
      double epsi = NaN;
      int Nmmdag = 0;
      double deltaS = PerformIntegration(I, 1, 5, 2, Nmmdag, epsi, 0);
      analyzerResults[count*4 + 0] = I;
      analyzerResults[count*4 + 1] = deltaS;
      analyzerResults[count*4 + 2] = epsi;
      analyzerResults[count*4 + 3] = Nmmdag;
      count++;
    }

    for (int I=1; I<=8; I++) {
      double epsi = NaN;
      int Nmmdag = 0;
      double deltaS = PerformIntegration(I, 2, 5, 2, Nmmdag, epsi, 0);
      analyzerResults[count*4 + 0] = I;
      analyzerResults[count*4 + 1] = deltaS;
      analyzerResults[count*4 + 2] = epsi;
      analyzerResults[count*4 + 3] = Nmmdag;
      count++;
    }
  }

  return true;
}


int AnalyzerObservableMultipleTimeScaleIntegration4::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservableMultipleTimeScaleIntegration4::getAnalyzerResultsCount() {
  if (MultiplePolynomFlag) {
    return 64;  
  } else {
    return 320;
  }
}
