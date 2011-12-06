#include "Propagator.h"

Propagator::Propagator(FermionMatrixOperations* fOps, double lam, double kap, double current, double c6, double c8, double c10, double lam6, double lam8, double lam10, int nf, double gam, bool sphMode, double sphZeta) {
  fermiOps = fOps;  
  neighbours = NULL;
  phiField = NULL;
  momentaField = NULL;
  forces = NULL;
  phiFieldOLD = NULL;
  momentaFieldOLD = NULL;
  Sold = NaN;
  SoldOLD = NaN;
  Sact = NaN;
  SbeforeProp = NaN;
  SafterProp = NaN;
  forceCount = 0;
  kappa = NaN;
  explicitCurrent = NaN;  
  ModelExtension_c6 = NaN;
  ModelExtension_c8 = NaN;
  ModelExtension_c10 = NaN;
  ModelExtension_lambda6 = NaN;
  ModelExtension_lambda8 = NaN;
  ModelExtension_lambda10 = NaN;  
  lambda = NaN;
  Nf = 0;
  gamma = NaN;
  SphericalMode = false;
  SphericalZeta = 0;  
  phiForceInternal = NULL;
  vectorInterim = NULL; 
  phiForceInternalPlanForward = NULL;
  phiForceInternalPlanBackward = NULL;  
  vectorInterimPlanForward = NULL;
  vectorInterimPlanBackward = NULL;
  momentaToInterimPlanForward = NULL;
  interimToMomentaPlanBackward = NULL;  
  internalToInterimPlanForward = NULL; 
  phiFieldToInterimPlanForward = NULL;
  phiForceFourierType = 0;
  phiForceFourierPara = 0;
  momentumMassesDetermined = false;
  tuneMode = false;
  phiTotalSphericalProjectedForceAnalysisPREC = NULL;
  phiTotalSphericalProjectedForceAnalysisGLOBAL = NULL;      
  phiTotalForceAnalysisPREC = NULL;
  phiTotalForceAnalysisGLOBAL = NULL;
  phiHiggsForceAnalysisPREC = NULL;
  phiHiggsForceAnalysisGLOBAL = NULL;
  phiFermionForceAnalysisPREC = NULL;
  phiFermionForceAnalysisGLOBAL = NULL;
  phiChangeAnalysisNoFACC = NULL;
  phiChangeAnalysisFACC = NULL;
  MomentaMasses = NULL;
  MomentumSpaceActionCoefficients = NULL;
}


void Propagator::ini(FermionMatrixOperations* fOps, double lam, double kap, double current, double c6, double c8, double c10, double lam6, double lam8, double lam10, int nf, double gam, bool sphMode, double sphZeta) {
  if (LogLevel>2) printf("Initializing Propagator with lambda = %1.3f, kappa = %1.3f, current = %1.3e, c6 = %1.3e, c8 = %1.3e, c10 = %1.3e, lam6 = %1.3e, lam8 = %1.3e, lam10 = %1.3e, Nf = %d, gamma = %1.3f, SphericalMode = %d, and SphericalZeta = %1.3e\n", lam, kap, current, c6, c8, c10, lam6, lam8, lam10, nf, gam, sphMode, sphZeta);  
  setKappa(kap);
  setCurrent(current);
  setModelExtensionParameterC6(c6);
  setModelExtensionParameterC8(c8);
  setModelExtensionParameterC10(c10);
  setModelExtensionParameterLambda6(lam6);
  setModelExtensionParameterLambda8(lam8);
  setModelExtensionParameterLambda10(lam10);
  setLambda(lam);
  setSphericalMode(sphMode, sphZeta);  
  setGamma(gam);  
  setNf(nf);  
  iniFields();
  iniAdditionalFields();
}


void Propagator::desiniFields() {
  delete[] neighbours;
  neighbours = NULL;
  Complex* c1 = (Complex*) phiField;
  Complex* c2 = (Complex*) phiFieldOLD;
  Complex* c3 = (Complex*) momentaField;
  Complex* c4 = (Complex*) momentaFieldOLD;
  Complex* c5 = (Complex*) phiForceInternal;
  Complex* c6 = (Complex*) vectorInterim;
  Complex* c7 = (Complex*) MomentaMasses; 
  Complex* c8 = (Complex*) MomentumSpaceActionCoefficients;
  fermiOps->destroyFermionVector(c1);
  fermiOps->destroyFermionVector(c2);
  fermiOps->destroyFermionVector(c3);
  fermiOps->destroyFermionVector(c4);
  fermiOps->destroyFermionVector(c5);
  fermiOps->destroyFermionVector(c6);
  fermiOps->destroyFermionVector(c7);
  fermiOps->destroyFermionVector(c8);  
  phiField = NULL;
  phiFieldOLD = NULL;
  momentaField = NULL;
  momentaFieldOLD = NULL;
  phiForceInternal = NULL;
  vectorInterim = NULL;
  MomentaMasses = NULL;
  MomentumSpaceActionCoefficients = NULL;
  
  fftw_destroy_plan(phiForceInternalPlanForward);
  fftw_destroy_plan(phiForceInternalPlanBackward);
  fftw_destroy_plan(vectorInterimPlanForward);
  fftw_destroy_plan(vectorInterimPlanBackward);
  fftw_destroy_plan(momentaToInterimPlanForward);
  fftw_destroy_plan(interimToMomentaPlanBackward);
  fftw_destroy_plan(internalToInterimPlanForward);
  fftw_destroy_plan(phiFieldToInterimPlanForward);
    
  delete phiTotalSphericalProjectedForceAnalysisPREC;
  delete phiTotalSphericalProjectedForceAnalysisGLOBAL;  
  delete phiTotalForceAnalysisPREC;
  delete phiTotalForceAnalysisGLOBAL;
  delete phiHiggsForceAnalysisPREC;
  delete phiHiggsForceAnalysisGLOBAL;
  delete phiChangeAnalysisNoFACC;
  delete phiChangeAnalysisFACC;
  
  phiTotalSphericalProjectedForceAnalysisPREC = NULL;
  phiTotalSphericalProjectedForceAnalysisGLOBAL = NULL;    
  phiTotalForceAnalysisPREC = NULL;
  phiTotalForceAnalysisGLOBAL = NULL;
  phiHiggsForceAnalysisPREC = NULL;
  phiHiggsForceAnalysisGLOBAL = NULL;
  phiChangeAnalysisNoFACC = NULL;
  phiChangeAnalysisFACC = NULL;

  int I;
  for (I=0; I<getFermionForceSubCategoryCount(); I++) {
    delete phiFermionForceAnalysisPREC[I];
    delete phiFermionForceAnalysisGLOBAL[I];
  }
  delete[] phiFermionForceAnalysisPREC;
  delete[] phiFermionForceAnalysisGLOBAL;
  phiFermionForceAnalysisPREC = NULL;
  phiFermionForceAnalysisGLOBAL = NULL;

  desiniAdditionalFields();
}


void Propagator::desiniForceCalculators() {
  int I;
  for (I=0; I<forceCount; I++) {
    delete forces[I];
    forces[I] = NULL;
  }
  delete[] forces;
  forces = NULL;
  forceCount = 0;
}


void Propagator::desini() {
  if (LogLevel>0) printf("Desinitializing Propagator...");
  desiniFields();
  desiniForceCalculators();
  
  if (LogLevel>0) printf("sucessfully.\n");      
}

  
void Propagator::killSlaves() {
  if (LogLevel>0) printf("Killing all slaves...\n");
  threadedExecute(Propagator_EXIT, 0);   
  if (LogLevel>0) printf("...successful!\n");
}
  
  
void Propagator::setKappa(double kap) {
  kappa = kap;
}


void Propagator::setCurrent(double current) {
  explicitCurrent = current;
  if (LogLevel > 1) printf("Setting explicit current in Propagator to %1.3e\n", explicitCurrent);
}


void Propagator::setModelExtensionParameterC6(double c6) {
  ModelExtension_c6 = c6;
  if (LogLevel > 1) printf("Setting ModelExtensionParameterC6 in Propagator to %1.3e\n", ModelExtension_c6);
}


void Propagator::setModelExtensionParameterC8(double c8) {
  ModelExtension_c8 = c8;
  if (LogLevel > 1) printf("Setting ModelExtensionParameterC8 in Propagator to %1.3e\n", ModelExtension_c8);
}


void Propagator::setModelExtensionParameterC10(double c10) {
  ModelExtension_c10 = c10;
  if (LogLevel > 1) printf("Setting ModelExtensionParameterC10 in Propagator to %1.3e\n", ModelExtension_c10);
}


void Propagator::setModelExtensionParameterLambda6(double lam6) {
  ModelExtension_lambda6 = lam6;
  if (LogLevel > 1) printf("Setting ModelExtensionParameterLambda6 in Propagator to %1.3e\n", ModelExtension_lambda6);
}


void Propagator::setModelExtensionParameterLambda8(double lam8) {
  ModelExtension_lambda8 = lam8;
  if (LogLevel > 1) printf("Setting ModelExtensionParameterLambda8 in Propagator to %1.3e\n", ModelExtension_lambda8);
}


void Propagator::setModelExtensionParameterLambda10(double lam10) {
  ModelExtension_lambda10 = lam10;
  if (LogLevel > 1) printf("Setting ModelExtensionParameterLambda10 in Propagator to %1.3e\n", ModelExtension_lambda10);
}


void Propagator::setGamma(double gam) {
  gamma = gam;
}


void Propagator::setLambda(double lam) {
  lambda = lam;
}


void Propagator::setSphericalMode(bool sphMode, double sphZeta) {
  SphericalMode = sphMode;
  SphericalZeta = abs(sphZeta);  
  if (LogLevel>2) printf("Spherical-Mode is set to %d with Zeta = %1.3e\n",SphericalMode,SphericalZeta);
}


void Propagator::checkValidityOfSetting() {
  if (isNaN(lambda)) {
    if  ((!SphericalMode) || (SphericalZeta!=0)) {
      printf("Incompatible setting: lambda=%1.3e, SphericalMode=%d, and SphericalZeta=%1.3e\n",lambda,SphericalMode,SphericalZeta);
      exit(0);
    }
  } else {
    if  ((SphericalMode) && (SphericalZeta<=0)) {
      printf("Incompatible setting: lambda=%1.3e, SphericalMode=%d, and SphericalZeta=%1.3e\n",lambda,SphericalMode,SphericalZeta);
      exit(0);
    }
  }
}


void Propagator::samplePhiMomentumField() {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  if (LogLevel>4) printf("Sampling phi momenta directly.\n");
  checkValidityOfSetting();
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  
  fermiOps->fillGaussRandomVector((Complex*)momentaField,2*L0*L1*L2*L3);
  double fac = 1/sqrt(gamma);
  if (phiForceFourierType>0) {
    if (SphericalMode) {    
      for (int I=0; I<L0*L1*L2*L3; I++) {
        double f = fac*sqrt(MomentaMasses[I]/(L0*L1*L2*L3));
	f /= sqrt(1.0 + SphericalZeta);
        vectorInterim[I][0] = f*momentaField[I][0];
        vectorInterim[I][1] = f*momentaField[I][1];
        vectorInterim[I][2] = f*momentaField[I][2];
        vectorInterim[I][3] = f*momentaField[I][3];
      }
      fftw_execute(interimToMomentaPlanBackward);  
    } else {
      for (int I=0; I<L0*L1*L2*L3; I++) {
        double f = fac*sqrt(MomentaMasses[I]);
        momentaField[I][0] *= f;
        momentaField[I][1] *= f;
        momentaField[I][2] *= f;
        momentaField[I][3] *= f;
      }
    }
  } else {
    if (SphericalMode) {      
      double f = fac / sqrt(1.0 + SphericalZeta);    
      for (int I=0; I<L0*L1*L2*L3; I++) {
        momentaField[I][0] *= f;
        momentaField[I][1] *= f;
        momentaField[I][2] *= f;
        momentaField[I][3] *= f;	
      }  
    } else {
      for (int I=0; I<L0*L1*L2*L3; I++) {
        momentaField[I][0] *= fac;
        momentaField[I][1] *= fac;
        momentaField[I][2] *= fac;
        momentaField[I][3] *= fac;
      }  
    }
  }
  Sold = NaN;
  Sact = NaN;  
  SbeforeProp = NaN;
  SafterProp = NaN;
  addPerformanceProfilingItem("Propagator::samplePhiMomentumField", performanceProfilerStartCycle, 0);  
}


void Propagator::multiplyPhiMomenta(double fac) {
  if (LogLevel>4) printf("Multiplying phi momenta with factor %f.\n",fac);
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  
  int I;
  for (I=0; I<L0*L1*L2*L3; I++) {
    momentaField[I][0] *= fac;
    momentaField[I][1] *= fac;
    momentaField[I][2] *= fac;
    momentaField[I][3] *= fac;
  }
}


void Propagator::iniFields() {
  int x0,x1,x2,x3,I;
  int count = 0;
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  
  if (LogLevel>2) printf("Initializing fields (%d pseudo fermion fields)...", forceCount);
  neighbours = new NeighboursType[L0*L1*L2*L3];
  
  for (x0=0; x0<L0; x0++) {
    for (x1=0; x1<L1; x1++) {
      for (x2=0; x2<L2; x2++) {
        for (x3=0; x3<L3; x3++) {
	  neighbours[count].neighbourIndex[0] = LPOS(x0+1,x1,x2,x3);
	  neighbours[count].neighbourIndex[1] = LPOS(x0-1,x1,x2,x3);
	  neighbours[count].neighbourIndex[2] = LPOS(x0,x1+1,x2,x3);
	  neighbours[count].neighbourIndex[3] = LPOS(x0,x1-1,x2,x3);
	  neighbours[count].neighbourIndex[4] = LPOS(x0,x1,x2+1,x3);
	  neighbours[count].neighbourIndex[5] = LPOS(x0,x1,x2-1,x3);
	  neighbours[count].neighbourIndex[6] = LPOS(x0,x1,x2,x3+1);
	  neighbours[count].neighbourIndex[7] = LPOS(x0,x1,x2,x3-1);
	  count++;
	}
      }
    }
  }
  
  phiField = (vector4D*) fermiOps->createFermionVector(2);
  phiFieldOLD = (vector4D*) fermiOps->createFermionVector(2);
  if (!(lambda==lambda)) {
    //Lambda == Infinity (i.e. NaN)
    for (I=0; I<L0*L1*L2*L3; I++) {
      phiField[I][0] = sqrt(Nf);
      phiField[I][1] = 0;
      phiField[I][2] = 0;
      phiField[I][3] = 0;
    }
  } else {
    //Lambda < Infinity
    double sqrtNf4 = sqrt(0.25*Nf);
    for (I=0; I<L0*L1*L2*L3; I++) {
      phiField[I][0] = sqrtNf4;
      phiField[I][1] = sqrtNf4;
      phiField[I][2] = sqrtNf4;
      phiField[I][3] = sqrtNf4;
    }
  }
  momentaField = (vector4D*) fermiOps->createFermionVector(2);
  momentaFieldOLD = (vector4D*) fermiOps->createFermionVector(2);
  
  iniPhiForceFourierData();    
  
  if (LogLevel>0) printf("sucessfully.\n");    
}


void Propagator::iniPhiForceFourierData() {
  int I;
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
    
  phiTotalSphericalProjectedForceAnalysisPREC = new LatticeMomentumBins(L0,L1,L2,L3);
  phiTotalSphericalProjectedForceAnalysisPREC->clearData();
  phiTotalSphericalProjectedForceAnalysisGLOBAL = new LatticeMomentumBins(L0,L1,L2,L3);
  phiTotalSphericalProjectedForceAnalysisGLOBAL->clearData();
  phiTotalForceAnalysisPREC = new LatticeMomentumBins(L0,L1,L2,L3);
  phiTotalForceAnalysisPREC->clearData();
  phiTotalForceAnalysisGLOBAL = new LatticeMomentumBins(L0,L1,L2,L3);
  phiTotalForceAnalysisGLOBAL->clearData();  
  phiHiggsForceAnalysisPREC = new LatticeMomentumBins(L0,L1,L2,L3);
  phiHiggsForceAnalysisPREC->clearData();
  phiHiggsForceAnalysisGLOBAL = new LatticeMomentumBins(L0,L1,L2,L3);
  phiHiggsForceAnalysisGLOBAL->clearData();  
  phiChangeAnalysisNoFACC = new LatticeMomentumBins(L0,L1,L2,L3);
  phiChangeAnalysisNoFACC->clearData();  
  phiChangeAnalysisFACC = new LatticeMomentumBins(L0,L1,L2,L3);
  phiChangeAnalysisFACC->clearData();  
  
  phiFermionForceAnalysisPREC = new LatticeMomentumBins*[getFermionForceSubCategoryCount()];
  phiFermionForceAnalysisGLOBAL = new LatticeMomentumBins*[getFermionForceSubCategoryCount()];
  for (I=0; I<getFermionForceSubCategoryCount(); I++) {
    phiFermionForceAnalysisPREC[I] = new LatticeMomentumBins(L0,L1,L2,L3);
    phiFermionForceAnalysisPREC[I]->clearData();
    phiFermionForceAnalysisGLOBAL[I] = new LatticeMomentumBins(L0,L1,L2,L3);
    phiFermionForceAnalysisGLOBAL[I]->clearData();  
  }
  
  phiForceInternal = (vector4D*) fermiOps->createFermionVector(2);
  vectorInterim = (vector4D*) fermiOps->createFermionVector(2);  
  MomentaMasses = (double*) fermiOps->createFermionVector(1);  
  MomentumSpaceActionCoefficients = (double*) fermiOps->createFermionVector(1);    
  
  for (I=0; I<L0*L1*L2*L3; I++) {
    MomentaMasses[I] = 1.0;  
    MomentaMasses[I+(L0*L1*L2*L3)] = 0.0;   
    double pSqr = phiTotalForceAnalysisGLOBAL->getLatMomSqrFromIndex(I);
    MomentumSpaceActionCoefficients[I]  = ModelExtension_c6 * pSqr*pSqr;
    MomentumSpaceActionCoefficients[I] += ModelExtension_c8 * pSqr*pSqr*pSqr;
    MomentumSpaceActionCoefficients[I] += ModelExtension_c10* pSqr*pSqr*pSqr*pSqr;    
    MomentumSpaceActionCoefficients[I+(L0*L1*L2*L3)] = NaN;
  }
  resetPhiTotalForceFourierComponentsPREC();
  resetPhiTotalForceFourierComponentsGLOBAL();
  resetPhiHiggsForceFourierComponentsPREC();
  resetPhiHiggsForceFourierComponentsGLOBAL();
  resetPhiFermionForceFourierComponentsPREC();
  resetPhiFermionForceFourierComponentsGLOBAL();
  resetPhiChangeFourierComponentsNoFACC();
  resetPhiChangeFourierComponentsFACC();
  
  int* n = new int[4];
  n[0] = L0;
  n[1] = L1;
  n[2] = L2;
  n[3] = L3;
  
  int rank = 4;
  int howmany = 2;
  int* inembed = new int[4];
  inembed[0] = L0;
  inembed[1] = L1;
  inembed[2] = L2;
  inembed[3] = L3;  
  int istride = 2;
  int idist = 1;
  
  int* onembed = new int[4];
  onembed[0] = L0;
  onembed[1] = L1;
  onembed[2] = L2;
  onembed[3] = L3;    
  int ostride = 2;
  int odist = 1;

  phiForceInternalPlanForward = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*) phiForceInternal, inembed,
                                            istride, idist,
	   			            (fftw_complex*)phiForceInternal, onembed, ostride, odist,
	  				    FFTW_FORWARD, FFTW_MEASURE);
				
				
  phiForceInternalPlanBackward = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*) phiForceInternal, inembed,
                                            istride, idist,
	   			            (fftw_complex*)phiForceInternal, onembed, ostride, odist,
	  				    FFTW_BACKWARD, FFTW_MEASURE);
				
					    
  vectorInterimPlanForward = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*) vectorInterim, inembed,
                                            istride, idist,
	   			            (fftw_complex*)vectorInterim, onembed, ostride, odist,
	  				    FFTW_FORWARD, FFTW_MEASURE);					   						

  vectorInterimPlanBackward = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*) vectorInterim, inembed,
                                            istride, idist,
	   			            (fftw_complex*)vectorInterim, onembed, ostride, odist,
	  				    FFTW_BACKWARD, FFTW_MEASURE);					   						


  momentaToInterimPlanForward = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*) momentaField, inembed,
                                            istride, idist,
	   			            (fftw_complex*)vectorInterim, onembed, ostride, odist,
	  				    FFTW_FORWARD, FFTW_MEASURE);					   						

  interimToMomentaPlanBackward = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*) vectorInterim, inembed,
                                            istride, idist,
	   			            (fftw_complex*)momentaField, onembed, ostride, odist,
	  				    FFTW_BACKWARD, FFTW_MEASURE);					   						
					    
  internalToInterimPlanForward  = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*) phiForceInternal, inembed,
                                            istride, idist,
	   			            (fftw_complex*)vectorInterim, onembed, ostride, odist,
	  				    FFTW_FORWARD, FFTW_MEASURE);					   						

  phiFieldToInterimPlanForward = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*) phiField, inembed,
                                            istride, idist,
	   			            (fftw_complex*)vectorInterim, onembed, ostride, odist,
	  				    FFTW_FORWARD, FFTW_MEASURE);
    						
  delete[] n;
  delete[] inembed;
  delete[] onembed;
}


void Propagator::resetPhiTotalSphericalProjectedForceFourierComponentsPREC() {
  phiTotalSphericalProjectedForceAnalysisPREC->clearData();
}


void Propagator::resetPhiTotalSphericalProjectedForceFourierComponentsGLOBAL() {
  phiTotalSphericalProjectedForceAnalysisGLOBAL->clearData();
}


void Propagator::resetPhiTotalForceFourierComponentsPREC() {
  phiTotalForceAnalysisPREC->clearData();
}


void Propagator::resetPhiTotalForceFourierComponentsGLOBAL() {
  phiTotalForceAnalysisGLOBAL->clearData();
}


void Propagator::resetPhiFermionForceFourierComponentsPREC() {
  int I;
  for (I=0; I<getFermionForceSubCategoryCount(); I++) {
    phiFermionForceAnalysisPREC[I]->clearData();
  }
}


void Propagator::resetPhiFermionForceFourierComponentsGLOBAL() {
  int I;
  for (I=0; I<getFermionForceSubCategoryCount(); I++) {
    phiFermionForceAnalysisGLOBAL[I]->clearData();
  }
}


void Propagator::resetPhiHiggsForceFourierComponentsPREC() {
  phiHiggsForceAnalysisPREC->clearData();
}


void Propagator::resetPhiHiggsForceFourierComponentsGLOBAL() {
  phiHiggsForceAnalysisGLOBAL->clearData();
}


void Propagator::resetPhiChangeFourierComponentsNoFACC() {
  phiChangeAnalysisNoFACC->clearData();
}


void Propagator::resetPhiChangeFourierComponentsFACC() {
  phiChangeAnalysisFACC->clearData();
}


double Propagator::optimalMomentumMass(double avgForce, double traLength, double latMomSqr, int i0, int i1, int i2, int i3) {
  if (phiForceFourierType == 1) {
    double normFac = 0.25*(avgForce)*(avgForce);
    return normFac;
  } else if (phiForceFourierType == 2) {
    double mF = phiForceFourierPara;
    double normFac = 1.0;
    double boostFac = (latMomSqr + sqr(mF)) / (16 + sqr(mF));
    return normFac*boostFac;
  } else if (phiForceFourierType == 3) {
    double mF = phiForceFourierPara;
    double normFac = 0.25*(avgForce)*(avgForce);
    double boostFac = (latMomSqr + sqr(mF)) / (16 + sqr(mF));
    return normFac*boostFac;
  } else if (phiForceFourierType == 4) {
    double mF = phiForceFourierPara;
    double normFac = 0.25*(avgForce)*(avgForce);
    double boostFac = 1.0;
    if ((i0==0) && (i1==0) && (i2==0)) boostFac = (latMomSqr + sqr(mF)) / (16 + sqr(mF));
    if ((i0==0) && (i1==0) && (i3==0)) boostFac = (latMomSqr + sqr(mF)) / (16 + sqr(mF));
    if ((i0==0) && (i2==0) && (i3==0)) boostFac = (latMomSqr + sqr(mF)) / (16 + sqr(mF));
    if ((i1==0) && (i2==0) && (i3==0)) boostFac = (latMomSqr + sqr(mF)) / (16 + sqr(mF));
    return normFac*boostFac;
  } else if (phiForceFourierType == 5) {
    double maxTraFac = ((int)phiForceFourierPara);
    double additionalBoost = 4.0*(phiForceFourierPara - maxTraFac);
    double boostFac = 1.0;
    if (latMomSqr<1.0) boostFac = 1.0 + ((1.0-latMomSqr)*additionalBoost);
    double normFac = 0.25*(avgForce)*(avgForce);
    double mass = normFac/(boostFac*boostFac);
    if ((1.0/sqrt(mass))>maxTraFac) mass = 1.0/(maxTraFac*maxTraFac);
    return mass;
  } else {
    printf("Unknown FACC - Type!!!\n");
    exit(0);
  }
}


void Propagator::calcPhiMomentumMasses(double traLength) {
  int i0,i1,i2,i3;

  if (LogLevel>2) printf("\nSetting Momenta Masses!!\n");

  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();

  if (SphericalMode) {
    if (!phiTotalSphericalProjectedForceAnalysisPREC->containsData()) {
      printf("ERROR: No Data available for Setting Momenta Masses (yN>0)!!!\n");
      if (phiForceFourierType == 2) {
        printf("...but ok, due to phiForceFourierType == 2!!!\n");
      } else {
        exit(0);
      }
    } else {
      phiTotalSphericalProjectedForceAnalysisPREC->getAverageVectorInflated((double*)vectorInterim);  
    }
  } else {
    if (!phiTotalForceAnalysisPREC->containsData()) {
      printf("ERROR: No Data available for Setting Momenta Masses (yN>0)!!!\n");
      if (phiForceFourierType == 2) {
        printf("...but ok, due to phiForceFourierType == 2!!!\n");
      } else {
        exit(0);
      }
    } else {
      phiTotalForceAnalysisPREC->getAverageVectorInflated((double*)vectorInterim);  
    }
  }

  int count = 0;
  for (i0=0; i0<L0; i0++) {
    for (i1=0; i1<L1; i1++) {
      for (i2=0; i2<L2; i2++) {
        for (i3=0; i3<L3; i3++) {
          double avgForce = ((double*)vectorInterim)[count];
          double latMomsqr = phiTotalForceAnalysisPREC->getLatMomSqrFromIndex(count);
          MomentaMasses[count] = optimalMomentumMass(avgForce, traLength, latMomsqr, i0, i1, i2, i3);
	  count++;
	}
      }
    }
  }
  momentumMassesDetermined = true;
  threadedExecute(Propagator_SETFOURIERACCELERATIONMASSES,0);
}


void Propagator::setPhiForceFourierType(int type, double para) {
  if (type != phiForceFourierType) changeOfFACCtype();
  if ((type != phiForceFourierType) || (para != phiForceFourierPara)) {
    if (phiForceFourierType==0) {
      if (LogLevel > 1) printf("Phi Force Fourier Mode activated, type=%d and para=%f.\n",type,para);
    } else if (type == 0) {
      if (LogLevel > 1) printf("Phi Force Fourier Mode deactivated, type=%d and para=%f.\n",type,para);    
    } else {
      if (LogLevel > 1) printf("Phi Force Fourier Mode changed, type=%d and para=%f.\n",type,para);        
    }
  }
  phiForceFourierType = type;
  phiForceFourierPara = para;
  
  threadedExecute(Propagator_SETFOURIERACCELERATIONTYPE,0);
} 


void Propagator::calcPhiForceFourierComponents(int forceType, int fermionForceSubCat, bool takeOverFFTresult) {
  calcPhiForceFourierComponents((double*) phiForceInternal, forceType, fermionForceSubCat, takeOverFFTresult);
}


void Propagator::calcPhiForceFourierComponents(double* dSdPhi, int forceType, int fermionForceSubCat, bool takeOverFFTresult) {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  if (LogLevel>4) printf("Computing Phi - Force - Fourier - Components...\n");

  if ((fermionForceSubCat!=0) && (forceType!=2)) {
    printf("ERROR: Force Sub-Cat (%d) not compatible with force-Selection %d!!!\n",fermionForceSubCat,forceType);
    exit(0);
  }

  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();

  if ((!tuneMode) && (SphericalMode) && ((forceType == 0) || ((fermiOps->getYN()==0) && (forceType==1)))) {
    double fac = 1.0 /sqrt(L0*L1*L2*L3);
    int count = 0;
    for (int I=0; I<L0*L1*L2*L3; I++) {      
      double ph0 = phiField[I][0];
      double ph1 = phiField[I][1];
      double ph2 = phiField[I][2];
      double ph3 = phiField[I][3];
      double dS0 = dSdPhi[count+0];
      double dS1 = dSdPhi[count+1];
      double dS2 = dSdPhi[count+2];
      double dS3 = dSdPhi[count+3];      
      double dp1 = -dS0*ph1 + dS1*ph0 - dS2*ph3 + dS3*ph2;     
      double dp2 = -dS0*ph2 + dS2*ph0 + dS1*ph3 - dS3*ph1;     
      double dp3 = -dS0*ph3 + dS3*ph0 - dS1*ph2 + dS2*ph1;     
      
      vectorInterim[I][0] = 0;
      vectorInterim[I][1] = fac*dp1;
      vectorInterim[I][2] = fac*dp2;
      vectorInterim[I][3] = fac*dp3;  
      count += 4;
    }
    fftw_execute(vectorInterimPlanForward);  
  
    phiTotalSphericalProjectedForceAnalysisPREC->addDataVectorFromfourierTrafoSPECIAL((Complex*) vectorInterim); 
    phiTotalSphericalProjectedForceAnalysisGLOBAL->addDataVectorFromfourierTrafoSPECIAL((Complex*) vectorInterim);
  }
  
  if (dSdPhi != ((double*)phiForceInternal)) {
    if (LogLevel>4) printf("Copying required...\n");  
    SSE_ZCopy(2*L0*L1*L2*L3, (Complex*) dSdPhi, 1, (Complex*) phiForceInternal, 1);
  }  

  vector4D* phiForceFFTresult = NULL;
  if (takeOverFFTresult) { 
    fftw_execute(phiForceInternalPlanForward);  
    double fac = 1.0 /sqrt(L0*L1*L2*L3);
    for (int I=0; I<L0*L1*L2*L3; I++) {
      phiForceInternal[I][0] *= fac;
      phiForceInternal[I][1] *= fac;
      phiForceInternal[I][2] *= fac;
      phiForceInternal[I][3] *= fac;  
    }
    phiForceFFTresult = phiForceInternal;
    
    if (dSdPhi != ((double*)phiForceInternal)) {
      if (LogLevel>4) printf("Copying required...\n");
      SSE_ZCopy(2*L0*L1*L2*L3, (Complex*) phiForceInternal, 1, (Complex*) dSdPhi, 1);
    }
  } else {
    fftw_execute(internalToInterimPlanForward);  
    double fac = 1.0 /sqrt(L0*L1*L2*L3);
    for (int I=0; I<L0*L1*L2*L3; I++) {
      vectorInterim[I][0] *= fac;
      vectorInterim[I][1] *= fac;
      vectorInterim[I][2] *= fac;
      vectorInterim[I][3] *= fac;  
    }
    phiForceFFTresult = vectorInterim;
  }
      
  if (!tuneMode) {
    if ((forceType == 0) || ((fermiOps->getYN()==0) && (forceType==1))) {
      phiTotalForceAnalysisPREC->addDataVectorFromfourierTrafoSPECIAL((Complex*) phiForceFFTresult);
      phiTotalForceAnalysisGLOBAL->addDataVectorFromfourierTrafoSPECIAL((Complex*) phiForceFFTresult);  
    } 
    if (forceType == 1) {
      phiHiggsForceAnalysisPREC->addDataVectorFromfourierTrafoSPECIAL((Complex*) phiForceFFTresult);
      phiHiggsForceAnalysisGLOBAL->addDataVectorFromfourierTrafoSPECIAL((Complex*) phiForceFFTresult);  
    } 
    if (forceType == 2) {
      phiFermionForceAnalysisPREC[fermionForceSubCat]->addDataVectorFromfourierTrafoSPECIAL((Complex*) phiForceFFTresult);
      phiFermionForceAnalysisGLOBAL[fermionForceSubCat]->addDataVectorFromfourierTrafoSPECIAL((Complex*) phiForceFFTresult);  
    } 
    if ((forceType<0) || (forceType>2)) {
      printf("CalcPhiForceFourierComponents: Unknown ForceType!!!\n");
      exit(0);
    }  
  }
  addPerformanceProfilingItem("Propagator::calcPhiForceFourierComponents", performanceProfilerStartCycle, 0);  
}


void Propagator::calcPhiChangeFourierComponents() {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  if (LogLevel>4) printf("Computing Phi - Change - Fourier - Components...\n");

  int I;
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();

  double fac = 1.0 /sqrt(L0*L1*L2*L3);
  for (I=0; I<4*L0*L1*L2*L3; I++) {
    ((double*)vectorInterim)[I] = fac*(((double*)phiField)[I] - ((double*)phiFieldOLD)[I]);  
  }
  if (LogLevel>4) {
    double ccc[4];  
    ccc[0] = 0;
    ccc[1] = 0;
    ccc[2] = 0;
    ccc[3] = 0;
    for (I=0; I<L0*L1*L2*L3; I++) {
      ccc[0] += vectorInterim[I][0];
      ccc[1] += vectorInterim[I][1];
      ccc[2] += vectorInterim[I][2];
      ccc[3] += vectorInterim[I][3];    
    }
    printf("PhiChange:%f %f %f %f\n",ccc[0],ccc[1],ccc[2],ccc[3]);
  }

  fftw_execute(vectorInterimPlanBackward);  
  
  if (!tuneMode) {
    if (momentumMassesDetermined) {
      phiChangeAnalysisFACC->addDataVectorFromfourierTrafoSPECIAL((Complex*) vectorInterim);  
    } else {
      phiChangeAnalysisNoFACC->addDataVectorFromfourierTrafoSPECIAL((Complex*) vectorInterim);    
    }
  }
  addPerformanceProfilingItem("Propagator::calcPhiChangeFourierComponents", performanceProfilerStartCycle, 0);  
}


void Propagator::savePhiTotalSphericalProjectedForceFourierComponentsPREC(char* fileName) {
  phiTotalSphericalProjectedForceAnalysisPREC->saveData(fileName);  
}


void Propagator::savePhiTotalSphericalProjectedForceFourierComponentsGLOBAL(char* fileName) {
  phiTotalSphericalProjectedForceAnalysisGLOBAL->saveData(fileName);  
}


void Propagator::savePhiTotalForceFourierComponentsPREC(char* fileName) {
  phiTotalForceAnalysisPREC->saveData(fileName);  
}


void Propagator::savePhiTotalForceFourierComponentsGLOBAL(char* fileName) {
  phiTotalForceAnalysisGLOBAL->saveData(fileName);  
}


void Propagator::savePhiHiggsForceFourierComponentsPREC(char* fileName) {
  phiHiggsForceAnalysisPREC->saveData(fileName);  
}


void Propagator::savePhiHiggsForceFourierComponentsGLOBAL(char* fileName) {
  phiHiggsForceAnalysisGLOBAL->saveData(fileName);  
}


void Propagator::savePhiFermionForceFourierComponentsPREC(char* fileName) {
  int I;
  for (I=0; I<getFermionForceSubCategoryCount(); I++) {
    char* fileNameNeu = new char[600];
    snprintf(fileNameNeu,600,"%sSubCat_%d.dat",fileName,I);
    phiFermionForceAnalysisPREC[I]->saveData(fileNameNeu);  
    delete[] fileNameNeu;
  }
}


void Propagator::savePhiFermionForceFourierComponentsGLOBAL(char* fileName) {
  int I;
  for (I=0; I<getFermionForceSubCategoryCount(); I++) {
    char* fileNameNeu = new char[600];
    snprintf(fileNameNeu,600,"%sSubCat_%d.dat",fileName,I);
    phiFermionForceAnalysisGLOBAL[I]->saveData(fileNameNeu);  
    delete[] fileNameNeu;
  }
}


void Propagator::savePhiChangeFourierComponentsNoFACC(char* fileName) {
  phiChangeAnalysisNoFACC->saveData(fileName);  
}


void Propagator::savePhiChangeFourierComponentsFACC(char* fileName) {
  phiChangeAnalysisFACC->saveData(fileName);  
}


void Propagator::writeMomentumMasses(char* fileName) {
  FILE* file = fopen(fileName,"w");
  int I;

  if (LogLevel>4) printf("Writing Momenta Masses to disk.\n");

  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  
  for (I=0; I<L0*L1*L2*L3; I++) {
    double latMomsqr = phiTotalForceAnalysisPREC->getLatMomSqrFromIndex(I);
    fprintf(file,"%1.15f %1.15f\n",latMomsqr,MomentaMasses[I]);
  }
  
  fclose(file);
}


void Propagator::loadPhiTotalSphericalProjectedForceFourierComponentsPREC(char* fileName) {
  phiTotalSphericalProjectedForceAnalysisPREC->loadData(fileName);  
}


void Propagator::loadPhiTotalSphericalProjectedForceFourierComponentsGLOBAL(char* fileName) {
  phiTotalSphericalProjectedForceAnalysisGLOBAL->loadData(fileName);  
}


void Propagator::loadPhiTotalForceFourierComponentsPREC(char* fileName) {
  phiTotalForceAnalysisPREC->loadData(fileName);  
}


void Propagator::loadPhiTotalForceFourierComponentsGLOBAL(char* fileName) {
  phiTotalForceAnalysisGLOBAL->loadData(fileName);  
}


void Propagator::loadPhiHiggsForceFourierComponentsPREC(char* fileName) {
  phiHiggsForceAnalysisPREC->loadData(fileName);  
}


void Propagator::loadPhiHiggsForceFourierComponentsGLOBAL(char* fileName) {
  phiHiggsForceAnalysisGLOBAL->loadData(fileName);  
}


void Propagator::loadPhiFermionForceFourierComponentsPREC(char* fileName) {
  int I;
  for (I=0; I<getFermionForceSubCategoryCount(); I++) {
    char* fileNameNeu = new char[600];
    snprintf(fileNameNeu,600,"%sSubCat_%d.dat",fileName,I);
    phiFermionForceAnalysisPREC[I]->loadData(fileNameNeu);  
    delete[] fileNameNeu;
  }
}


void Propagator::loadPhiFermionForceFourierComponentsGLOBAL(char* fileName) {
  int I;
  for (I=0; I<getFermionForceSubCategoryCount(); I++) {
    char* fileNameNeu = new char[600];
    snprintf(fileNameNeu,600,"%sSubCat_%d.dat",fileName,I);
    phiFermionForceAnalysisGLOBAL[I]->loadData(fileNameNeu);  
    delete[] fileNameNeu;
  }
}
 

void Propagator::loadPhiChangeFourierComponentsNoFACC(char* fileName) {
  phiChangeAnalysisNoFACC->loadData(fileName);  
}


void Propagator::loadPhiChangeFourierComponentsFACC(char* fileName) {
  phiChangeAnalysisFACC->loadData(fileName);  
}


void Propagator::measure(double& measurePhiNorm, double& measureStaggeredPhiNorm, double& avgNorm, double& sigmaNorm) {
  int I1,I2,I3,I4;
  int count = 0;
  vector4D measurePhi, measureStaggeredPhi;
  vector4D measurePhiThirdMoment;
  double dummy;
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  
  measurePhi[0] = 0;
  measurePhi[1] = 0;
  measurePhi[2] = 0;
  measurePhi[3] = 0;
  measureStaggeredPhi[0] = 0;
  measureStaggeredPhi[1] = 0;
  measureStaggeredPhi[2] = 0;
  measureStaggeredPhi[3] = 0;
  measurePhiThirdMoment[0] = 0;
  measurePhiThirdMoment[1] = 0;
  measurePhiThirdMoment[2] = 0;
  measurePhiThirdMoment[3] = 0;  
  
  avgNorm = 0;
  sigmaNorm = 0;
  
  count = 0;
  for (I1=0; I1<L0; I1++) {
    for (I2=0; I2<L1; I2++) {
      for (I3=0; I3<L2; I3++) {
        for (I4=0; I4<L3; I4++) {
          measurePhi[0] += phiField[count][0];
          measurePhi[1] += phiField[count][1];
          measurePhi[2] += phiField[count][2];
          measurePhi[3] += phiField[count][3];
	  
	  avgNorm += dummy = sqrt( phiField[count][0]*phiField[count][0] 
	                          +phiField[count][1]*phiField[count][1]
	                          +phiField[count][2]*phiField[count][2]
	                          +phiField[count][3]*phiField[count][3]);
          sigmaNorm += dummy*dummy;
	  
	  measurePhiThirdMoment[0] += dummy*dummy*phiField[count][0];
	  measurePhiThirdMoment[1] += dummy*dummy*phiField[count][1];
	  measurePhiThirdMoment[2] += dummy*dummy*phiField[count][2];
	  measurePhiThirdMoment[3] += dummy*dummy*phiField[count][3];
	  
	  double stagFac = 1;
	  if (((I1+I2+I3+I4) % 2) == 1) stagFac = -1;
          measureStaggeredPhi[0] += stagFac*phiField[count][0];
          measureStaggeredPhi[1] += stagFac*phiField[count][1];
          measureStaggeredPhi[2] += stagFac*phiField[count][2];
          measureStaggeredPhi[3] += stagFac*phiField[count][3];
	  count++;
	}
      }
    }
  }
  double norm = L0*L1*L2*L3;
  measurePhi[0] = measurePhi[0]/norm;
  measurePhi[1] = measurePhi[1]/norm;
  measurePhi[2] = measurePhi[2]/norm;
  measurePhi[3] = measurePhi[3]/norm;
  measureStaggeredPhi[0] = measureStaggeredPhi[0]/norm;
  measureStaggeredPhi[1] = measureStaggeredPhi[1]/norm;
  measureStaggeredPhi[2] = measureStaggeredPhi[2]/norm;
  measureStaggeredPhi[3] = measureStaggeredPhi[3]/norm;
  measurePhiThirdMoment[0] = measurePhiThirdMoment[0]/norm;
  measurePhiThirdMoment[1] = measurePhiThirdMoment[1]/norm;
  measurePhiThirdMoment[2] = measurePhiThirdMoment[2]/norm;
  measurePhiThirdMoment[3] = measurePhiThirdMoment[3]/norm;

  avgNorm /= norm;
  sigmaNorm = sqrt(sigmaNorm/norm - avgNorm*avgNorm);
  
  measurePhiNorm = measurePhi[0]*measurePhi[0]
                 + measurePhi[1]*measurePhi[1]
                 + measurePhi[2]*measurePhi[2]
                 + measurePhi[3]*measurePhi[3];
  measureStaggeredPhiNorm = measureStaggeredPhi[0]*measureStaggeredPhi[0]
                          + measureStaggeredPhi[1]*measureStaggeredPhi[1]
                          + measureStaggeredPhi[2]*measureStaggeredPhi[2]
                          + measureStaggeredPhi[3]*measureStaggeredPhi[3];
		 
  measurePhiNorm = sqrt(measurePhiNorm);
  measureStaggeredPhiNorm = sqrt(measureStaggeredPhiNorm);
  
  if (LogLevel>1) {
    printf("m = %1.15f\n", measurePhiNorm);
    printf("s = %1.15f\n", measureStaggeredPhiNorm);
  }
  if (LogLevel>1) printf("\n");
}


void Propagator::calcPhiDerivativesOfSPhi() {
  calcPhiDerivativesOfSPhi((double*) phiForceInternal);
}


void Propagator::calcPhiDerivativesOfSPhi(double* dSdPhi) {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  int count = 0;
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();  
  double lamVal = lambda;
  if (isNaN(lambda)) {
    lamVal = 0;
  }
  
  for (int pos=0; pos<L0*L1*L2*L3; pos++) {
    double phiSqr = phiField[pos][0]*phiField[pos][0] 
                  + phiField[pos][1]*phiField[pos][1] 
                  + phiField[pos][2]*phiField[pos][2] 
                  + phiField[pos][3]*phiField[pos][3];

    int mu;  
    for (mu=0; mu<4; mu++) {
      double Phi = phiField[pos][mu];
      double neighbourSum = phiField[neighbours[pos].neighbourIndex[0]][mu] + phiField[neighbours[pos].neighbourIndex[1]][mu]
                          + phiField[neighbours[pos].neighbourIndex[2]][mu] + phiField[neighbours[pos].neighbourIndex[3]][mu]
                          + phiField[neighbours[pos].neighbourIndex[4]][mu] + phiField[neighbours[pos].neighbourIndex[5]][mu]
                          + phiField[neighbours[pos].neighbourIndex[6]][mu] + phiField[neighbours[pos].neighbourIndex[7]][mu];
  		    
      dSdPhi[count] = 2 * Phi - 2*kappa*neighbourSum + 4*lamVal*(phiSqr-Nf)*Phi;
      
      count++;
    }    
  }
  
  if (explicitCurrent != 0) {
    count = 0;
    for (int pos=0; pos<L0*L1*L2*L3; pos++) {
      dSdPhi[count] -= explicitCurrent;
      count += 4;    
    }
  }
  
  if ((ModelExtension_lambda6 != 0) ||(ModelExtension_lambda8 != 0) ||(ModelExtension_lambda10 != 0)) {
    count = 0;
    for (int pos=0; pos<L0*L1*L2*L3; pos++) {
      double phiSqr = phiField[pos][0]*phiField[pos][0] 
                    + phiField[pos][1]*phiField[pos][1] 
                    + phiField[pos][2]*phiField[pos][2] 
                    + phiField[pos][3]*phiField[pos][3];

      int mu;  
      for (mu=0; mu<4; mu++) {
        double Phi = phiField[pos][mu];
  		 
        double dummy = phiSqr*phiSqr;
        dSdPhi[count] += 6*ModelExtension_lambda6*dummy*Phi;
        dummy *= phiSqr;
        dSdPhi[count] += 8*ModelExtension_lambda8*dummy*Phi;
        dummy *= phiSqr;
        dSdPhi[count] += 10*ModelExtension_lambda10*dummy*Phi;
      
        count++;
      }    
    }
  }

  if ((ModelExtension_c6 != 0) ||(ModelExtension_c8 != 0) ||(ModelExtension_c10 != 0)) {
    fftw_execute(phiFieldToInterimPlanForward);  
    double norm = 2.0 / (L0*L1*L2*L3);
    for (int pos=0; pos<L0*L1*L2*L3; pos++) {
      double fac = MomentumSpaceActionCoefficients[pos] * norm;
      vectorInterim[pos][0] *= fac;
      vectorInterim[pos][1] *= fac;
      vectorInterim[pos][2] *= fac;
      vectorInterim[pos][3] *= fac;
    }
    
    fftw_execute(vectorInterimPlanBackward);  
    count = 0;
    for (int pos=0; pos<L0*L1*L2*L3; pos++) {
      dSdPhi[count+0] += vectorInterim[pos][0];
      dSdPhi[count+1] += vectorInterim[pos][1];
      dSdPhi[count+2] += vectorInterim[pos][2];
      dSdPhi[count+3] += vectorInterim[pos][3];
      count += 4;
    }
  }
  
  addPerformanceProfilingItem("Propagator::calcPhiDerivativesOfSPhi", performanceProfilerStartCycle, 0);  
}


void Propagator::calcFullPhiDerivatives(double propTol, int forceSelect) {
  calcFullPhiDerivatives((double*) phiForceInternal, propTol, forceSelect);
}


void Propagator::calcFullPhiDerivatives(double* dSdPhi, double propTol, int forceSelect) {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  int I,I2;
  
  if ((forceSelect ==3) && (dSdPhi == ((double*)phiForceInternal))) {
    printf("CalcFullPhiDerivatives: ERROR ==> dSdPhi == PhiForceInternal and forceSelect=3!!!\n");
    exit(0);
  }
  
  if ((forceSelect == 0) || (forceSelect == 1)) {
    calcPhiDerivativesOfSPhi(dSdPhi);
  } else {
    for (I=0; I<4*L0*L1*L2*L3; I++) {
      dSdPhi[I] = 0;
    }
  }
  
  if (forceSelect == 3) {
    calcPhiDerivativesOfSPhi((double*) phiForceInternal);
  }
  
  if ((forceSelect == 0) || (forceSelect == 2) || (forceSelect == 3)) {
    if (fermiOps->getYN()>0) {
      threadedExecute(Propagator_FORCE, propTol);
      for (I=0; I<forceCount; I++) {
        double* dummy_dS = (double*) forces[I]->getdSdPhi();
        for (I2=0; I2<4*L0*L1*L2*L3; I2++) dSdPhi[I2] += 0.5 * dummy_dS[I2];
      }
      if (forceSelect == 3) {
        double* dummy_dS = (double*) phiForceInternal;
        for (I2=0; I2<4*L0*L1*L2*L3; I2++) dummy_dS[I2] += dSdPhi[I2];
        if (phiForceFourierType>0) {
          calcPhiForceFourierComponents(dummy_dS, 0, getFermionForceSubCategory(propTol), !SphericalMode);
        }  
      }
    }  
  }
    
  if (phiForceFourierType>0) {
    if (forceSelect == 3) forceSelect = 2;
    calcPhiForceFourierComponents(dSdPhi, forceSelect, getFermionForceSubCategory(propTol), !SphericalMode);
  }  
  addPerformanceProfilingItem("Propagator::calcFullPhiDerivatives", performanceProfilerStartCycle, 0);  
}


double Propagator::calcPhiAction() {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  // ... Phi - Contribution  for Lambda < Infinity
  double S = 0;
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  double lamVal = lambda;
  if (isNaN(lambda)) {
    lamVal = 0;
  }
  
  for (int pos=0; pos<L0*L1*L2*L3; pos++) {
    double phiSqr = phiField[pos][0]*phiField[pos][0] 
                  + phiField[pos][1]*phiField[pos][1] 
                  + phiField[pos][2]*phiField[pos][2] 
                  + phiField[pos][3]*phiField[pos][3];
  
    S += phiSqr;
    S += lamVal * (phiSqr - Nf) * (phiSqr - Nf);
    S -= explicitCurrent * phiField[pos][0];
    S += ModelExtension_lambda6  * phiSqr*phiSqr*phiSqr;
    S += ModelExtension_lambda8  * phiSqr*phiSqr*phiSqr*phiSqr;
    S += ModelExtension_lambda10 * phiSqr*phiSqr*phiSqr*phiSqr*phiSqr;
    

    int mu;  
    for (mu=0; mu<4; mu++) {
      double Phi = phiField[pos][mu];
      double neighbourSum = phiField[neighbours[pos].neighbourIndex[0]][mu] + phiField[neighbours[pos].neighbourIndex[1]][mu]
                          + phiField[neighbours[pos].neighbourIndex[2]][mu] + phiField[neighbours[pos].neighbourIndex[3]][mu]
                          + phiField[neighbours[pos].neighbourIndex[4]][mu] + phiField[neighbours[pos].neighbourIndex[5]][mu]
                          + phiField[neighbours[pos].neighbourIndex[6]][mu] + phiField[neighbours[pos].neighbourIndex[7]][mu];
  		    
      S -= kappa * Phi * neighbourSum;
    }
  }
  
  if ((ModelExtension_c6 != 0) ||(ModelExtension_c8 != 0) ||(ModelExtension_c10 != 0)) {
    fftw_execute(phiFieldToInterimPlanForward);  
    double momS = 0;
    
    for (int pos=0; pos<L0*L1*L2*L3; pos++) {
      double phiMomSqr = vectorInterim[pos][0]*vectorInterim[pos][0] 
                       + vectorInterim[pos][1]*vectorInterim[pos][1] 
                       + vectorInterim[pos][2]*vectorInterim[pos][2] 
                       + vectorInterim[pos][3]*vectorInterim[pos][3];
    
      momS += phiMomSqr * MomentumSpaceActionCoefficients[pos];
    }
    momS /= (L0*L1*L2*L3);
    S +=  momS;
  }
  
  addPerformanceProfilingItem("Propagator::calcPhiAction", performanceProfilerStartCycle, 0);
  return S;
}


double Propagator::calcPhiMomentumAction() {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  double S = 0;
  Complex dummy;  
  checkValidityOfSetting();
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  
  dummy.x = 0;
  double* p = &(momentaField[0][0]);
  int count = 0;
  if (phiForceFourierType>0) {
    if (!SphericalMode) {
      for (int I=0; I<L0*L1*L2*L3; I++) {
        double f = 1.0 / MomentaMasses[I];
        dummy.x += f*p[count+0]*p[count+0];
        dummy.x += f*p[count+1]*p[count+1];
        dummy.x += f*p[count+2]*p[count+2];
        dummy.x += f*p[count+3]*p[count+3];    
        count += 4;
      }  
    } else {
      fftw_execute(momentaToInterimPlanForward);
      p = &(vectorInterim[0][0]);
      double fac = (1.0 + SphericalZeta) / (L0*L1*L2*L3);      
      for (int I=0; I<L0*L1*L2*L3; I++) {
        double f = fac / MomentaMasses[I];
        dummy.x += f*p[count+0]*p[count+0];
        dummy.x += f*p[count+1]*p[count+1];
        dummy.x += f*p[count+2]*p[count+2];
        dummy.x += f*p[count+3]*p[count+3];   
        count += 4;
      }  
    }
  } else {  
    if (!SphericalMode) {
      for (int I=0; I<L0*L1*L2*L3; I++) {
        dummy.x += p[count+0]*p[count+0];
        dummy.x += p[count+1]*p[count+1];
        dummy.x += p[count+2]*p[count+2];
        dummy.x += p[count+3]*p[count+3];    
        count += 4;
      }
    } else {
      double fac = (1.0 + SphericalZeta);
      for (int I=0; I<L0*L1*L2*L3; I++) {
        dummy.x += fac*p[count+0]*p[count+0];
        dummy.x += fac*p[count+1]*p[count+1];
        dummy.x += fac*p[count+2]*p[count+2];
        dummy.x += fac*p[count+3]*p[count+3];    
        count += 4;
      }    
    }
  }
  S = 0.5*gamma*dummy.x;
  addPerformanceProfilingItem("Propagator::calcPhiMomentumAction", performanceProfilerStartCycle, 0);  
  return S;
}


void Propagator::savePhiFields() {
  if (LogLevel>4) printf("Saving Phi fields.\n");
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
//  cblas_zcopy(2*N*N*N*N, phiField, 1, phiFieldOLD, 1);
  SSE_ZCopy(2*L0*L1*L2*L3, (Complex*)(&(phiField[0][0])), 1, (Complex*)(&(phiFieldOLD[0][0])), 1);
//  cblas_zcopy(2*N*N*N*N, momentaField, 1, momentaFieldOLD, 1);
  SSE_ZCopy(2*L0*L1*L2*L3, (Complex*)(&(momentaField[0][0])), 1, (Complex*)(&(momentaFieldOLD[0][0])), 1);
  SoldOLD =Sold;
}


void Propagator::getCopyOfSavedPhi(double* copy) {
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  SSE_ZCopy(2*L0*L1*L2*L3, (Complex*)(&(phiFieldOLD[0][0])), 1, (Complex*)(&(copy[0])), 1);
}


void Propagator::restorePhiFields() {
  if (LogLevel>4) printf("Restoring Phi fields.\n");
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
//  cblas_zcopy(2*N*N*N*N, phiFieldOLD, 1, phiField, 1);
  SSE_ZCopy(2*L0*L1*L2*L3, (Complex*)(&(phiFieldOLD[0][0])), 1, (Complex*)(&(phiField[0][0])), 1);
//  cblas_zcopy(2*N*N*N*N, momentaFieldOLD, 1, momentaField, 1);
  SSE_ZCopy(2*L0*L1*L2*L3, (Complex*)(&(momentaFieldOLD[0][0])), 1, (Complex*)(&(momentaField[0][0])), 1);
  Sold = SoldOLD;
}


void Propagator::MiniPhiMomentumStep(double eps) {
  MiniPhiMomentumStep((double*) phiForceInternal, eps) ;
}


void Propagator::MiniPhiMomentumStep(double* dSdPhi, double eps) {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  checkValidityOfSetting();

  if (!SphericalMode) {
    for (int I=0; I<4*L0*L1*L2*L3; I++) ((double*)momentaField)[I] -=  eps * dSdPhi[I];
  } else {
    int count = 0;
    double ph0,ph1,ph2,ph3;
    double dS0,dS1,dS2,dS3;
    double dp1,dp2,dp3;
    double epsFac1 = eps * 1.0/(1.0 + SphericalZeta);
    double epsFac2 = eps * SphericalZeta/(1.0 + SphericalZeta);    
    for (int I=0; I<L0*L1*L2*L3; I++) {
      ph0 = phiField[I][0];
      ph1 = phiField[I][1];
      ph2 = phiField[I][2];
      ph3 = phiField[I][3];
      dS0 = dSdPhi[count+0];
      dS1 = dSdPhi[count+1];
      dS2 = dSdPhi[count+2];
      dS3 = dSdPhi[count+3];
      
      dp1 = -dS0*ph1 + dS1*ph0 - dS2*ph3 + dS3*ph2;     
      dp2 = -dS0*ph2 + dS2*ph0 + dS1*ph3 - dS3*ph1;     
      dp3 = -dS0*ph3 + dS3*ph0 - dS1*ph2 + dS2*ph1;     

      momentaField[I][0] -= epsFac2*dS0;
      momentaField[I][1] -= epsFac2*dS1 + epsFac1 * dp1;
      momentaField[I][2] -= epsFac2*dS2 + epsFac1 * dp2;
      momentaField[I][3] -= epsFac2*dS3 + epsFac1 * dp3;
            
      count += 4;
    }
  }
  addPerformanceProfilingItem("Propagator::MiniPhiMomentumStep", performanceProfilerStartCycle, 0);  
}


void Propagator::MiniPhiStep(double eps) {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  checkValidityOfSetting();
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  double epsGamma = eps * gamma;
  
  if (!SphericalMode) {
    if (phiForceFourierType>0) {
      int count = 0;
      double f = 1.0 / sqrt(L0*L1*L2*L3);
      for (int I=0; I<L0*L1*L2*L3; I++) {
        double fac = f / MomentaMasses[I];
        ((double*)vectorInterim)[count+0] = fac*(((double*)momentaField)[count+0]);  
        ((double*)vectorInterim)[count+1] = fac*(((double*)momentaField)[count+1]);  
        ((double*)vectorInterim)[count+2] = fac*(((double*)momentaField)[count+2]);  
        ((double*)vectorInterim)[count+3] = fac*(((double*)momentaField)[count+3]);  
        count += 4;
      }
      fftw_execute(vectorInterimPlanBackward);  
      for (int I=0; I<4*L0*L1*L2*L3; I++) ((double*)phiField)[I] += epsGamma * ((double*)vectorInterim)[I];  
    } else {
      for (int I=0; I<4*L0*L1*L2*L3; I++) ((double*)phiField)[I] += epsGamma * ((double*)momentaField)[I];  
    }
  } else {
    vector4D* momenta = momentaField;
    if (phiForceFourierType>0) {
      int count = 0;
      double f = 1.0 / (L0*L1*L2*L3);      
      fftw_execute(momentaToInterimPlanForward);      
      for (int I=0; I<L0*L1*L2*L3; I++) {
        double fac = f / MomentaMasses[I];
        ((double*)vectorInterim)[count+0] = fac*(((double*)vectorInterim)[count+0]);  
        ((double*)vectorInterim)[count+1] = fac*(((double*)vectorInterim)[count+1]);  
        ((double*)vectorInterim)[count+2] = fac*(((double*)vectorInterim)[count+2]);  
        ((double*)vectorInterim)[count+3] = fac*(((double*)vectorInterim)[count+3]);  
        count += 4;
      }
      fftw_execute(vectorInterimPlanBackward);  
      momenta = vectorInterim;
    }

    double ph0,ph1,ph2,ph3; 
    Quat q;
    Quat q2;
    double epsGammaZetaHalf = 0.5*eps * gamma * SphericalZeta;
    for (int I=0; I<L0*L1*L2*L3; I++) {    
      phiField[I][0] += epsGammaZetaHalf * momenta[I][0];  
      phiField[I][1] += epsGammaZetaHalf * momenta[I][1];  
      phiField[I][2] += epsGammaZetaHalf * momenta[I][2];  
      phiField[I][3] += epsGammaZetaHalf * momenta[I][3];      
    
      ph0 = phiField[I][0];
      ph1 = phiField[I][1];
      ph2 = phiField[I][2];
      ph3 = phiField[I][3];
      q.x0 = 0;
      q.x1 = epsGamma*momenta[I][1];
      q.x2 = epsGamma*momenta[I][2];
      q.x3 = epsGamma*momenta[I][3];

      q2 = exp(q);
      
      phiField[I][0] =  q2.x0*ph0 - q2.x1*ph1 - q2.x2*ph2 - q2.x3*ph3;
      phiField[I][1] =  q2.x0*ph1 + q2.x1*ph0 + q2.x2*ph3 - q2.x3*ph2;
      phiField[I][2] =  q2.x0*ph2 + q2.x2*ph0 - q2.x1*ph3 + q2.x3*ph1;
      phiField[I][3] =  q2.x0*ph3 + q2.x3*ph0 + q2.x1*ph2 - q2.x2*ph1;      
      
      phiField[I][0] += epsGammaZetaHalf * momenta[I][0];  
      phiField[I][1] += epsGammaZetaHalf * momenta[I][1];  
      phiField[I][2] += epsGammaZetaHalf * momenta[I][2];  
      phiField[I][3] += epsGammaZetaHalf * momenta[I][3];      
    }
        
    if (isNaN(lambda)) {
      double norm, fac;
      for (int I=0; I<L0*L1*L2*L3; I++) {
        norm = phiField[I][0]*phiField[I][0]
             + phiField[I][1]*phiField[I][1]
             + phiField[I][2]*phiField[I][2]
             + phiField[I][3]*phiField[I][3];
        fac = sqrt(Nf / norm);
      
        phiField[I][0] *= fac;
        phiField[I][1] *= fac;
        phiField[I][2] *= fac;
        phiField[I][3] *= fac;
      }
    }
  }
  
  phiFieldWasChanged();
  addPerformanceProfilingItem("Propagator::MiniPhiStep", performanceProfilerStartCycle, 0);  
}


void Propagator::LeapOmelyanPhiPropagation(int iterations, double epsilon, double lambda, double rho, double theta, double mu) {
  int I2;  
  double epsilonHalf = 0.5 * epsilon;
  double epsilonLambda = lambda * epsilon;
  double epsilonRho = rho * epsilon;
  double epsilonMu = mu * epsilon;
  double epsilonTheta = theta * epsilon;
  double epsilonHalfOneMinusTwoLambdaPlusTheta = 0.5*(1.0 - 2.0*(lambda + theta))*epsilon;
  double epsilonOneMinusTwoMuPlusRho = (1.0 - 2.0*(mu + rho))*epsilon;
  
  double twoEpsilonLambda = 2.0 * lambda * epsilon;
  double epsilonOneMinusTwoLambda = (1.0 - 2.0*lambda)*epsilon;

  calcFullPhiDerivatives(0, 1); //calculates also omegaForces
  MiniPhiMomentumStep(epsilonLambda);

  for (I2=0; I2<iterations; I2++) {
    if ((lambda==0.5) && (theta==0) && (rho==0) && (mu==0)) {
      //Normaler Leap-Frog Schritt

      MiniPhiStep(epsilon);

    } else if ((theta==0) && (rho==0) && (mu==0)) {
      //Omelyan - Schritt mit Order 2

      MiniPhiStep(epsilonHalf);

      calcFullPhiDerivatives(0, 1);
      MiniPhiMomentumStep(epsilonOneMinusTwoLambda);

      MiniPhiStep(epsilonHalf);

    } else {
      //Omelyan - Schritt mit Order 4

      MiniPhiStep(epsilonRho);

      calcFullPhiDerivatives(0, 1);
      MiniPhiMomentumStep(epsilonTheta);

      MiniPhiStep(epsilonMu);

      calcFullPhiDerivatives(0, 1);
      MiniPhiMomentumStep(epsilonHalfOneMinusTwoLambdaPlusTheta);

      MiniPhiStep(epsilonOneMinusTwoMuPlusRho);

      calcFullPhiDerivatives(0, 1);
      MiniPhiMomentumStep(epsilonHalfOneMinusTwoLambdaPlusTheta);

      MiniPhiStep(epsilonMu);
      
      calcFullPhiDerivatives(0, 1);
      MiniPhiMomentumStep(epsilonTheta);

      MiniPhiStep(epsilonRho);
    }

    calcFullPhiDerivatives(0, 1);
    MiniPhiMomentumStep(twoEpsilonLambda);
  }
  //Bring momenta on integer step numbers again
  MiniPhiMomentumStep(-epsilonLambda);
}


void Propagator::LeapFrogPhiPropagation(int iterations, double epsilon) {
  double lambda = 0.5;
  LeapOmelyanPhiPropagation(iterations, epsilon, lambda, 0, 0, 0);
}


void Propagator::OmelyanO2PhiPropagation(int iterations, double epsilon) {
  double lambda = 0.1931833275037836;
  LeapOmelyanPhiPropagation(iterations, epsilon, lambda, 0, 0, 0);
}


void Propagator::OmelyanO4PhiPropagation(int iterations, double epsilon) {
  double lambda = 0.08398315262876693;
  double rho    = 0.2539785108410595;
  double theta  = 0.6822365335719091;
  double mu     = -0.03230286765269967;
  LeapOmelyanPhiPropagation(iterations, epsilon, lambda, rho, theta, mu);
}


bool Propagator::LeapFrogMarkovStep(int iterations, double epsilon, double propTOL, double finalTOL) {
  double lambda = 0.5;
  return LeapOmelyanMarkovStep(iterations, epsilon, propTOL, finalTOL, lambda, 0, 0, 0);
}


bool Propagator::OmelyanO2MarkovStep(int iterations, double epsilon, double propTOL, double finalTOL) {
  double lambda = 0.1931833275037836;
  return LeapOmelyanMarkovStep(iterations, epsilon, propTOL, finalTOL, lambda, 0, 0, 0);
}


bool Propagator::OmelyanO4MarkovStep(int iterations, double epsilon, double propTOL, double finalTOL) {
  double lambda = 0.08398315262876693;
  double rho    = 0.2539785108410595;
  double theta  = 0.6822365335719091;
  double mu     = -0.03230286765269967;
  return LeapOmelyanMarkovStep(iterations, epsilon, propTOL, finalTOL, lambda, rho, theta, mu);
}



bool Propagator::multiTimeScaleMarkovStep(int LevelCount, int* iterations, double epsilon, int* integrators, int* subPolNr, double* propTOL, double finalTOL) {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  if (LevelCount<=0) {
    printf("ERROR: LevelCount <= 0 in MultiTimeScaleMarkovStep.\n");
    exit(0);
  }

  double* lambda = new double[LevelCount];
  double* rho = new double[LevelCount];
  double* theta = new double[LevelCount];
  double* mu = new double[LevelCount];

  int I;
  for (I=0; I<LevelCount; I++) {
    if (integrators[I] == 0) {
      lambda[I] = 0.5;
      rho[I]    = 0;
      theta[I]  = 0;
      mu[I]     = 0;
    } else if (integrators[I] == 1) {
      lambda[I] = 0.1931833275037836;
      rho[I]    = 0;
      theta[I]  = 0;
      mu[I]     = 0;
    } else if (integrators[I] == 2) {
      lambda[I] = 0.08398315262876693;
      rho[I]    = 0.2539785108410595;
      theta[I]  = 0.6822365335719091;
      mu[I]     = -0.03230286765269967;
    } else {
      printf("ERROR: Unknown Integrator in MultiTimeScaleMarkovStep.\n");
      exit(0);
    }
  }

  bool accepted = LeapOmelyanMultiTimeScaleMarkovStep(LevelCount-1, iterations, epsilon, subPolNr, propTOL, finalTOL, lambda, rho, theta, mu);

  delete[] lambda;
  delete[] rho;
  delete[] theta;
  delete[] mu;
  
  addPerformanceProfilingItem("Propagator::multiTimeScaleMarkovStep", performanceProfilerStartCycle, 0);
  return accepted;
}


void Propagator::multiTimeScalePropagation(int LevelCount, int* iterations, double epsilon, int* integrators, int* subPolNr, double* propTOL, double finalTOL) {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  if (LevelCount<=0) {
    printf("ERROR: LevelCount <= 0 in MultiTimeScalePropagation.\n");
    exit(0);
  }

  double* lambda = new double[LevelCount];
  double* rho = new double[LevelCount];
  double* theta = new double[LevelCount];
  double* mu = new double[LevelCount];

  int I;
  for (I=0; I<LevelCount; I++) {
    if (integrators[I] == 0) {
      lambda[I] = 0.5;
      rho[I]    = 0;
      theta[I]  = 0;
      mu[I]     = 0;
    } else if (integrators[I] == 1) {
      lambda[I] = 0.1931833275037836;
      rho[I]    = 0;
      theta[I]  = 0;
      mu[I]     = 0;
    } else if (integrators[I] == 2) {
      lambda[I] = 0.08398315262876693;
      rho[I]    = 0.2539785108410595;
      theta[I]  = 0.6822365335719091;
      mu[I]     = -0.03230286765269967;
    } else {
      printf("ERROR: Unknown Integrator in MultiTimeScaleMarkovStep.\n");
      exit(0);
    }
  }

  LeapOmelyanMultiTimeScalePropagation(LevelCount-1, iterations, epsilon, subPolNr, propTOL, finalTOL, lambda, rho, theta, mu);

  delete[] lambda;
  delete[] rho;
  delete[] theta;
  delete[] mu;
  addPerformanceProfilingItem("Propagator::multiTimeScalePropagation", performanceProfilerStartCycle, 0);
}


void Propagator::improvePreconditioningParametersFAST() {
  bool PrecUse;
  double PrecMold, PrecSold;
  fermiOps->getPreconditionerParameter(PrecUse, PrecMold, PrecSold);  

  if (!PrecUse) return;
  if (forceCount<=0) return;
  if (fermiOps->getYN()<=0) return;
  
  double avgPhi, avgStagPhi, avgNorm, sigmaNorm;
  measure(avgPhi, avgStagPhi, avgNorm, sigmaNorm);
  
  fermiOps->setPreconditioner(true, avgNorm, 0.0);
  
  double PrecMact, PrecSact;
  fermiOps->getPreconditionerParameter(PrecUse, PrecMact, PrecSact);
  if (LogLevel>1) printf("Fast Preconditioning Determination sets to: PrecM = %1.3f, PrecS = %1.3f\n", PrecMact, PrecSact); 
  PreconditionerWasChanged();
}
