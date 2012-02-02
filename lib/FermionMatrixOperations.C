#include "FermionMatrixOperations.h"


void FermionMatrixOperations::ini() {
  InputVectorFourierTransform = NULL;
  InterimVectorFourierTransform = NULL;

  xFFT = NULL;
  xFFTusage = false;
  xFFT_DistributedFFT_ThreadCount = 1;
  
  int I;
  for (I=0; I<FFTPlanDBMax; I++) {
    FFTPlanDataBase[I] = NULL;
    InputVectorDataBase[I] = NULL;
    OutputVectorDataBase[I] = NULL;
    BackwardForwardDataBase[I] = ExtremeFFT4D_Forward;
  }
  FFTDataBaseCounter = 0;
  
  Inverse_rest = NULL;
  Inverse_b = NULL;
  Inverse_p = NULL;
  Inverse_s = NULL;
  NeubergerDiracOpApplicationEWData = NULL;
  NeubergerDiracOpApplicationEWDataRescaled = NULL;
  NeubergerDiracPlusExplicitMassOpApplicationEWData = NULL;
  NeubergerDiracPlusExplicitMassOpApplicationEWDataRescaled = NULL;
  NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWData = NULL;
  NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWDataRescaled = NULL;
  NeubergerDiracOpApplicationSinPStdData = NULL;
  NeubergerDiracOpApplicationSinPDagData = NULL;  
  NeubergerWithChiStaticFermionInverseOpApplicationMainDiagonalData = NULL;
  NeubergerWithChiStaticFermionInverseOpApplicationSubDiagonalData = NULL;
  NeubergerWithChiStaticMdagMFermionInverseOpApplicationMainDiagonalData = NULL;
  NeubergerWithChiStaticMdagMFermionInverseOpApplicationSubDiagonalData = NULL;
  NeubergerWithChiStaticMdagMFermionInverseOpApplicationSSEAuxData = NULL;
  QPreconditionerDiagonalData = NULL;
  QPreconditionerDiagonalDataNormalized = NULL;
  inverseQPreconditionerDiagonalData = NULL;
  inverseQPreconditionerDiagonalDataNormalized = NULL;
  RPreconditionerDiagonalData = NULL;
  RPreconditionerDiagonalDataNormalized = NULL;
  inverseRPreconditionerDiagonalData = NULL;
  inverseRPreconditionerDiagonalDataNormalized = NULL;  
  RSQRPreconditionerDiagonalData = NULL;
  RSQRPreconditionerDiagonalDataNormalized = NULL;
  inverseRSQRPreconditionerDiagonalData = NULL;
  inverseRSQRPreconditionerDiagonalDataNormalized = NULL;  
  halfMomentumForwardFFTfactors = NULL;
  halfMomentumBackwardFFTfactors = NULL;
  
  Index_PlusPiPiPiPi = NULL;
  Index_PlusPiPiPiPiXtrSize = NULL;
  Index_PiModes = NULL;
  
  Preconditioner_M = NaN;
  Preconditioner_S = NaN;
  Preconditioner_Usage = false;
  QPreconditioner_mu = NaN;
  QPreconditioner_beta = NaN;
  QPreconditioner_Usage = false;
  RPreconditioner_m = NaN;
  RPreconditioner_f = NaN;
  RPreconditioner_Usage = false;  
  
  Parameter_r = 0.5;
  Parameter_rho = 1.0;
  Parameter_yN = 0.0;
  Parameter_MassSplit = 1.0;
  Parameter_ExplicitMass = 0.0;
  Parameter_AntiPeriodicBoundaryConditionInTime = false;
  NeubergerDiracOp = new NeubergerMatrix();
  WilsonDiracOp = new WilsonMatrix();
  oneDimSizeL0 = 0;
  oneDimSizeL1 = 0;
  oneDimSizeL2 = 0;
  oneDimSizeL3 = 0;
  vectorLength = 0;
  vectorLengthXtrSize = 0;
  
  threadedOps = NULL;
  DistributedInputVectorFourierTransform = NULL;  
  DistributedInterimVectorFourierTransform = NULL;    
  DistributedInverse_rest = NULL;
  DistributedInverse_b = NULL;
  DistributedInverse_p = NULL;
  DistributedInverse_s = NULL;

  Distributed_NeubergerDiracOpApplicationEWData = NULL;
  Distributed_NeubergerDiracOpApplicationEWDataRescaled = NULL;
  Distributed_NeubergerDiracPlusExplicitMassOpApplicationEWData = NULL;
  Distributed_NeubergerDiracPlusExplicitMassOpApplicationEWDataRescaled = NULL;
  Distributed_NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWData = NULL;
  Distributed_NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWDataRescaled = NULL;
  Distributed_NeubergerDiracOpApplicationSinPStdData = NULL;
  Distributed_NeubergerDiracOpApplicationSinPDagData = NULL;  
  Distributed_NeubergerWithChiStaticFermionInverseOpApplicationMainDiagonalData = NULL;
  Distributed_NeubergerWithChiStaticFermionInverseOpApplicationSubDiagonalData = NULL;
  Distributed_NeubergerWithChiStaticMdagMFermionInverseOpApplicationMainDiagonalData = NULL;
  Distributed_NeubergerWithChiStaticMdagMFermionInverseOpApplicationSubDiagonalData = NULL;
  Distributed_NeubergerWithChiStaticMdagMFermionInverseOpApplicationSSEAuxData = NULL;
  Distributed_QPreconditionerDiagonalData = NULL;
  Distributed_QPreconditionerDiagonalDataNormalized = NULL;
  Distributed_inverseQPreconditionerDiagonalData = NULL;
  Distributed_inverseQPreconditionerDiagonalDataNormalized = NULL;
  Distributed_RPreconditionerDiagonalData = NULL;
  Distributed_RPreconditionerDiagonalDataNormalized = NULL;
  Distributed_inverseRPreconditionerDiagonalData = NULL;
  Distributed_inverseRPreconditionerDiagonalDataNormalized = NULL;  
  Distributed_RSQRPreconditionerDiagonalData = NULL;
  Distributed_RSQRPreconditionerDiagonalDataNormalized = NULL;
  Distributed_inverseRSQRPreconditionerDiagonalData = NULL;
  Distributed_inverseRSQRPreconditionerDiagonalDataNormalized = NULL;
  Distributed_halfMomentumForwardFFTfactors = NULL;
  Distributed_halfMomentumBackwardFFTfactors = NULL;
    
  Distributed_Index_PlusPiPiPiPi = NULL;
  Distributed_Index_PlusPiPiPiPiXtrSize = NULL;
  Distributed_Index_PiModes = NULL;  
}


void FermionMatrixOperations::setSize(int sizeL0, int sizeL1, int sizeL2, int sizeL3) {
  if (sizeL0 <= 0) sizeL0 = 0;
  if (sizeL1 <= 0) sizeL1 = 0;
  if (sizeL2 <= 0) sizeL2 = 0;
  if (sizeL3 <= 0) sizeL3 = 0;
  if ((oneDimSizeL0 == sizeL0) && (oneDimSizeL1 == sizeL1) && (oneDimSizeL2 == sizeL2) && (oneDimSizeL3 == sizeL3)) return;    
  
  oneDimSizeL0 = sizeL0;
  oneDimSizeL1 = sizeL1;
  oneDimSizeL2 = sizeL2;
  oneDimSizeL3 = sizeL3;
    
  timeDirection = 0;
  oneDimSizeLargest = oneDimSizeL0;
  if (oneDimSizeL1 > oneDimSizeLargest) { 
    timeDirection = 1;
    oneDimSizeLargest = oneDimSizeL1;
  }
  if (oneDimSizeL2 > oneDimSizeLargest) {
    timeDirection = 2;
    oneDimSizeLargest = oneDimSizeL2;
  }
  if (oneDimSizeL3 > oneDimSizeLargest) {
    timeDirection = 3;
    oneDimSizeLargest = oneDimSizeL3;
  }
   
  vectorLength = sizeL0*sizeL1*sizeL2*sizeL3*8;
  vectorLengthXtrSize = (sizeL3+xtraSize3)*(sizeL2+xtraSize2)*(sizeL1+xtraSize1)*sizeL0*8;
  destroyFermionVector(InputVectorFourierTransform);
  InputVectorFourierTransform = createFermionVector();
  destroyFermionVector(InterimVectorFourierTransform);
  InterimVectorFourierTransform = createFermionVector();
  
  clearTrafoPlans();

  delete xFFT;
  xFFT = NULL;

  destroyFermionVector(Inverse_rest);
  Inverse_rest = createFermionVector();
  destroyFermionVector(Inverse_b);
  Inverse_b = createFermionVector();
  destroyFermionVector(Inverse_p);
  Inverse_p = createFermionVector();
  destroyFermionVector(Inverse_s);
  Inverse_s = createFermionVector();

  destroyFermionVector(NeubergerDiracOpApplicationEWData);
  NeubergerDiracOpApplicationEWData = createFermionVector(1);
  destroyFermionVector(NeubergerDiracOpApplicationEWDataRescaled);
  NeubergerDiracOpApplicationEWDataRescaled = createFermionVector(1);
  destroyFermionVector(NeubergerDiracPlusExplicitMassOpApplicationEWData);
  NeubergerDiracPlusExplicitMassOpApplicationEWData = createFermionVector(1);
  destroyFermionVector(NeubergerDiracPlusExplicitMassOpApplicationEWDataRescaled);
  NeubergerDiracPlusExplicitMassOpApplicationEWDataRescaled = createFermionVector(1);

  destroyFermionVector(NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWData);
  NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWData = createFermionVector(1);
  destroyFermionVector(NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWDataRescaled);
  NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWDataRescaled = createFermionVector(1);

  destroyFermionVector(QPreconditionerDiagonalData);
  QPreconditionerDiagonalData = createFermionVector(1);
  destroyFermionVector(QPreconditionerDiagonalDataNormalized);
  QPreconditionerDiagonalDataNormalized = createFermionVector(1);
  destroyFermionVector(inverseQPreconditionerDiagonalData);
  inverseQPreconditionerDiagonalData = createFermionVector(1);
  destroyFermionVector(inverseQPreconditionerDiagonalDataNormalized);
  inverseQPreconditionerDiagonalDataNormalized = createFermionVector(1);

  destroyFermionVector(RPreconditionerDiagonalData);
  RPreconditionerDiagonalData = createFermionVector(1);
  destroyFermionVector(RPreconditionerDiagonalDataNormalized);
  RPreconditionerDiagonalDataNormalized = createFermionVector(1);
  destroyFermionVector(inverseRPreconditionerDiagonalData);
  inverseRPreconditionerDiagonalData = createFermionVector(1);
  destroyFermionVector(inverseRPreconditionerDiagonalDataNormalized);
  inverseRPreconditionerDiagonalDataNormalized = createFermionVector(1);

  destroyFermionVector(RSQRPreconditionerDiagonalData);
  RSQRPreconditionerDiagonalData = createFermionVector(1);
  destroyFermionVector(RSQRPreconditionerDiagonalDataNormalized);
  RSQRPreconditionerDiagonalDataNormalized = createFermionVector(1);
  destroyFermionVector(inverseRSQRPreconditionerDiagonalData);
  inverseRSQRPreconditionerDiagonalData = createFermionVector(1);
  destroyFermionVector(inverseRSQRPreconditionerDiagonalDataNormalized);
  inverseRSQRPreconditionerDiagonalDataNormalized = createFermionVector(1);
  destroyFermionVector(halfMomentumForwardFFTfactors);
  halfMomentumForwardFFTfactors = createSuperAlignedComplex(oneDimSizeL0+oneDimSizeL1+oneDimSizeL2+oneDimSizeL3);
  destroyFermionVector(halfMomentumBackwardFFTfactors);
  halfMomentumBackwardFFTfactors = createSuperAlignedComplex(oneDimSizeL0+oneDimSizeL1+oneDimSizeL2+oneDimSizeL3);

  destroyFermionVector(NeubergerDiracOpApplicationSinPStdData);
  destroyFermionVector(NeubergerDiracOpApplicationSinPDagData);
  NeubergerDiracOpApplicationSinPStdData = createSuperAlignedComplex(2*4*128);
  NeubergerDiracOpApplicationSinPDagData = createSuperAlignedComplex(2*4*128);
  
  destroyFermionVector(NeubergerWithChiStaticFermionInverseOpApplicationMainDiagonalData);
  destroyFermionVector(NeubergerWithChiStaticFermionInverseOpApplicationSubDiagonalData);
  NeubergerWithChiStaticFermionInverseOpApplicationMainDiagonalData = createFermionVector(1);
  NeubergerWithChiStaticFermionInverseOpApplicationSubDiagonalData = createFermionVector(1);
  destroyFermionVector(NeubergerWithChiStaticMdagMFermionInverseOpApplicationMainDiagonalData);
  destroyFermionVector(NeubergerWithChiStaticMdagMFermionInverseOpApplicationSubDiagonalData);
  NeubergerWithChiStaticMdagMFermionInverseOpApplicationMainDiagonalData = createFermionVector(1);
  NeubergerWithChiStaticMdagMFermionInverseOpApplicationSubDiagonalData = createFermionVector(1);
  destroyFermionVector(NeubergerWithChiStaticMdagMFermionInverseOpApplicationSSEAuxData);
  NeubergerWithChiStaticMdagMFermionInverseOpApplicationSSEAuxData = createFermionVector(2);
  destroySuperAlignedLongInt(Index_PlusPiPiPiPi);
  Index_PlusPiPiPiPi = (long int*) createSuperAlignedComplex(1024);
  destroySuperAlignedLongInt(Index_PlusPiPiPiPiXtrSize);
  Index_PlusPiPiPiPiXtrSize = (long int*) createSuperAlignedComplex(1024);
  destroySuperAlignedLongInt(Index_PiModes);
  Index_PiModes = (long int*) createSuperAlignedComplex(15);
  
  setSizeDistributed(sizeL0, sizeL1, sizeL2, sizeL3);  
  
  generateOpApplicationData();  
}


void FermionMatrixOperations::setSizeDistributed(int sizeL0, int sizeL1, int sizeL2, int sizeL3) {
  if (threadedOps == NULL) return;

  destroyDistributedComplexVector(DistributedInputVectorFourierTransform);
  DistributedInputVectorFourierTransform = createDistributedFermionVector();  
  destroyDistributedComplexVector(DistributedInterimVectorFourierTransform);
  DistributedInterimVectorFourierTransform = createDistributedFermionVector();
  destroyDistributedComplexVector(DistributedInverse_rest);
  DistributedInverse_rest = createDistributedFermionVector();   
  destroyDistributedComplexVector(DistributedInverse_b);
  DistributedInverse_b = createDistributedFermionVector();    
  destroyDistributedComplexVector(DistributedInverse_p);
  DistributedInverse_p = createDistributedFermionVector();    
  destroyDistributedComplexVector(DistributedInverse_s);
  DistributedInverse_s = createDistributedFermionVector();
  
  destroyDistributedComplexVector(Distributed_NeubergerDiracOpApplicationEWData);
  Distributed_NeubergerDiracOpApplicationEWData = createDistributedUniformLatticeComplexVector(1);
  destroyDistributedComplexVector(Distributed_NeubergerDiracOpApplicationEWDataRescaled);
  Distributed_NeubergerDiracOpApplicationEWDataRescaled = createDistributedUniformLatticeComplexVector(1);
  destroyDistributedComplexVector(Distributed_NeubergerDiracPlusExplicitMassOpApplicationEWData);
  Distributed_NeubergerDiracPlusExplicitMassOpApplicationEWData = createDistributedUniformLatticeComplexVector(1);
  destroyDistributedComplexVector(Distributed_NeubergerDiracPlusExplicitMassOpApplicationEWDataRescaled);
  Distributed_NeubergerDiracPlusExplicitMassOpApplicationEWDataRescaled = createDistributedUniformLatticeComplexVector(1);

  destroyDistributedComplexVector(Distributed_NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWData);
  Distributed_NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWData = createDistributedUniformLatticeComplexVector(1);
  destroyDistributedComplexVector(Distributed_NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWDataRescaled);
  Distributed_NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWDataRescaled = createDistributedUniformLatticeComplexVector(1);  

  destroyDistributedComplexVector(Distributed_NeubergerDiracOpApplicationSinPStdData);
  Distributed_NeubergerDiracOpApplicationSinPStdData = threadedOps->allocateDistibutedMemory(2*4*128);
  destroyDistributedComplexVector(Distributed_NeubergerDiracOpApplicationSinPDagData);  
  Distributed_NeubergerDiracOpApplicationSinPDagData = threadedOps->allocateDistibutedMemory(2*4*128);
  
  destroyDistributedComplexVector(Distributed_NeubergerWithChiStaticFermionInverseOpApplicationMainDiagonalData);
  Distributed_NeubergerWithChiStaticFermionInverseOpApplicationMainDiagonalData = createDistributedUniformLatticeComplexVector(1);
  destroyDistributedComplexVector(Distributed_NeubergerWithChiStaticFermionInverseOpApplicationSubDiagonalData);
  Distributed_NeubergerWithChiStaticFermionInverseOpApplicationSubDiagonalData = createDistributedUniformLatticeComplexVector(1);
  destroyDistributedComplexVector(Distributed_NeubergerWithChiStaticMdagMFermionInverseOpApplicationMainDiagonalData);
  Distributed_NeubergerWithChiStaticMdagMFermionInverseOpApplicationMainDiagonalData = createDistributedUniformLatticeComplexVector(1);
  destroyDistributedComplexVector(Distributed_NeubergerWithChiStaticMdagMFermionInverseOpApplicationSubDiagonalData);
  Distributed_NeubergerWithChiStaticMdagMFermionInverseOpApplicationSubDiagonalData = createDistributedUniformLatticeComplexVector(1);  
  
  destroyDistributedComplexVector(Distributed_NeubergerWithChiStaticMdagMFermionInverseOpApplicationSSEAuxData);
  Distributed_NeubergerWithChiStaticMdagMFermionInverseOpApplicationSSEAuxData = createDistributedUniformLatticeComplexVector(2);
  
  destroyDistributedComplexVector(Distributed_QPreconditionerDiagonalData);
  Distributed_QPreconditionerDiagonalData = createDistributedUniformLatticeComplexVector(1);
  destroyDistributedComplexVector(Distributed_QPreconditionerDiagonalDataNormalized);
  Distributed_QPreconditionerDiagonalDataNormalized = createDistributedUniformLatticeComplexVector(1);
  destroyDistributedComplexVector(Distributed_inverseQPreconditionerDiagonalData);
  Distributed_inverseQPreconditionerDiagonalData = createDistributedUniformLatticeComplexVector(1);
  destroyDistributedComplexVector(Distributed_inverseQPreconditionerDiagonalDataNormalized);
  Distributed_inverseQPreconditionerDiagonalDataNormalized = createDistributedUniformLatticeComplexVector(1);
  
  destroyDistributedComplexVector(Distributed_RPreconditionerDiagonalData);
  Distributed_RPreconditionerDiagonalData = createDistributedUniformLatticeComplexVector(1);
  destroyDistributedComplexVector(Distributed_RPreconditionerDiagonalDataNormalized);
  Distributed_RPreconditionerDiagonalDataNormalized = createDistributedUniformLatticeComplexVector(1);
  destroyDistributedComplexVector(Distributed_inverseRPreconditionerDiagonalData);
  Distributed_inverseRPreconditionerDiagonalData = createDistributedUniformLatticeComplexVector(1);
  destroyDistributedComplexVector(Distributed_inverseRPreconditionerDiagonalDataNormalized);  
  Distributed_inverseRPreconditionerDiagonalDataNormalized = createDistributedUniformLatticeComplexVector(1);
  
  destroyDistributedComplexVector(Distributed_RSQRPreconditionerDiagonalData);
  Distributed_RSQRPreconditionerDiagonalData = createDistributedUniformLatticeComplexVector(1);
  destroyDistributedComplexVector(Distributed_RSQRPreconditionerDiagonalDataNormalized);
  Distributed_RSQRPreconditionerDiagonalDataNormalized = createDistributedUniformLatticeComplexVector(1);
  destroyDistributedComplexVector(Distributed_inverseRSQRPreconditionerDiagonalData);
  Distributed_inverseRSQRPreconditionerDiagonalData = createDistributedUniformLatticeComplexVector(1);
  destroyDistributedComplexVector(Distributed_inverseRSQRPreconditionerDiagonalDataNormalized);    
  Distributed_inverseRSQRPreconditionerDiagonalDataNormalized = createDistributedUniformLatticeComplexVector(1);

  destroyDistributedComplexVector(Distributed_halfMomentumForwardFFTfactors);
  Distributed_halfMomentumForwardFFTfactors = threadedOps->allocateDistibutedMemory(oneDimSizeL0+oneDimSizeL1+oneDimSizeL2+oneDimSizeL3);
  destroyDistributedComplexVector(Distributed_halfMomentumBackwardFFTfactors);
  Distributed_halfMomentumBackwardFFTfactors = threadedOps->allocateDistibutedMemory(oneDimSizeL0+oneDimSizeL1+oneDimSizeL2+oneDimSizeL3);

  destroyDistributedComplexVector(Distributed_Index_PlusPiPiPiPi);
  Distributed_Index_PlusPiPiPiPi = threadedOps->allocateDistibutedMemory(1024);
  destroyDistributedComplexVector(Distributed_Index_PlusPiPiPiPiXtrSize);
  Distributed_Index_PlusPiPiPiPiXtrSize = threadedOps->allocateDistibutedMemory(1024);  
  destroyDistributedComplexVector(Distributed_Index_PiModes);  
  Distributed_Index_PiModes = threadedOps->allocateDistibutedMemory(15);
}


void FermionMatrixOperations::setDiracParameters(double rho, double r) {
  if ((Parameter_rho != rho) || (Parameter_r != r)) {
    Parameter_rho = rho;
    Parameter_r = r;
    generateOpApplicationData();      
  }
}


void FermionMatrixOperations::getDiracParameters(double& rho, double& r) {
  rho = Parameter_rho;
  r = Parameter_r;
}


void FermionMatrixOperations::setPreconditioner(bool usePrecon, double m, double s) {
  m = fabs(m);
  s = fabs(s);
  if (m<0.01) m = 0.01;
  if ((Preconditioner_Usage != usePrecon) || (Preconditioner_M != m) || (Preconditioner_S != s)) {
    if ((LogLevel>1) && (usePrecon)) printf("Setting P - Preconditioning on with m= %1.3f and s= %1.3f\n",m,s);
    if ((LogLevel>1) && (!usePrecon)) printf("Setting P - Preconditioning off.\n");
    Preconditioner_Usage = usePrecon;  
    Preconditioner_M = m;
    Preconditioner_S = s;
    generateOpApplicationData();      
  }
}


void FermionMatrixOperations::setQPreconditioner(bool useqPrecon, double mu, double beta) {
  mu = fabs(mu);
  beta = fabs(beta);
  if (mu<0.01) mu = 0.01;
  if (beta<0.01) beta = 0.01;
  if ((QPreconditioner_Usage != useqPrecon) || (QPreconditioner_mu != mu) || (QPreconditioner_beta != beta)) {
    if ((LogLevel>1) && (useqPrecon)) printf("Setting Q - Preconditioning on with mu= %1.3f and beta= %1.3f\n",mu,beta);
    if ((LogLevel>1) && (!useqPrecon)) printf("Setting Q - Preconditioning off.\n");
    QPreconditioner_Usage = useqPrecon;  
    QPreconditioner_mu = mu;
    QPreconditioner_beta = beta;
    generateOpApplicationData();      
  }
}


void FermionMatrixOperations::setRPreconditioner(bool userPrecon, double m, double f) {
  m = fabs(m);
  f = fabs(f);
  if (m<0.01) m = 0.01;
  if (f<0.01) f = 0.01;
  if ((RPreconditioner_Usage != userPrecon) || (RPreconditioner_m != m) || (RPreconditioner_f != f)) {
    if ((LogLevel>1) && (userPrecon)) printf("Setting R - Preconditioning on with m= %1.3f and f= %1.3f\n",m,f);
    if ((LogLevel>1) && (!userPrecon)) printf("Setting R - Preconditioning off.\n");
    RPreconditioner_Usage = userPrecon;  
    RPreconditioner_m = m;
    RPreconditioner_f = f;
    generateOpApplicationData();      
  }
}


void FermionMatrixOperations::setMassSplitRatio(double split) {
  if (split <= 0) {
    printf("Mass split cannot be %e\n",split);
    exit(0);
  }
  if (split != Parameter_MassSplit) {
    Parameter_MassSplit = split;  
    generateOpApplicationData();      
    if (LogLevel>1) printf("Set Mass Split ratio to %f\n",Parameter_MassSplit);
  }
}


void FermionMatrixOperations::setExplicitMass(double mass) {
  if (mass != Parameter_ExplicitMass) {
    Parameter_ExplicitMass = mass;  
    generateOpApplicationData();      
    if (LogLevel>1) printf("Set Explicit mass to %f\n",Parameter_ExplicitMass);
  }
}


void FermionMatrixOperations::setAntiPeriodicBoundaryConditionsInTime(bool aPerBCinTime) {
  if (aPerBCinTime != Parameter_AntiPeriodicBoundaryConditionInTime) {
    Parameter_AntiPeriodicBoundaryConditionInTime = aPerBCinTime;  
    generateOpApplicationData();      
    if (LogLevel>1) printf("Set Anti-Periodic Boundary Conditions in time to %d\n",Parameter_AntiPeriodicBoundaryConditionInTime);
  }
}


double FermionMatrixOperations::getMassSplitRatio() {
  return Parameter_MassSplit;
}


double FermionMatrixOperations::getExplicitMass() {
  return Parameter_ExplicitMass;
}


void FermionMatrixOperations::setYukawaCoupling(double yN) {
  if (Parameter_yN != yN) {
    Parameter_yN = yN;
    generateOpApplicationData();      
  }
}


double FermionMatrixOperations::getYukawaCoupling() {
  return Parameter_yN;
}


void FermionMatrixOperations::clearTrafoPlans() {
  int I;
  
  for (I=0; I<FFTPlanDBMax; I++) {
    fftw_destroy_plan(FFTPlanDataBase[I]);
    InputVectorDataBase[I] = NULL;
    OutputVectorDataBase[I] = NULL;
  }
  FFTDataBaseCounter = 0;
}


void FermionMatrixOperations::desini() {
  if (LogLevel>0) printf("Desinitializing FermionMatrixOperations...");

  deactivateMultiThreadedOps();

  destroyFermionVector(InputVectorFourierTransform);
  destroyFermionVector(InterimVectorFourierTransform);  
  destroyFermionVector(Inverse_rest);
  destroyFermionVector(Inverse_b);
  destroyFermionVector(Inverse_p);
  destroyFermionVector(Inverse_s);  

  destroyFermionVector(QPreconditionerDiagonalData);
  QPreconditionerDiagonalData = NULL;
  destroyFermionVector(QPreconditionerDiagonalDataNormalized);
  QPreconditionerDiagonalDataNormalized = NULL;
  destroyFermionVector(inverseQPreconditionerDiagonalData);
  inverseQPreconditionerDiagonalData = NULL;
  destroyFermionVector(inverseQPreconditionerDiagonalDataNormalized);
  inverseQPreconditionerDiagonalDataNormalized = NULL;
  
  destroyFermionVector(RPreconditionerDiagonalData);
  RPreconditionerDiagonalData = NULL;
  destroyFermionVector(RPreconditionerDiagonalDataNormalized);
  RPreconditionerDiagonalDataNormalized = NULL;
  destroyFermionVector(inverseRPreconditionerDiagonalData);
  inverseRPreconditionerDiagonalData = NULL;
  destroyFermionVector(inverseRPreconditionerDiagonalDataNormalized);
  inverseRPreconditionerDiagonalDataNormalized = NULL;

  destroyFermionVector(RSQRPreconditionerDiagonalData);
  RSQRPreconditionerDiagonalData = NULL;
  destroyFermionVector(RSQRPreconditionerDiagonalDataNormalized);
  RSQRPreconditionerDiagonalDataNormalized = NULL;
  destroyFermionVector(inverseRSQRPreconditionerDiagonalData);
  inverseRSQRPreconditionerDiagonalData = NULL;
  destroyFermionVector(inverseRSQRPreconditionerDiagonalDataNormalized);
  inverseRSQRPreconditionerDiagonalDataNormalized = NULL;
  
  destroyFermionVector(NeubergerDiracOpApplicationEWData);
  NeubergerDiracOpApplicationEWData = NULL;
  destroyFermionVector(NeubergerDiracOpApplicationEWDataRescaled);
  NeubergerDiracOpApplicationEWDataRescaled = NULL;
  destroyFermionVector(NeubergerDiracPlusExplicitMassOpApplicationEWData);
  NeubergerDiracPlusExplicitMassOpApplicationEWData = NULL;
  destroyFermionVector(NeubergerDiracPlusExplicitMassOpApplicationEWDataRescaled);
  NeubergerDiracPlusExplicitMassOpApplicationEWDataRescaled = NULL;

  destroyFermionVector(NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWData);
  NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWData = NULL;
  destroyFermionVector(NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWDataRescaled);
  NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWDataRescaled = NULL;

  destroyFermionVector(NeubergerDiracOpApplicationSinPStdData);
  destroyFermionVector(NeubergerDiracOpApplicationSinPDagData);  
  NeubergerDiracOpApplicationSinPStdData = NULL;
  NeubergerDiracOpApplicationSinPDagData = NULL;
  destroyFermionVector(NeubergerWithChiStaticFermionInverseOpApplicationMainDiagonalData);
  destroyFermionVector(NeubergerWithChiStaticFermionInverseOpApplicationSubDiagonalData);
  NeubergerWithChiStaticFermionInverseOpApplicationMainDiagonalData = NULL;
  NeubergerWithChiStaticFermionInverseOpApplicationSubDiagonalData = NULL;
  destroyFermionVector(NeubergerWithChiStaticMdagMFermionInverseOpApplicationMainDiagonalData);
  destroyFermionVector(NeubergerWithChiStaticMdagMFermionInverseOpApplicationSubDiagonalData);
  NeubergerWithChiStaticMdagMFermionInverseOpApplicationMainDiagonalData = NULL;
  NeubergerWithChiStaticMdagMFermionInverseOpApplicationSubDiagonalData = NULL;
  destroyFermionVector(NeubergerWithChiStaticMdagMFermionInverseOpApplicationSSEAuxData);
  NeubergerWithChiStaticMdagMFermionInverseOpApplicationSSEAuxData = NULL;

  destroyFermionVector(halfMomentumForwardFFTfactors);
  halfMomentumForwardFFTfactors = NULL;
  destroyFermionVector(halfMomentumBackwardFFTfactors);
  halfMomentumBackwardFFTfactors = NULL;

  destroySuperAlignedLongInt(Index_PlusPiPiPiPi);
  Index_PlusPiPiPiPi = NULL;
  destroySuperAlignedLongInt(Index_PlusPiPiPiPiXtrSize);
  Index_PlusPiPiPiPiXtrSize = NULL;
  destroySuperAlignedLongInt(Index_PiModes);
  Index_PiModes = NULL;
printf("1\n");  

  delete xFFT;
  xFFT = NULL;
  
  delete NeubergerDiracOp;
  delete WilsonDiracOp;

  NeubergerDiracOp = NULL;
  WilsonDiracOp = NULL;
  clearTrafoPlans();
  if (LogLevel>0) printf("sucessfully.\n");      
}


FermionMatrixOperations::FermionMatrixOperations(int sizeL0, int sizeL1, int sizeL2, int sizeL3, double rho, double r, double yN) {
  if (LogLevel>1) printf("Initializing FermionMatrixOperations with L=%dx%dx%dx%d, rho=%f, r=%f, and yN=%f...\n", sizeL0, sizeL1, sizeL2, sizeL3, rho, r, yN);
  ini();

  setSize(sizeL0, sizeL1, sizeL2, sizeL3);  
  setDiracParameters(rho,r);
  setYukawaCoupling(yN);
}


FermionMatrixOperations::FermionMatrixOperations() {
  ini();
}


FermionMatrixOperations::~FermionMatrixOperations() {
  desini();
}


void FermionMatrixOperations::activateMultiThreadedOps(int OpMode, bool enforceCoreTie) {
  deactivateMultiThreadedOps();
  threadedOps = new MultiThreadedOperations(OpMode, enforceCoreTie);

  setSizeDistributed(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
  generateDistributedOpApplicationData();
}


void FermionMatrixOperations::deactivateMultiThreadedOps() {  
  destroyDistributedComplexVector(DistributedInputVectorFourierTransform);
  destroyDistributedComplexVector(DistributedInterimVectorFourierTransform);
  destroyDistributedComplexVector(DistributedInverse_rest);
  destroyDistributedComplexVector(DistributedInverse_b);
  destroyDistributedComplexVector(DistributedInverse_p);
  destroyDistributedComplexVector(DistributedInverse_s);

  destroyDistributedComplexVector(Distributed_NeubergerDiracOpApplicationEWData);
  destroyDistributedComplexVector(Distributed_NeubergerDiracOpApplicationEWDataRescaled);
  destroyDistributedComplexVector(Distributed_NeubergerDiracPlusExplicitMassOpApplicationEWData);
  destroyDistributedComplexVector(Distributed_NeubergerDiracPlusExplicitMassOpApplicationEWDataRescaled);
  destroyDistributedComplexVector(Distributed_NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWData);
  destroyDistributedComplexVector(Distributed_NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWDataRescaled);
  destroyDistributedComplexVector(Distributed_NeubergerDiracOpApplicationSinPStdData);
  destroyDistributedComplexVector(Distributed_NeubergerDiracOpApplicationSinPDagData);  
  destroyDistributedComplexVector(Distributed_NeubergerWithChiStaticFermionInverseOpApplicationMainDiagonalData);
  destroyDistributedComplexVector(Distributed_NeubergerWithChiStaticFermionInverseOpApplicationSubDiagonalData);
  destroyDistributedComplexVector(Distributed_NeubergerWithChiStaticMdagMFermionInverseOpApplicationMainDiagonalData);
  destroyDistributedComplexVector(Distributed_NeubergerWithChiStaticMdagMFermionInverseOpApplicationSubDiagonalData);
  destroyDistributedComplexVector(Distributed_NeubergerWithChiStaticMdagMFermionInverseOpApplicationSSEAuxData);
  destroyDistributedComplexVector(Distributed_QPreconditionerDiagonalData);
  destroyDistributedComplexVector(Distributed_QPreconditionerDiagonalDataNormalized);
  destroyDistributedComplexVector(Distributed_inverseQPreconditionerDiagonalData);
  destroyDistributedComplexVector(Distributed_inverseQPreconditionerDiagonalDataNormalized);
  destroyDistributedComplexVector(Distributed_RPreconditionerDiagonalData);
  destroyDistributedComplexVector(Distributed_RPreconditionerDiagonalDataNormalized);
  destroyDistributedComplexVector(Distributed_inverseRPreconditionerDiagonalData);
  destroyDistributedComplexVector(Distributed_inverseRPreconditionerDiagonalDataNormalized);  
  destroyDistributedComplexVector(Distributed_RSQRPreconditionerDiagonalData);
  destroyDistributedComplexVector(Distributed_RSQRPreconditionerDiagonalDataNormalized);
  destroyDistributedComplexVector(Distributed_inverseRSQRPreconditionerDiagonalData);
  destroyDistributedComplexVector(Distributed_inverseRSQRPreconditionerDiagonalDataNormalized);
  destroyDistributedComplexVector(Distributed_halfMomentumForwardFFTfactors);
  destroyDistributedComplexVector(Distributed_halfMomentumBackwardFFTfactors);
    
  destroyDistributedComplexVector(Distributed_Index_PlusPiPiPiPi);
  destroyDistributedComplexVector(Distributed_Index_PlusPiPiPiPiXtrSize);
  destroyDistributedComplexVector(Distributed_Index_PiModes);  

  delete threadedOps;
  threadedOps = NULL;
}


Complex* FermionMatrixOperations::createFermionVector(int localSize) {
  int size = localSize * (oneDimSizeL3 + xtraSize3) * (oneDimSizeL2 + xtraSize2) * (oneDimSizeL1 + xtraSize1) * oneDimSizeL0;
  return (Complex*) createSuperAlignedComplex(size, 128);
}


Complex* FermionMatrixOperations::createFermionVector() {
  return createFermionVector(8);
}


void FermionMatrixOperations::destroyFermionVector(Complex* &v) {
  destroySuperAlignedComplex(v);
  v = NULL;
}


void FermionMatrixOperations::zeroDistributedFermionVector(DistributedMemoryObject* v) {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::zeroDistributedFermionVector: threadedOps==NULL\n");
    exit(0);
  }
  threadedOps->zeroFermionVector(v, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
}


void FermionMatrixOperations::zeroFermionVector(Complex* v) {
  int I;
  for (I=0; I<vectorLengthXtrSize; I++) {
    v[I].x = 0;
    v[I].y = 0;
  }
}


DistributedMemoryObject* FermionMatrixOperations::createDistributedFermionVector() {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::createDistributedFermionVector: threadedOps==NULL\n");
    exit(0);
  }
  return threadedOps->allocateDistibutedFermionVector(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
}


DistributedMemoryObject* FermionMatrixOperations::createDistributedUniformLatticeComplexVector(int localSize) {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::createDistributedUniformLatticeComplexVector: threadedOps==NULL\n");
    exit(0);
  }
  return threadedOps->allocateDistibutedMemory(localSize * oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3);
}


void FermionMatrixOperations::destroyDistributedComplexVector(DistributedMemoryObject* &v) {
  delete v;
  v = NULL;
}


void FermionMatrixOperations::copyFermionVectorToDistributedFermionVector(Complex* vec, DistributedMemoryObject* memObj) {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::copyFermionVectorToDistributedFermionVector: threadedOps==NULL\n");
    exit(0);
  }
  threadedOps->copyFermionVectorToDistributedFermionVector(vec, memObj, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
}


void FermionMatrixOperations::copyDistributedFermionVectorToFermionVector(DistributedMemoryObject* memObj, Complex* vec) {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::copyDistributedFermionVectorToFermionVector: threadedOps==NULL\n");
    exit(0);
  }
  threadedOps->copyDistributedFermionVectorToFermionVector(memObj, vec, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
}


void FermionMatrixOperations::copyPhiFieldUniformlyToDistributedPhiField(double* phiField, DistributedMemoryObject* memObj) {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::copyPhiFieldUniformlyToDistributedPhiField: threadedOps==NULL\n");
    exit(0);
  }
  threadedOps->copyComplexVectorToUniformlyDistributedMemory((Complex*)phiField, memObj, 2*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3);
}


bool FermionMatrixOperations::isMultiThreadedOpsActivated() {
  return (threadedOps!=NULL);
}


MultiThreadedOperations* FermionMatrixOperations::getMultiThreadedOps() {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::getsMultiThreadedOps: threadedOps==NULL\n");
    exit(0);
  }
  return threadedOps;
}


void FermionMatrixOperations::planeWaveFermionVector(Complex* v, double p0, double p1, double p2, double p3) {
  double Arg0 = 2*p0*pi/oneDimSizeL0;
  double Arg1 = 2*p1*pi/oneDimSizeL1;
  double Arg2 = 2*p2*pi/oneDimSizeL2;
  double Arg3 = 2*p3*pi/oneDimSizeL3;
  double norm = 1.0/sqrt(oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3);

  int count = 0;
  for (int I0=0; I0<oneDimSizeL0; I0++) {
    double c0 = cos(Arg0*I0);
    for (int I1=0; I1<oneDimSizeL1; I1++) {
      double c1 = cos(Arg1*I1);
      for (int I2=0; I2<oneDimSizeL2; I2++) {
        double c2 = cos(Arg2*I2);
        for (int I3=0; I3<oneDimSizeL3; I3++) {
          double c3 = cos(Arg3*I3);
          for (int i=0; i<8; i++) {
	    v[i+count].x = 0;
	    v[i+count].y = 0;
	  }
	  v[count].x = norm*c0*c1*c2*c3;
	  count += 8;
        }
      }
    }
  }
  transformToXtraSizeArray(v, v);
}


int FermionMatrixOperations::getVectorLength() {
  return vectorLength;
}


int FermionMatrixOperations::getVectorLengthXtrSize() {
  return vectorLengthXtrSize;
}


int FermionMatrixOperations::get1DSizeL0() {
  return oneDimSizeL0;
}


int FermionMatrixOperations::get1DSizeL1() {
  return oneDimSizeL1;
}


int FermionMatrixOperations::get1DSizeL2() {
  return oneDimSizeL2;
}


int FermionMatrixOperations::get1DSizeL3() {
  return oneDimSizeL3;
}


int FermionMatrixOperations::get1DSizeLargest() {
  return oneDimSizeLargest;
}


int FermionMatrixOperations::getTimeDirection() {
  return timeDirection;
}


double FermionMatrixOperations::getYN() {
  return Parameter_yN;
}


void FermionMatrixOperations::generateOpApplicationData() {
  int I0,I1,I2,I3;
  vector4D p, pPi;
  int count = 0;
  Complex nup;
  Complex nupPi;
  double twoRho = 2*Parameter_rho;
  Complex m00, m10, m01, m11, Gamma;
  double fac = 1.0/(oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3);
  double m = Preconditioner_M;
  double s = Preconditioner_S;
  int N0half = oneDimSizeL0 / 2;
  int N1half = oneDimSizeL1 / 2;
  int N2half = oneDimSizeL2 / 2;
  int N3half = oneDimSizeL3 / 2;


  for (I0=0; I0<oneDimSizeL0; I0++) {
    p[0] = 2*I0*pi/oneDimSizeL0;
    if ((Parameter_AntiPeriodicBoundaryConditionInTime) && (timeDirection==0)) p[0] += pi/oneDimSizeL0;

    NeubergerDiracOpApplicationSinPStdData[0*2*128+2*(oneDimSizeL0-I0-1)+0].x = +sin(p[0]);
    NeubergerDiracOpApplicationSinPStdData[0*2*128+2*(oneDimSizeL0-I0-1)+0].y = +sin(p[0]);
    NeubergerDiracOpApplicationSinPStdData[0*2*128+2*(oneDimSizeL0-I0-1)+1].x = -sin(p[0]);
    NeubergerDiracOpApplicationSinPStdData[0*2*128+2*(oneDimSizeL0-I0-1)+1].y = +sin(p[0]);
    
    NeubergerDiracOpApplicationSinPDagData[0*2*128+2*(oneDimSizeL0-I0-1)+0].x = -sin(p[0]);
    NeubergerDiracOpApplicationSinPDagData[0*2*128+2*(oneDimSizeL0-I0-1)+0].y = -sin(p[0]);
    NeubergerDiracOpApplicationSinPDagData[0*2*128+2*(oneDimSizeL0-I0-1)+1].x = +sin(p[0]);
    NeubergerDiracOpApplicationSinPDagData[0*2*128+2*(oneDimSizeL0-I0-1)+1].y = -sin(p[0]);
  }
  for (I1=0; I1<oneDimSizeL1; I1++) {
    p[1] = 2*I1*pi/oneDimSizeL1;
    if ((Parameter_AntiPeriodicBoundaryConditionInTime) && (timeDirection==1)) p[1] += pi/oneDimSizeL1;

    NeubergerDiracOpApplicationSinPStdData[1*2*128+2*(oneDimSizeL1-I1-1)+0].x = +sin(p[1]);
    NeubergerDiracOpApplicationSinPStdData[1*2*128+2*(oneDimSizeL1-I1-1)+0].y = +sin(p[1]);
    NeubergerDiracOpApplicationSinPStdData[1*2*128+2*(oneDimSizeL1-I1-1)+1].x = -sin(p[1]);
    NeubergerDiracOpApplicationSinPStdData[1*2*128+2*(oneDimSizeL1-I1-1)+1].y = +sin(p[1]);
    
    NeubergerDiracOpApplicationSinPDagData[1*2*128+2*(oneDimSizeL1-I1-1)+0].x = -sin(p[1]);
    NeubergerDiracOpApplicationSinPDagData[1*2*128+2*(oneDimSizeL1-I1-1)+0].y = -sin(p[1]);
    NeubergerDiracOpApplicationSinPDagData[1*2*128+2*(oneDimSizeL1-I1-1)+1].x = +sin(p[1]);
    NeubergerDiracOpApplicationSinPDagData[1*2*128+2*(oneDimSizeL1-I1-1)+1].y = -sin(p[1]);
  }
  for (I2=0; I2<oneDimSizeL2; I2++) {
    p[2] = 2*I2*pi/oneDimSizeL2;
    if ((Parameter_AntiPeriodicBoundaryConditionInTime) && (timeDirection==2)) p[2] += pi/oneDimSizeL2;

    NeubergerDiracOpApplicationSinPStdData[2*2*128+2*(oneDimSizeL2-I2-1)+0].x = +sin(p[2]);
    NeubergerDiracOpApplicationSinPStdData[2*2*128+2*(oneDimSizeL2-I2-1)+0].y = +sin(p[2]);
    NeubergerDiracOpApplicationSinPStdData[2*2*128+2*(oneDimSizeL2-I2-1)+1].x = -sin(p[2]);
    NeubergerDiracOpApplicationSinPStdData[2*2*128+2*(oneDimSizeL2-I2-1)+1].y = +sin(p[2]);
   
    NeubergerDiracOpApplicationSinPDagData[2*2*128+2*(oneDimSizeL2-I2-1)+0].x = -sin(p[2]);
    NeubergerDiracOpApplicationSinPDagData[2*2*128+2*(oneDimSizeL2-I2-1)+0].y = -sin(p[2]);
    NeubergerDiracOpApplicationSinPDagData[2*2*128+2*(oneDimSizeL2-I2-1)+1].x = +sin(p[2]);
    NeubergerDiracOpApplicationSinPDagData[2*2*128+2*(oneDimSizeL2-I2-1)+1].y = -sin(p[2]);
  }
  for (I3=0; I3<oneDimSizeL3; I3++) {
    p[3] = 2*I3*pi/oneDimSizeL3;
    if ((Parameter_AntiPeriodicBoundaryConditionInTime) && (timeDirection==3)) p[3] += pi/oneDimSizeL3;

    NeubergerDiracOpApplicationSinPStdData[3*2*128+2*(oneDimSizeL3-I3-1)+0].x = +sin(p[3]);
    NeubergerDiracOpApplicationSinPStdData[3*2*128+2*(oneDimSizeL3-I3-1)+0].y = +sin(p[3]);
    NeubergerDiracOpApplicationSinPStdData[3*2*128+2*(oneDimSizeL3-I3-1)+1].x = -sin(p[3]);
    NeubergerDiracOpApplicationSinPStdData[3*2*128+2*(oneDimSizeL3-I3-1)+1].y = +sin(p[3]);
    
    NeubergerDiracOpApplicationSinPDagData[3*2*128+2*(oneDimSizeL3-I3-1)+0].x = -sin(p[3]);
    NeubergerDiracOpApplicationSinPDagData[3*2*128+2*(oneDimSizeL3-I3-1)+0].y = -sin(p[3]);
    NeubergerDiracOpApplicationSinPDagData[3*2*128+2*(oneDimSizeL3-I3-1)+1].x = +sin(p[3]);
    NeubergerDiracOpApplicationSinPDagData[3*2*128+2*(oneDimSizeL3-I3-1)+1].y = -sin(p[3]);
  }
  
  for (I0=0; I0<oneDimSizeLargest; I0++) {
    halfMomentumForwardFFTfactors[I0] = exp((-I0*pi/oneDimSizeLargest)*ComplexI);
    halfMomentumBackwardFFTfactors[I0] = exp((I0*pi/oneDimSizeLargest)*ComplexI);
  }
 
  for (I0=0; I0<oneDimSizeL0; I0++) {
    p[0] = 2*I0*pi/oneDimSizeL0;
    if ((Parameter_AntiPeriodicBoundaryConditionInTime) && (timeDirection==0)) p[0] += pi/oneDimSizeL0;

    pPi[0] = p[0] + pi;
    for (I1=0; I1<oneDimSizeL1; I1++) {
      p[1] = 2*I1*pi/oneDimSizeL1;
      if ((Parameter_AntiPeriodicBoundaryConditionInTime) && (timeDirection==1)) p[1] += pi/oneDimSizeL1;

      pPi[1] = p[1] + pi;
      for (I2=0; I2<oneDimSizeL2; I2++) {
        p[2] = 2*I2*pi/oneDimSizeL2;
        if ((Parameter_AntiPeriodicBoundaryConditionInTime) && (timeDirection==2)) p[2] += pi/oneDimSizeL2;

        pPi[2] = p[2] + pi;
        for (I3=0; I3<oneDimSizeL3; I3++) {
          p[3] = 2*I3*pi/oneDimSizeL3;
          if ((Parameter_AntiPeriodicBoundaryConditionInTime) && (timeDirection==3)) p[3] += pi/oneDimSizeL3;
          pPi[3] = p[3] + pi;

          nup = NeubergerDiracOp->analyticalEigenvalue(p);
          nupPi = NeubergerDiracOp->analyticalEigenvalue(pPi);
	  
	  double psqr = 4*(sqr(sin(p[0]/2)) + sqr(sin(p[1]/2)) + sqr(sin(p[2]/2)) + sqr(sin(p[3]/2)));
	  if ((!isNaN(QPreconditioner_beta)) && (!isNaN(QPreconditioner_mu))) {
            double Qd1 = exp(log(16+1E-10)*QPreconditioner_beta);
	    double Qd2 = exp(log(sqr(QPreconditioner_mu))*QPreconditioner_beta);
	    double Qd3 = exp(log(psqr+1E-10)*QPreconditioner_beta);
	    QPreconditionerDiagonalData[count].x = (Qd1 + Qd2) / (Qd3 + Qd2);
            QPreconditionerDiagonalData[count].x -= 0.15*psqr/16.0;
	    QPreconditionerDiagonalData[count].y = 0;

            QPreconditionerDiagonalDataNormalized[count].x = fac * QPreconditionerDiagonalData[count].x;
            QPreconditionerDiagonalDataNormalized[count].y = 0;
            inverseQPreconditionerDiagonalData[count].x = 1.0 / QPreconditionerDiagonalData[count].x;
            inverseQPreconditionerDiagonalData[count].y = 0;
            inverseQPreconditionerDiagonalDataNormalized[count].y = fac / QPreconditionerDiagonalData[count].x;
            inverseQPreconditionerDiagonalDataNormalized[count].y = 0;
	  }

	  
          NeubergerDiracOpApplicationEWData[count] = fac * nup;
          NeubergerDiracOpApplicationEWDataRescaled[count] = nup;
          NeubergerDiracPlusExplicitMassOpApplicationEWData[count] = fac * (nup + Parameter_ExplicitMass);
          NeubergerDiracPlusExplicitMassOpApplicationEWDataRescaled[count] = nup + Parameter_ExplicitMass;

          Complex ewp(1.0,0);
          if (Parameter_ExplicitMass != 0) {
            if (fabs(nup.x-twoRho)>10000*eps) {
              ewp = ewp - ((1/twoRho)*nup);
              ewp = (nup+Parameter_ExplicitMass) / ewp;
            } else {
              ewp = nup+Parameter_ExplicitMass;
	    }
          } else {
            if (fabs(nup.x-twoRho)>10000*eps) {
              ewp = ewp - ((1/twoRho)*nup);
              ewp = nup / ewp;
            } else {
              ewp = nup;
	    }
          }

          NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWData[count] = fac*ewp;
          NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWDataRescaled[count] = ewp;


	  if ((!isNaN(RPreconditioner_m)) && (!isNaN(RPreconditioner_f))) {
	    Complex e1 = ewp;
	    Complex e2 = ewp;
            if (fabs(nup.x-twoRho)>10000*eps) {
              e1.x += Parameter_yN * RPreconditioner_m;
              e2.x += Parameter_MassSplit * Parameter_yN * RPreconditioner_m;	      
	    }
            e1.x = 1.0/sqrt(sqrt(e1.x*e1.x + e1.y*e1.y) / RPreconditioner_f);
            e1.y = 0;
            e2.x = 1.0/sqrt(sqrt(e2.x*e2.x + e2.y*e2.y) / RPreconditioner_f);
            e2.y = 0;

            RPreconditionerDiagonalData[count].x = e1.x;
            RPreconditionerDiagonalData[count].y = e2.x;    
            RPreconditionerDiagonalDataNormalized[count].x = fac * RPreconditionerDiagonalData[count].x;
            RPreconditionerDiagonalDataNormalized[count].y = fac * RPreconditionerDiagonalData[count].y;	    
            inverseRPreconditionerDiagonalData[count].x = 1.0 / RPreconditionerDiagonalData[count].x;
            inverseRPreconditionerDiagonalData[count].y = 1.0 / RPreconditionerDiagonalData[count].y;	    
            inverseRPreconditionerDiagonalDataNormalized[count].x = fac / RPreconditionerDiagonalData[count].x; 
            inverseRPreconditionerDiagonalDataNormalized[count].y = fac / RPreconditionerDiagonalData[count].y;
	    
	    e1.x = e1.x * e1.x;
            e1.y = 0;
	    e2.x = e2.x * e2.x;
            e2.y = 0;
	    
            RSQRPreconditionerDiagonalData[count].x = e1.x;
            RSQRPreconditionerDiagonalData[count].y = e2.x;    
            RSQRPreconditionerDiagonalDataNormalized[count].x = fac * RPreconditionerDiagonalData[count].x;
            RSQRPreconditionerDiagonalDataNormalized[count].y = fac * RPreconditionerDiagonalData[count].y;	    
            inverseRSQRPreconditionerDiagonalData[count].x = 1.0 / RPreconditionerDiagonalData[count].x;
            inverseRSQRPreconditionerDiagonalData[count].y = 1.0 / RPreconditionerDiagonalData[count].y;   
            inverseRSQRPreconditionerDiagonalDataNormalized[count].x = fac / RPreconditionerDiagonalData[count].x; 
            inverseRSQRPreconditionerDiagonalDataNormalized[count].y = fac / RPreconditionerDiagonalData[count].y;
          }


          if (Parameter_ExplicitMass != 0) {
            m00 = Parameter_yN*m*(nup - twoRho) - twoRho*(nup+Parameter_ExplicitMass);
            m01 = Parameter_yN*s*(adj(nupPi) - twoRho);
            m10 = Parameter_yN*s*(nup - twoRho);
            m11 = Parameter_yN*m*(adj(nupPi) - twoRho) - twoRho*adj(nupPi+Parameter_ExplicitMass);
          } else {
            m00 = Parameter_yN*m*(nup - twoRho) - twoRho*nup;
            m01 = Parameter_yN*s*(adj(nupPi) - twoRho);
            m10 = Parameter_yN*s*(nup - twoRho);
            m11 = Parameter_yN*m*(adj(nupPi) - twoRho) - twoRho*adj(nupPi);
          }

          Gamma = Complex(fac,0) / (m00*m11 - m10*m01);

	  NeubergerWithChiStaticFermionInverseOpApplicationMainDiagonalData[count] = Gamma*m11;
	  NeubergerWithChiStaticFermionInverseOpApplicationSubDiagonalData[count] = (-1.0)*Gamma*m10;

          if (Parameter_ExplicitMass != 0) {
            m00.x = sqr(Parameter_yN)*(m*m+s*s)*normSqr(nup-twoRho) - 2*twoRho*Parameter_yN*m*((adj(nup+Parameter_ExplicitMass)*(nup-twoRho)).x) + sqr(twoRho)*normSqr(nup+Parameter_ExplicitMass);
            m00.y = 0;
            m01 = 2*m*s*sqr(Parameter_yN)*adj(nup - twoRho)*adj(nupPi-twoRho) - twoRho*Parameter_yN*s*adj(nup+Parameter_ExplicitMass)*adj(nupPi-twoRho) - twoRho*Parameter_yN*s*adj(nup-twoRho)*adj(nupPi+Parameter_ExplicitMass);
            m10 = 2*m*s*sqr(Parameter_yN)*(nupPi - twoRho)*(nup-twoRho) - twoRho*Parameter_yN*s*(nupPi+Parameter_ExplicitMass)*(nup-twoRho) - twoRho*Parameter_yN*s*(nupPi-twoRho)*(nup+Parameter_ExplicitMass);
            m11.x = sqr(Parameter_yN)*(m*m+s*s)*normSqr(nupPi-twoRho) - 2*twoRho*Parameter_yN*m*((adj(nupPi+Parameter_ExplicitMass)*(nupPi-twoRho)).x) + sqr(twoRho)*normSqr(nupPi+Parameter_ExplicitMass);
            m11.y = 0;
          } else {
            m00.x = sqr(Parameter_yN)*(m*m+s*s)*normSqr(nup-twoRho) - 2*twoRho*Parameter_yN*m*((adj(nup)*(nup-twoRho)).x) + sqr(twoRho)*normSqr(nup);
            m00.y = 0;
            m01 = 2*m*s*sqr(Parameter_yN)*adj(nup - twoRho)*adj(nupPi-twoRho) - twoRho*Parameter_yN*s*adj(nup)*adj(nupPi-twoRho) - twoRho*Parameter_yN*s*adj(nup-twoRho)*adj(nupPi);
            m10 = 2*m*s*sqr(Parameter_yN)*(nupPi - twoRho)*(nup-twoRho) - twoRho*Parameter_yN*s*(nupPi)*(nup-twoRho) - twoRho*Parameter_yN*s*(nupPi-twoRho)*(nup);
            m11.x = sqr(Parameter_yN)*(m*m+s*s)*normSqr(nupPi-twoRho) - 2*twoRho*Parameter_yN*m*((adj(nupPi)*(nupPi-twoRho)).x) + sqr(twoRho)*normSqr(nupPi);
            m11.y = 0;
          }
          Gamma = Complex(fac,0) / (m00*m11 - m10*m01);
	  NeubergerWithChiStaticMdagMFermionInverseOpApplicationMainDiagonalData[count] = Gamma*m11;
	  NeubergerWithChiStaticMdagMFermionInverseOpApplicationSubDiagonalData[count] = (-1.0)*Gamma*m10;

	  double normP = sqr(NeubergerDiracOpApplicationSinPStdData[0*2*128+2*(oneDimSizeL0-I0-1)+0].x);
	  normP += sqr(NeubergerDiracOpApplicationSinPStdData[1*2*128+2*(oneDimSizeL1-I1-1)+0].x);
	  normP += sqr(NeubergerDiracOpApplicationSinPStdData[2*2*128+2*(oneDimSizeL2-I2-1)+0].x);
	  normP += sqr(NeubergerDiracOpApplicationSinPStdData[3*2*128+2*(oneDimSizeL3-I3-1)+0].x);
	  normP = sqrt(normP);

	  
	  if (normP>100*eps) {
  	    NeubergerDiracOpApplicationEWData[count].y /= normP;
  	    NeubergerDiracOpApplicationEWDataRescaled[count].y /= normP;
  	    NeubergerDiracPlusExplicitMassOpApplicationEWData[count].y /= normP;
  	    NeubergerDiracPlusExplicitMassOpApplicationEWDataRescaled[count].y /= normP;
            NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWData[count].y /= normP;
            NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWDataRescaled[count].y /= normP;
	    NeubergerWithChiStaticFermionInverseOpApplicationMainDiagonalData[count].y /= normP;
	    NeubergerWithChiStaticFermionInverseOpApplicationSubDiagonalData[count].y /= normP;
	    NeubergerWithChiStaticMdagMFermionInverseOpApplicationMainDiagonalData[count].y /= normP;
	    NeubergerWithChiStaticMdagMFermionInverseOpApplicationSubDiagonalData[count].y /= normP;
	  } else {
  	    NeubergerDiracOpApplicationEWData[count].y = 0;
  	    NeubergerDiracOpApplicationEWDataRescaled[count].y = 0;
  	    NeubergerDiracPlusExplicitMassOpApplicationEWData[count].y = 0;
  	    NeubergerDiracPlusExplicitMassOpApplicationEWDataRescaled[count].y = 0;
            NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWData[count].y = 0;
            NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWDataRescaled[count].y = 0;
	    NeubergerWithChiStaticFermionInverseOpApplicationMainDiagonalData[count].y = 0;
	    NeubergerWithChiStaticFermionInverseOpApplicationSubDiagonalData[count].y = 0;
	    NeubergerWithChiStaticMdagMFermionInverseOpApplicationMainDiagonalData[count].y = 0;
	    NeubergerWithChiStaticMdagMFermionInverseOpApplicationSubDiagonalData[count].y = 0;
	  }
          count++;
	}
      }
    }
  }

  for (I0=0; I0<oneDimSizeL0; I0++) Index_PlusPiPiPiPi[4*(oneDimSizeL0-1-I0+0)]   = 8*16*((I0+N0half)%oneDimSizeL0)*oneDimSizeL3*oneDimSizeL2 *oneDimSizeL1;
  for (I1=0; I1<oneDimSizeL1; I1++) Index_PlusPiPiPiPi[4*(oneDimSizeL1-1-I1+128)] = 8*16*((I1+N1half)%oneDimSizeL1)*oneDimSizeL3*oneDimSizeL2;
  for (I2=0; I2<oneDimSizeL2; I2++) Index_PlusPiPiPiPi[4*(oneDimSizeL2-1-I2+256)] = 8*16*((I2+N2half)%oneDimSizeL2)*oneDimSizeL3;
  for (I3=0; I3<oneDimSizeL3; I3++) Index_PlusPiPiPiPi[4*(oneDimSizeL3-1-I3+384)] = 8*16*((I3+N3half)%oneDimSizeL3);

  for (I0=0; I0<oneDimSizeL0; I0++) Index_PlusPiPiPiPiXtrSize[4*(oneDimSizeL0-1-I0+0)]   = 8*16*((I0+N0half)%oneDimSizeL0)*(oneDimSizeL3+xtraSize3)*(oneDimSizeL2+xtraSize2) *(oneDimSizeL1+xtraSize1);
  for (I1=0; I1<oneDimSizeL1; I1++) Index_PlusPiPiPiPiXtrSize[4*(oneDimSizeL1-1-I1+128)] = 8*16*((I1+N1half)%oneDimSizeL1)*(oneDimSizeL3+xtraSize3)*(oneDimSizeL2+xtraSize2);
  for (I2=0; I2<oneDimSizeL2; I2++) Index_PlusPiPiPiPiXtrSize[4*(oneDimSizeL2-1-I2+256)] = 8*16*((I2+N2half)%oneDimSizeL2)*(oneDimSizeL3+xtraSize3);
  for (I3=0; I3<oneDimSizeL3; I3++) Index_PlusPiPiPiPiXtrSize[4*(oneDimSizeL3-1-I3+384)] = 8*16*((I3+N3half)%oneDimSizeL3);
  
  count = 0;
  int countPi = 0;
  for (I0=oneDimSizeL0-1; I0>=N0half; I0--) {
    countPi += Index_PlusPiPiPiPi[4*(I0+0)];
    for (I1=oneDimSizeL1-1; I1>=0; I1--) {
      countPi += Index_PlusPiPiPiPi[4*(I1+128)];
      for (I2=oneDimSizeL2-1; I2>=0; I2--) {
        countPi += Index_PlusPiPiPiPi[4*(I2+256)];
        for (I3=oneDimSizeL3-1; I3>=0; I3--) {
          countPi += Index_PlusPiPiPiPi[4*(I3+384)];

          NeubergerWithChiStaticMdagMFermionInverseOpApplicationSSEAuxData[4*count+0].x = NeubergerWithChiStaticMdagMFermionInverseOpApplicationMainDiagonalData[count].x;	  
          NeubergerWithChiStaticMdagMFermionInverseOpApplicationSSEAuxData[4*count+0].y = NeubergerWithChiStaticMdagMFermionInverseOpApplicationMainDiagonalData[count].x;	  
          NeubergerWithChiStaticMdagMFermionInverseOpApplicationSSEAuxData[4*count+1]   = NeubergerWithChiStaticMdagMFermionInverseOpApplicationSubDiagonalData[count];
          NeubergerWithChiStaticMdagMFermionInverseOpApplicationSSEAuxData[4*count+2].x = NeubergerWithChiStaticMdagMFermionInverseOpApplicationMainDiagonalData[countPi/128].x;	  
          NeubergerWithChiStaticMdagMFermionInverseOpApplicationSSEAuxData[4*count+2].y = NeubergerWithChiStaticMdagMFermionInverseOpApplicationMainDiagonalData[countPi/128].x;	  
          NeubergerWithChiStaticMdagMFermionInverseOpApplicationSSEAuxData[4*count+3]   = adj(NeubergerWithChiStaticMdagMFermionInverseOpApplicationSubDiagonalData[countPi/128]);
	  
	  count++;
          countPi -= Index_PlusPiPiPiPi[4*(I3+384)];
	}
        countPi -= Index_PlusPiPiPiPi[4*(I2+256)];
      }
      countPi -= Index_PlusPiPiPiPi[4*(I1+128)];
    }
    countPi -= Index_PlusPiPiPiPi[4*(I0+0)];
  }

  Index_PiModes[0] = 8*(oneDimSizeL3/2);
  Index_PiModes[1] = 8*(oneDimSizeL2/2)*(oneDimSizeL3+xtraSize3);
  Index_PiModes[2] = Index_PiModes[1]+Index_PiModes[0];
  Index_PiModes[3] = 8*(oneDimSizeL1/2)*(oneDimSizeL2+xtraSize2)*(oneDimSizeL3+xtraSize3);
  Index_PiModes[4] = Index_PiModes[3]+Index_PiModes[0];
  Index_PiModes[5] = Index_PiModes[3]+Index_PiModes[1];
  Index_PiModes[6] = Index_PiModes[3]+Index_PiModes[1]+Index_PiModes[0];
  Index_PiModes[7] = 8*(oneDimSizeL0/2)*(oneDimSizeL1+xtraSize1)*(oneDimSizeL2+xtraSize2)*(oneDimSizeL3+xtraSize3);
  Index_PiModes[8] = Index_PiModes[7]+Index_PiModes[0];
  Index_PiModes[9] = Index_PiModes[7]+Index_PiModes[1];
  Index_PiModes[10] = Index_PiModes[7]+Index_PiModes[1]+Index_PiModes[0];
  Index_PiModes[11] = Index_PiModes[7]+Index_PiModes[3];
  Index_PiModes[12] = Index_PiModes[7]+Index_PiModes[3]+Index_PiModes[0];
  Index_PiModes[13] = Index_PiModes[7]+Index_PiModes[3]+Index_PiModes[1];
  Index_PiModes[14] = Index_PiModes[7]+Index_PiModes[3]+Index_PiModes[1]+Index_PiModes[0];
  
  if (Parameter_AntiPeriodicBoundaryConditionInTime) {
    for (int I=0; I<15; I++) {
      Index_PiModes[I] = -100;  //No Pi-Modes!!!
    }
  }

  generateDistributedOpApplicationData();
}


void FermionMatrixOperations::generateDistributedOpApplicationData() {
  if (threadedOps==NULL) return;

  threadedOps->copyComplexVectorToUniformlyDistributedMemory(NeubergerDiracOpApplicationEWData, Distributed_NeubergerDiracOpApplicationEWData, 1*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3); 
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(NeubergerDiracOpApplicationEWDataRescaled, Distributed_NeubergerDiracOpApplicationEWDataRescaled, 1*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3);
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(NeubergerDiracPlusExplicitMassOpApplicationEWData, Distributed_NeubergerDiracPlusExplicitMassOpApplicationEWData, 1*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3); 
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(NeubergerDiracPlusExplicitMassOpApplicationEWDataRescaled, Distributed_NeubergerDiracPlusExplicitMassOpApplicationEWDataRescaled, 1*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3);
 
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWData, Distributed_NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWData, 1*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3); 
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWDataRescaled, Distributed_NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWDataRescaled, 1*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3); 

  threadedOps->copyComplexVectorToUniformlyDistributedMemory(NeubergerDiracOpApplicationSinPStdData, Distributed_NeubergerDiracOpApplicationSinPStdData, 2*4*128); 
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(NeubergerDiracOpApplicationSinPDagData, Distributed_NeubergerDiracOpApplicationSinPDagData, 2*4*128); 

  threadedOps->copyComplexVectorToUniformlyDistributedMemory(NeubergerWithChiStaticFermionInverseOpApplicationMainDiagonalData, Distributed_NeubergerWithChiStaticFermionInverseOpApplicationMainDiagonalData, 1*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3);  
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(NeubergerWithChiStaticFermionInverseOpApplicationSubDiagonalData, Distributed_NeubergerWithChiStaticFermionInverseOpApplicationSubDiagonalData, 1*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3); 
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(NeubergerWithChiStaticMdagMFermionInverseOpApplicationMainDiagonalData, Distributed_NeubergerWithChiStaticMdagMFermionInverseOpApplicationMainDiagonalData, 1*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3); 
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(NeubergerWithChiStaticMdagMFermionInverseOpApplicationSubDiagonalData, Distributed_NeubergerWithChiStaticMdagMFermionInverseOpApplicationSubDiagonalData, 1*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3); 

  threadedOps->copyComplexVectorToUniformlyDistributedMemory(NeubergerWithChiStaticMdagMFermionInverseOpApplicationSSEAuxData, Distributed_NeubergerWithChiStaticMdagMFermionInverseOpApplicationSSEAuxData, 2*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3); 
  
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(QPreconditionerDiagonalData, Distributed_QPreconditionerDiagonalData, 1*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3); 
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(QPreconditionerDiagonalDataNormalized, Distributed_QPreconditionerDiagonalDataNormalized, 1*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3); 
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(inverseQPreconditionerDiagonalData, Distributed_inverseQPreconditionerDiagonalData, 1*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3); 
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(inverseQPreconditionerDiagonalDataNormalized, Distributed_inverseQPreconditionerDiagonalDataNormalized, 1*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3); 
  
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(RPreconditionerDiagonalData, Distributed_RPreconditionerDiagonalData, 1*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3); 
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(RPreconditionerDiagonalDataNormalized, Distributed_RPreconditionerDiagonalDataNormalized, 1*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3); 
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(inverseRPreconditionerDiagonalData, Distributed_inverseRPreconditionerDiagonalData, 1*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3); 
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(inverseRPreconditionerDiagonalDataNormalized, Distributed_inverseRPreconditionerDiagonalDataNormalized, 1*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3); 
  
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(RSQRPreconditionerDiagonalData, Distributed_RSQRPreconditionerDiagonalData, 1*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3); 
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(RSQRPreconditionerDiagonalDataNormalized, Distributed_RSQRPreconditionerDiagonalDataNormalized, 1*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3); 
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(inverseRSQRPreconditionerDiagonalData, Distributed_inverseRSQRPreconditionerDiagonalData, 1*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3); 
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(inverseRSQRPreconditionerDiagonalDataNormalized,Distributed_inverseRSQRPreconditionerDiagonalDataNormalized, 1*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3); 
  
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(halfMomentumForwardFFTfactors, Distributed_halfMomentumForwardFFTfactors, oneDimSizeL0+oneDimSizeL1+oneDimSizeL2+oneDimSizeL3); 
  threadedOps->copyComplexVectorToUniformlyDistributedMemory(halfMomentumBackwardFFTfactors, Distributed_halfMomentumBackwardFFTfactors, oneDimSizeL0+oneDimSizeL1+oneDimSizeL2+oneDimSizeL3); 

  threadedOps->copyComplexVectorToUniformlyDistributedMemory((Complex*)Index_PlusPiPiPiPi, Distributed_Index_PlusPiPiPiPi, 1024); 
  threadedOps->copyComplexVectorToUniformlyDistributedMemory((Complex*)Index_PlusPiPiPiPiXtrSize, Distributed_Index_PlusPiPiPiPiXtrSize, 1024); 
  threadedOps->copyComplexVectorToUniformlyDistributedMemory((Complex*)Index_PiModes, Distributed_Index_PiModes, 15); 
}


void FermionMatrixOperations::fillGaussRandomVector(Complex* v, int count) {
  int I;
  
  if (count<0) count = vectorLength;
  for (I=0; I<count; I++) {
    AdvancedGaussZufall(AdvancedSeed, v[I].x, v[I].y);
  }
}


void FermionMatrixOperations::multiplyFermionVectorWithTimeIndexedComplexScalars(Complex* input, Complex* output, Complex* scalars) {
  int ind[4];
  int xtrAdd1 = xtraSize1*8*(oneDimSizeL3+xtraSize3)*(oneDimSizeL2+xtraSize2);
  int xtrAdd2 = xtraSize2*8*(oneDimSizeL3+xtraSize3);
  int xtrAdd3 = xtraSize3*8;
  int index = 0;
  for (ind[0]=0; ind[0]<oneDimSizeL0; ind[0]++) {
    for (ind[1]=0; ind[1]<oneDimSizeL1; ind[1]++) {
      for (ind[2]=0; ind[2]<oneDimSizeL2; ind[2]++) {
        for (ind[3]=0; ind[3]<oneDimSizeL3; ind[3]++) {
          double fx = scalars[ind[timeDirection]].x;
          double fy = scalars[ind[timeDirection]].y;
          for (int i=0; i<8; i++) {
            double dummy = input[index].x*fx - input[index].y*fy;
            output[index].y = input[index].x*fy + input[index].y*fx;
            output[index].x = dummy; 
            index++;
          }
        }
        index += xtrAdd3; 
      }
      index += xtrAdd2;
    }
    index += xtrAdd1;
  }
}


void FermionMatrixOperations::performFFT(Complex* input, Complex* output, bool forw) {
  if ((Parameter_AntiPeriodicBoundaryConditionInTime) && (forw)) {
    multiplyFermionVectorWithTimeIndexedComplexScalars(input, output, halfMomentumForwardFFTfactors);

    fftw_plan plan = getFFTPlan(output, output, forw);
    fftw_execute(plan);   
    plan = NULL;
  }

  if ((Parameter_AntiPeriodicBoundaryConditionInTime) && (!forw)) {
    fftw_plan plan = getFFTPlan(input, output, forw);
    fftw_execute(plan);   
    plan = NULL;

    multiplyFermionVectorWithTimeIndexedComplexScalars(output, output, halfMomentumBackwardFFTfactors);
  }

  if (!Parameter_AntiPeriodicBoundaryConditionInTime) {
    if (xFFTusage) {
      if (xFFT==NULL) {
        xFFT = new ExtremeFFT4D(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, 8, 1, NULL, true);  
      }
      xFFT->FastFourierTrafo(input, output, forw);
    } else {
      fftw_plan plan = getFFTPlan(input, output, forw);
      fftw_execute(plan);   
      plan = NULL;
    }
  }
}


void FermionMatrixOperations::performDistributedFFT(DistributedMemoryObject* input, DistributedMemoryObject* output, bool forw) {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::performDistributedFFT: threadedOps==NULL\n");
    exit(0);
  }

  if ((Parameter_AntiPeriodicBoundaryConditionInTime) && (forw)) {
    threadedOps->multiplyFermionVectorWithTimeIndexedComplexScalars(input, output, Distributed_halfMomentumForwardFFTfactors, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);

    threadedOps->perform_FFTWFourierTransformationOfFermionVector(output, output, forw, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
  }

  if ((Parameter_AntiPeriodicBoundaryConditionInTime) && (!forw)) {
    threadedOps->perform_FFTWFourierTransformationOfFermionVector(input, output, forw, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);

    threadedOps->multiplyFermionVectorWithTimeIndexedComplexScalars(output, output, Distributed_halfMomentumBackwardFFTfactors, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
  }

  if (!Parameter_AntiPeriodicBoundaryConditionInTime) {
    if (xFFTusage) {
      threadedOps->perform_xFFTFourierTransformationOfFermionVector(input, output, forw, xFFT_DistributedFFT_ThreadCount, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3); 
    } else {
      threadedOps->perform_FFTWFourierTransformationOfFermionVector(input, output, forw, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
    }
  }  
}


void FermionMatrixOperations::tuneDistributedFFT(int tuneLevel) {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::tuneDistributedFFT: threadedOps==NULL\n");
    exit(0);
  }
  if (xFFTusage) {
    threadedOps->tune_xFFTFourierTransformationOfFermionVector(DistributedInverse_s, DistributedInverse_p, xFFT_DistributedFFT_ThreadCount, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, tuneLevel);
  } else {
  } 
}


void FermionMatrixOperations::tuneDistributedFFT(DistributedMemoryObject* input, DistributedMemoryObject* output, int tuneLevel) {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::tuneDistributedFFT: threadedOps==NULL\n");
    exit(0);
  }
  if (xFFTusage) {
    threadedOps->tune_xFFTFourierTransformationOfFermionVector(input, output, xFFT_DistributedFFT_ThreadCount, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, tuneLevel);
  } else {
  } 
}


void FermionMatrixOperations::tuneDistributedFFT(DistributedMemoryObject* input, DistributedMemoryObject* output, int tuneLevel, char* &fftPlanDescriptor) {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::tuneDistributedFFT: threadedOps==NULL\n");
    exit(0);
  }
  if (xFFTusage) {
    threadedOps->tune_xFFTFourierTransformationOfFermionVector(input, output, xFFT_DistributedFFT_ThreadCount, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, tuneLevel, fftPlanDescriptor);
  } else {
  }  
}


bool FermionMatrixOperations::setDistributedFFTPlan(char* fftPlanDescriptor) {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::setDistributedFFTPlan: threadedOps==NULL\n");
    exit(0);
  }
  if (xFFTusage) {
    return threadedOps->setFFTPlan_xFFTFourierTransformationOfFermionVector(xFFT_DistributedFFT_ThreadCount, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, fftPlanDescriptor);
  } else {
    return true;
  }  
}


fftw_plan FermionMatrixOperations::getFFTPlan(Complex* input, Complex* output, bool forw) {
  int I;
  for (I=0; I<FFTPlanDBMax; I++) {
    if ((InputVectorDataBase[I]==input) && (OutputVectorDataBase[I]==output) && (BackwardForwardDataBase[I]==forw)) {
      return FFTPlanDataBase[I];
    }
  }
  
  if (LogLevel>2) printf("Generating FFT-plan. Forward = %d\n", forw);  
  
  int* n = new int[4];
  n[0] = oneDimSizeL0;
  n[1] = oneDimSizeL1;
  n[2] = oneDimSizeL2;
  n[3] = oneDimSizeL3;
  
  int rank = 4;
  int howmany = 8;
  int* inembed = new int[4];
  inembed[0] = oneDimSizeL0;
  inembed[1] = oneDimSizeL1+xtraSize1;
  inembed[2] = oneDimSizeL2+xtraSize2;
  inembed[3] = oneDimSizeL3+xtraSize3;  
  int istride = 8;
  int idist = 1;
  
  int* onembed = new int[4];
  onembed[0] = oneDimSizeL0;
  onembed[1] = oneDimSizeL1+xtraSize1;
  onembed[2] = oneDimSizeL2+xtraSize2;
  onembed[3] = oneDimSizeL3+xtraSize3;    
  int ostride = 8;
  int odist = 1;

  if (FFTPlanDataBase[FFTDataBaseCounter] != NULL) {
    fftw_destroy_plan(FFTPlanDataBase[FFTDataBaseCounter]);
    FFTPlanDataBase[FFTDataBaseCounter] = NULL;
  }
  InputVectorDataBase[FFTDataBaseCounter] = input;
  OutputVectorDataBase[FFTDataBaseCounter] = output;
  BackwardForwardDataBase[FFTDataBaseCounter] = forw;
  
  Complex* Save1 = createFermionVector();
  Complex* Save2 = createFermionVector();
//  cblas_zcopy(vectorLengthXtrSize, input, 1, Save1, 1);
  SSE_ZCopy(vectorLengthXtrSize, input, 1, Save1, 1);
//  cblas_zcopy(vectorLengthXtrSize, output, 1, Save2, 1);
  SSE_ZCopy(vectorLengthXtrSize, output, 1, Save2, 1);

  int index = FFTDataBaseCounter;
  int measureFlag = FFTW_MEASURE | FFTW_EXHAUSTIVE;
  if (DebugMode) measureFlag = FFTW_MEASURE;
  if (forw == ExtremeFFT4D_Forward) {
    FFTPlanDataBase[index] = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*)input, inembed,
                                                istride, idist,
	   			                (fftw_complex*)output, onembed, ostride, odist,
	  				        FFTW_FORWARD, measureFlag);
  } else {
    FFTPlanDataBase[index] = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*)input, inembed,
                                                istride, idist,
	   				        (fftw_complex*)output, onembed, ostride, odist,
	  				        FFTW_BACKWARD, measureFlag);
  }
  FFTDataBaseCounter++;
  if (FFTDataBaseCounter>=FFTPlanDBMax) FFTDataBaseCounter = 0;
    
//  cblas_zcopy(vectorLengthXtrSize, Save1, 1, input, 1);
  SSE_ZCopy(vectorLengthXtrSize, Save1, 1, input, 1);
//  cblas_zcopy(vectorLengthXtrSize, Save2, 1, output, 1);
  SSE_ZCopy(vectorLengthXtrSize, Save2, 1, output, 1);
  destroyFermionVector(Save1);
  destroyFermionVector(Save2);
  delete[] n;
  delete[] inembed;
  delete[] onembed;
  return FFTPlanDataBase[index];
}


void FermionMatrixOperations::constructNeubergerWithXiFermionMatrix(ComplexMatrix& FermionMatrix, vector4D* phiField, double massSplitF) {
  if (Parameter_AntiPeriodicBoundaryConditionInTime) {
    printf("Error in FermionMatrixOperations::constructNeubergerWithXiFermionMatrix: Anti-Periodic Boundary Conditions not implemented yet!\n");
    exit(0);
  }

  if ((FermionMatrix.matrixSize!=vectorLength)) {
    FermionMatrix = ComplexMatrix(vectorLength);  
  }
  if ((NeubergerDiracOp->matrixSize!=vectorLength) || (NeubergerDiracOp->getRho()!=Parameter_rho) || (NeubergerDiracOp->getR()!=Parameter_r)) {
    delete NeubergerDiracOp;
    NeubergerDiracOp = new NeubergerMatrix(Parameter_rho, Parameter_r, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, 2);
  }

  ComplexMatrix Btilde(8);
  ComplexMatrix Dtilde(8);
  ComplexMatrix dummy(4);
  ComplexMatrix dummy8(8);
  int countZ = 0;
  int countS = 0;
  int z0,z1,z2,z3;
  int s0,s1,s2,s3;
  double f = massSplitF;
  
  if (LogLevel>2) printf("Constructing (Neuberger with Xi-fields) fermion matrix...");
  for (z0=0; z0<oneDimSizeL0; z0++) {
    for (z1=0; z1<oneDimSizeL1; z1++) {
      for (z2=0; z2<oneDimSizeL2; z2++) {
        for (z3=0; z3<oneDimSizeL3; z3++) {
          dummy.matrix[0][0].x = dummy.matrix[1][1].x = dummy.matrix[2][2].x = dummy.matrix[3][3].x = phiField[countZ][0];
          dummy.matrix[0][0].y = dummy.matrix[1][1].y = phiField[countZ][3];
	  dummy.matrix[2][2].y = dummy.matrix[3][3].y =-phiField[countZ][3];
	  Btilde.insertMatrix(dummy,0,0);
	  
          dummy.matrix[0][0].x = dummy.matrix[1][1].x = dummy.matrix[2][2].x = dummy.matrix[3][3].x = f*phiField[countZ][0];
          dummy.matrix[0][0].y = dummy.matrix[1][1].y =-f*phiField[countZ][3];
	  dummy.matrix[2][2].y = dummy.matrix[3][3].y = f*phiField[countZ][3];
	  Btilde.insertMatrix(dummy,4,4);
	  	  
          dummy.matrix[0][0].x = dummy.matrix[1][1].x = dummy.matrix[2][2].x = dummy.matrix[3][3].x = 0.5*(1-f)*phiField[countZ][2];
          dummy.matrix[0][0].y = dummy.matrix[1][1].y = dummy.matrix[2][2].y = dummy.matrix[3][3].y = 0.5*(1-f)*phiField[countZ][1];
          dummy.matrix[0][0].x += 0.5*(1+f)*phiField[countZ][2];
          dummy.matrix[0][0].y += 0.5*(1+f)*phiField[countZ][1];
          dummy.matrix[1][1].x += 0.5*(1+f)*phiField[countZ][2];
          dummy.matrix[1][1].y += 0.5*(1+f)*phiField[countZ][1];
          dummy.matrix[2][2].x -= 0.5*(1+f)*phiField[countZ][2];
          dummy.matrix[2][2].y -= 0.5*(1+f)*phiField[countZ][1];
          dummy.matrix[3][3].x -= 0.5*(1+f)*phiField[countZ][2];
          dummy.matrix[3][3].y -= 0.5*(1+f)*phiField[countZ][1];
	  Btilde.insertMatrix(dummy,0,4);

          dummy.matrix[0][0].x = dummy.matrix[1][1].x = dummy.matrix[2][2].x = dummy.matrix[3][3].x = 0.5*(1-f)*phiField[countZ][2];
          dummy.matrix[0][0].y = dummy.matrix[1][1].y = dummy.matrix[2][2].y = dummy.matrix[3][3].y = -0.5*(1-f)*phiField[countZ][1];
          dummy.matrix[0][0].x += -0.5*(1+f)*phiField[countZ][2];
          dummy.matrix[0][0].y += 0.5*(1+f)*phiField[countZ][1];
          dummy.matrix[1][1].x += -0.5*(1+f)*phiField[countZ][2];
          dummy.matrix[1][1].y += 0.5*(1+f)*phiField[countZ][1];
          dummy.matrix[2][2].x -= -0.5*(1+f)*phiField[countZ][2];
          dummy.matrix[2][2].y -= 0.5*(1+f)*phiField[countZ][1];
          dummy.matrix[3][3].x -= -0.5*(1+f)*phiField[countZ][2];
          dummy.matrix[3][3].y -= 0.5*(1+f)*phiField[countZ][1];
	  Btilde.insertMatrix(dummy,4,0);

	  
	  countS = 0;
          for (s0=0; s0<oneDimSizeL0; s0++) {
            for (s1=0; s1<oneDimSizeL1; s1++) {
              for (s2=0; s2<oneDimSizeL2; s2++) {
                for (s3=0; s3<oneDimSizeL3; s3++) {
                  Dtilde = NeubergerDiracOp->getSubMatrix(countZ*8,countS*8,8);
		  dummy8 = Parameter_yN*Btilde*Dtilde - (2*Parameter_rho)*Dtilde;
		  if (countZ == countS) {
 		    dummy8 = dummy8 - (2*Parameter_rho*Parameter_yN)*Btilde;
                    dummy8.matrix[0][0].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    dummy8.matrix[1][1].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    dummy8.matrix[2][2].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    dummy8.matrix[3][3].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    dummy8.matrix[4][4].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    dummy8.matrix[5][5].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    dummy8.matrix[6][6].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    dummy8.matrix[7][7].x -= 2*Parameter_rho*Parameter_ExplicitMass;
		  }
		  
		  FermionMatrix.insertMatrix(dummy8,countZ*8,countS*8);

                  countS++;		
		}
              }
	    }
	  }

          countZ++;
	}
      }
    }
  }

  if (LogLevel>2) printf("ready.\n");
}


void FermionMatrixOperations::constructNeubergerWithXiFermionMatrix(ComplexMatrix& FermionMatrix, vector4D* phiField) {
  if (Parameter_AntiPeriodicBoundaryConditionInTime) {
    printf("Error in FermionMatrixOperations::constructNeubergerWithXiFermionMatrix: Anti-Periodic Boundary Conditions not implemented yet!\n");
    exit(0);
  }

  if ((FermionMatrix.matrixSize!=vectorLength)) {
    FermionMatrix = ComplexMatrix(vectorLength);  
  }
  if ((NeubergerDiracOp->matrixSize!=vectorLength) || (NeubergerDiracOp->getRho()!=Parameter_rho) || (NeubergerDiracOp->getR()!=Parameter_r)) {
    delete NeubergerDiracOp;
    NeubergerDiracOp = new NeubergerMatrix(Parameter_rho, Parameter_r, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, 2);
  }

  ComplexMatrix Btilde(8);
  ComplexMatrix Dtilde(8);
  ComplexMatrix dummy(4);
  ComplexMatrix dummy8(8);
  int countZ = 0;
  int countS = 0;
  int z0,z1,z2,z3;
  int s0,s1,s2,s3;
  
  if (LogLevel>2) printf("Constructing (Neuberger with Xi-fields) fermion matrix...");
  for (z0=0; z0<oneDimSizeL0; z0++) {
    for (z1=0; z1<oneDimSizeL1; z1++) {
      for (z2=0; z2<oneDimSizeL2; z2++) {
        for (z3=0; z3<oneDimSizeL3; z3++) {
          dummy.matrix[0][0].x = dummy.matrix[1][1].x = dummy.matrix[2][2].x = dummy.matrix[3][3].x = phiField[countZ][0];
          dummy.matrix[0][0].y = dummy.matrix[1][1].y = phiField[countZ][3];
	  dummy.matrix[2][2].y = dummy.matrix[3][3].y =-phiField[countZ][3];
	  Btilde.insertMatrix(dummy,0,0);
	  
          dummy.matrix[0][0].y = dummy.matrix[1][1].y =-phiField[countZ][3];
	  dummy.matrix[2][2].y = dummy.matrix[3][3].y = phiField[countZ][3];
	  Btilde.insertMatrix(dummy,4,4);
	  
          dummy.matrix[0][0].x = dummy.matrix[1][1].x = phiField[countZ][2];
	  dummy.matrix[2][2].x = dummy.matrix[3][3].x =-phiField[countZ][2];
          dummy.matrix[0][0].y = dummy.matrix[1][1].y = phiField[countZ][1];
	  dummy.matrix[2][2].y = dummy.matrix[3][3].y =-phiField[countZ][1];
	  Btilde.insertMatrix(dummy,0,4);

          dummy.matrix[0][0].x = dummy.matrix[1][1].x =-phiField[countZ][2];
	  dummy.matrix[2][2].x = dummy.matrix[3][3].x = phiField[countZ][2];
	  Btilde.insertMatrix(dummy,4,0);
	  
	  countS = 0;
          for (s0=0; s0<oneDimSizeL0; s0++) {
            for (s1=0; s1<oneDimSizeL1; s1++) {
              for (s2=0; s2<oneDimSizeL2; s2++) {
                for (s3=0; s3<oneDimSizeL3; s3++) {
                  Dtilde = NeubergerDiracOp->getSubMatrix(countZ*8,countS*8,8);
		  dummy8 = Parameter_yN*Btilde*Dtilde - (2*Parameter_rho)*Dtilde;
		  if (countZ == countS) {
 		    dummy8 = dummy8 - (2*Parameter_rho*Parameter_yN)*Btilde;
                    dummy8.matrix[0][0].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    dummy8.matrix[1][1].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    dummy8.matrix[2][2].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    dummy8.matrix[3][3].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    dummy8.matrix[4][4].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    dummy8.matrix[5][5].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    dummy8.matrix[6][6].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    dummy8.matrix[7][7].x -= 2*Parameter_rho*Parameter_ExplicitMass;
		  }
		  
		  FermionMatrix.insertMatrix(dummy8,countZ*8,countS*8);

                  countS++;		
		}
              }
	    }
	  }

          countZ++;
	}
      }
    }
  }

  if (LogLevel>2) printf("ready.\n");
}


void FermionMatrixOperations::constructNeubergerWithOutXiFermionMatrix(ComplexMatrix& FermionMatrix, vector4D* phiField) {
  if (Parameter_AntiPeriodicBoundaryConditionInTime) {
    printf("Error in FermionMatrixOperations::constructNeubergerWithOutXiFermionMatrix: Anti-Periodic Boundary Conditions not implemented yet!\n");
    exit(0);
  }

  if ((FermionMatrix.matrixSize!=vectorLength)) {
    FermionMatrix = ComplexMatrix(vectorLength);  
  }
  if ((NeubergerDiracOp->matrixSize!=vectorLength) || (NeubergerDiracOp->getRho()!=Parameter_rho) || (NeubergerDiracOp->getR()!=Parameter_r)) {
    delete NeubergerDiracOp;
    NeubergerDiracOp = new NeubergerMatrix(Parameter_rho, Parameter_r, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, 2);
  }

  ComplexMatrix Btilde(8);
  ComplexMatrix Dtilde(8);
  ComplexMatrix dummy(4);
  int countZ = 0;
  int countS = 0;
  int z0,z1,z2,z3;
  int s0,s1,s2,s3;
  
  if (LogLevel>2) printf("Constructing (Neuberger without Xi-fields) fermion matrix...");
  
  for (z0=0; z0<oneDimSizeL0; z0++) {
    for (z1=0; z1<oneDimSizeL1; z1++) {
      for (z2=0; z2<oneDimSizeL2; z2++) {
        for (z3=0; z3<oneDimSizeL3; z3++) {
          dummy.matrix[0][0].x = dummy.matrix[1][1].x = dummy.matrix[2][2].x = dummy.matrix[3][3].x = phiField[countZ][0];
          dummy.matrix[0][0].y = dummy.matrix[1][1].y = phiField[countZ][3];
	  dummy.matrix[2][2].y = dummy.matrix[3][3].y =-phiField[countZ][3];
	  Btilde.insertMatrix(dummy,0,0);
	  
          dummy.matrix[0][0].y = dummy.matrix[1][1].y =-phiField[countZ][3];
	  dummy.matrix[2][2].y = dummy.matrix[3][3].y = phiField[countZ][3];
	  Btilde.insertMatrix(dummy,4,4);
	  
          dummy.matrix[0][0].x = dummy.matrix[1][1].x = phiField[countZ][2];
	  dummy.matrix[2][2].x = dummy.matrix[3][3].x =-phiField[countZ][2];
          dummy.matrix[0][0].y = dummy.matrix[1][1].y = phiField[countZ][1];
	  dummy.matrix[2][2].y = dummy.matrix[3][3].y =-phiField[countZ][1];
	  Btilde.insertMatrix(dummy,0,4);

          dummy.matrix[0][0].x = dummy.matrix[1][1].x =-phiField[countZ][2];
	  dummy.matrix[2][2].x = dummy.matrix[3][3].x = phiField[countZ][2];
	  Btilde.insertMatrix(dummy,4,0);
	  
	  countS = 0;
          for (s0=0; s0<oneDimSizeL0; s0++) {
            for (s1=0; s1<oneDimSizeL1; s1++) {
              for (s2=0; s2<oneDimSizeL2; s2++) {
                for (s3=0; s3<oneDimSizeL3; s3++) {
                  Dtilde = NeubergerDiracOp->getSubMatrix(countZ*8,countS*8,8);
		  if (countZ == countS) {
 		    Dtilde = Dtilde + (Parameter_yN*Btilde);
                    Dtilde.matrix[0][0].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    Dtilde.matrix[1][1].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    Dtilde.matrix[2][2].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    Dtilde.matrix[3][3].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    Dtilde.matrix[4][4].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    Dtilde.matrix[5][5].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    Dtilde.matrix[6][6].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    Dtilde.matrix[7][7].x -= 2*Parameter_rho*Parameter_ExplicitMass;
		  }
		  
		  FermionMatrix.insertMatrix(Dtilde,countZ*8,countS*8);

                  countS++;		
		}
              }
	    }
	  }

          countZ++;
	}
      }
    }
  }
  if (LogLevel>2) printf("ready.\n");
}



void FermionMatrixOperations::constructNeubergerWithOutXiFermionMatrix(ComplexMatrix& FermionMatrix, vector4D* phiField, double massSplitF) {
  if (Parameter_AntiPeriodicBoundaryConditionInTime) {
    printf("Error in FermionMatrixOperations::constructNeubergerWithOutXiFermionMatrix: Anti-Periodic Boundary Conditions not implemented yet!\n");
    exit(0);
  }

  if ((FermionMatrix.matrixSize!=vectorLength)) {
    FermionMatrix = ComplexMatrix(vectorLength);  
  }
  if ((NeubergerDiracOp->matrixSize!=vectorLength) || (NeubergerDiracOp->getRho()!=Parameter_rho) || (NeubergerDiracOp->getR()!=Parameter_r)) {
    delete NeubergerDiracOp;
    NeubergerDiracOp = new NeubergerMatrix(Parameter_rho, Parameter_r, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, 2);
  }

  ComplexMatrix Btilde(8);
  ComplexMatrix Dtilde(8);
  ComplexMatrix dummy(4);
  int countZ = 0;
  int countS = 0;
  int z0,z1,z2,z3;
  int s0,s1,s2,s3;
  double f = massSplitF;
  
  if (LogLevel>2) printf("Constructing (Neuberger without Xi-fields) fermion matrix...");
  
  for (z0=0; z0<oneDimSizeL0; z0++) {
    for (z1=0; z1<oneDimSizeL1; z1++) {
      for (z2=0; z2<oneDimSizeL2; z2++) {
        for (z3=0; z3<oneDimSizeL3; z3++) {
          dummy.matrix[0][0].x = dummy.matrix[1][1].x = dummy.matrix[2][2].x = dummy.matrix[3][3].x = phiField[countZ][0];
          dummy.matrix[0][0].y = dummy.matrix[1][1].y = phiField[countZ][3];
	  dummy.matrix[2][2].y = dummy.matrix[3][3].y =-phiField[countZ][3];
	  Btilde.insertMatrix(dummy,0,0);
	  
          dummy.matrix[0][0].x = dummy.matrix[1][1].x = dummy.matrix[2][2].x = dummy.matrix[3][3].x = f*phiField[countZ][0];
          dummy.matrix[0][0].y = dummy.matrix[1][1].y =-f*phiField[countZ][3];
	  dummy.matrix[2][2].y = dummy.matrix[3][3].y = f*phiField[countZ][3];
	  Btilde.insertMatrix(dummy,4,4);
	  	  
          dummy.matrix[0][0].x = dummy.matrix[1][1].x = dummy.matrix[2][2].x = dummy.matrix[3][3].x = 0.5*(1-f)*phiField[countZ][2];
          dummy.matrix[0][0].y = dummy.matrix[1][1].y = dummy.matrix[2][2].y = dummy.matrix[3][3].y = 0.5*(1-f)*phiField[countZ][1];
          dummy.matrix[0][0].x += 0.5*(1+f)*phiField[countZ][2];
          dummy.matrix[0][0].y += 0.5*(1+f)*phiField[countZ][1];
          dummy.matrix[1][1].x += 0.5*(1+f)*phiField[countZ][2];
          dummy.matrix[1][1].y += 0.5*(1+f)*phiField[countZ][1];
          dummy.matrix[2][2].x -= 0.5*(1+f)*phiField[countZ][2];
          dummy.matrix[2][2].y -= 0.5*(1+f)*phiField[countZ][1];
          dummy.matrix[3][3].x -= 0.5*(1+f)*phiField[countZ][2];
          dummy.matrix[3][3].y -= 0.5*(1+f)*phiField[countZ][1];
	  Btilde.insertMatrix(dummy,0,4);

          dummy.matrix[0][0].x = dummy.matrix[1][1].x = dummy.matrix[2][2].x = dummy.matrix[3][3].x = 0.5*(1-f)*phiField[countZ][2];
          dummy.matrix[0][0].y = dummy.matrix[1][1].y = dummy.matrix[2][2].y = dummy.matrix[3][3].y = -0.5*(1-f)*phiField[countZ][1];
          dummy.matrix[0][0].x += -0.5*(1+f)*phiField[countZ][2];
          dummy.matrix[0][0].y += 0.5*(1+f)*phiField[countZ][1];
          dummy.matrix[1][1].x += -0.5*(1+f)*phiField[countZ][2];
          dummy.matrix[1][1].y += 0.5*(1+f)*phiField[countZ][1];
          dummy.matrix[2][2].x -= -0.5*(1+f)*phiField[countZ][2];
          dummy.matrix[2][2].y -= 0.5*(1+f)*phiField[countZ][1];
          dummy.matrix[3][3].x -= -0.5*(1+f)*phiField[countZ][2];
          dummy.matrix[3][3].y -= 0.5*(1+f)*phiField[countZ][1];
	  Btilde.insertMatrix(dummy,4,0);	
		  
	  countS = 0;
          for (s0=0; s0<oneDimSizeL0; s0++) {
            for (s1=0; s1<oneDimSizeL1; s1++) {
              for (s2=0; s2<oneDimSizeL2; s2++) {
                for (s3=0; s3<oneDimSizeL3; s3++) {
                  Dtilde = NeubergerDiracOp->getSubMatrix(countZ*8,countS*8,8);
		  if (countZ == countS) {
 		    Dtilde = Dtilde + (Parameter_yN*Btilde);
                    Dtilde.matrix[0][0].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    Dtilde.matrix[1][1].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    Dtilde.matrix[2][2].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    Dtilde.matrix[3][3].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    Dtilde.matrix[4][4].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    Dtilde.matrix[5][5].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    Dtilde.matrix[6][6].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    Dtilde.matrix[7][7].x -= 2*Parameter_rho*Parameter_ExplicitMass;
		  }
		  
		  FermionMatrix.insertMatrix(Dtilde,countZ*8,countS*8);

                  countS++;		
		}
              }
	    }
	  }

          countZ++;
	}
      }
    }
  }
  if (LogLevel>2) printf("ready.\n");
}


void FermionMatrixOperations::constructWilsonWithOutXiFermionMatrix(ComplexMatrix& FermionMatrix, vector4D* phiField) {
  if (Parameter_AntiPeriodicBoundaryConditionInTime) {
    printf("Error in FermionMatrixOperations::constructWilsonWithOutXiFermionMatrix: Anti-Periodic Boundary Conditions not implemented yet!\n");
    exit(0);
  }

  if ((FermionMatrix.matrixSize!=vectorLength)) {
    FermionMatrix = ComplexMatrix(vectorLength);  
  }
  if ((WilsonDiracOp->matrixSize!=vectorLength) || (WilsonDiracOp->getR()!=Parameter_r)) {
    delete WilsonDiracOp;
    WilsonDiracOp = new WilsonMatrix(Parameter_r, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, 2);
  }

  ComplexMatrix Btilde(8);
  ComplexMatrix Dtilde(8);
  ComplexMatrix dummy(4);
  int countZ = 0;
  int countS = 0;
  int z0,z1,z2,z3;
  int s0,s1,s2,s3;
  
  if (LogLevel>2) printf("Constructing (Wilson without Xi-fields) fermion matrix...");
  
  for (z0=0; z0<oneDimSizeL0; z0++) {
    for (z1=0; z1<oneDimSizeL1; z1++) {
      for (z2=0; z2<oneDimSizeL2; z2++) {
        for (z3=0; z3<oneDimSizeL3; z3++) {
          dummy.matrix[0][0].x = dummy.matrix[1][1].x = dummy.matrix[2][2].x = dummy.matrix[3][3].x = phiField[countZ][0];
          dummy.matrix[0][0].y = dummy.matrix[1][1].y = phiField[countZ][3];
	  dummy.matrix[2][2].y = dummy.matrix[3][3].y =-phiField[countZ][3];
	  Btilde.insertMatrix(dummy,0,0);
	  
          dummy.matrix[0][0].y = dummy.matrix[1][1].y =-phiField[countZ][3];
	  dummy.matrix[2][2].y = dummy.matrix[3][3].y = phiField[countZ][3];
	  Btilde.insertMatrix(dummy,4,4);
	  
          dummy.matrix[0][0].x = dummy.matrix[1][1].x = phiField[countZ][2];
	  dummy.matrix[2][2].x = dummy.matrix[3][3].x =-phiField[countZ][2];
          dummy.matrix[0][0].y = dummy.matrix[1][1].y = phiField[countZ][1];
	  dummy.matrix[2][2].y = dummy.matrix[3][3].y =-phiField[countZ][1];
	  Btilde.insertMatrix(dummy,0,4);

          dummy.matrix[0][0].x = dummy.matrix[1][1].x =-phiField[countZ][2];
	  dummy.matrix[2][2].x = dummy.matrix[3][3].x = phiField[countZ][2];
	  Btilde.insertMatrix(dummy,4,0);
	  
	  countS = 0;
          for (s0=0; s0<oneDimSizeL0; s0++) {
            for (s1=0; s1<oneDimSizeL1; s1++) {
              for (s2=0; s2<oneDimSizeL2; s2++) {
                for (s3=0; s3<oneDimSizeL3; s3++) {
                  Dtilde = WilsonDiracOp->getSubMatrix(countZ*8,countS*8,8);
		  if (countZ == countS) {
 		    Dtilde = Dtilde + (Parameter_yN*Btilde);
                    Dtilde.matrix[0][0].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    Dtilde.matrix[1][1].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    Dtilde.matrix[2][2].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    Dtilde.matrix[3][3].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    Dtilde.matrix[4][4].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    Dtilde.matrix[5][5].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    Dtilde.matrix[6][6].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    Dtilde.matrix[7][7].x -= 2*Parameter_rho*Parameter_ExplicitMass;
		  }
		  
		  FermionMatrix.insertMatrix(Dtilde,countZ*8,countS*8);

                  countS++;		
		}
              }
	    }
	  }

          countZ++;
	}
      }
    }
  }
  if (LogLevel>2) printf("ready.\n");
}


void FermionMatrixOperations::constructWilsonWithXiFermionMatrix(ComplexMatrix& FermionMatrix, vector4D* phiField) {
  if ((FermionMatrix.matrixSize!=vectorLength)) {
    FermionMatrix = ComplexMatrix(vectorLength);  
  }
  if ((WilsonDiracOp->matrixSize!=vectorLength) || (WilsonDiracOp->getR()!=Parameter_r)) {
    delete WilsonDiracOp;
    WilsonDiracOp = new WilsonMatrix(Parameter_r, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, 2);
  }

  ComplexMatrix Btilde(8);
  ComplexMatrix Dtilde(8);
  ComplexMatrix dummy(4);
  ComplexMatrix dummy8(8);
  int countZ = 0;
  int countS = 0;
  int z0,z1,z2,z3;
  int s0,s1,s2,s3;
  
  if (LogLevel>2) printf("Constructing (Wilson with Xi-fields) fermion matrix...");
  
  for (z0=0; z0<oneDimSizeL0; z0++) {
    for (z1=0; z1<oneDimSizeL1; z1++) {
      for (z2=0; z2<oneDimSizeL2; z2++) {
        for (z3=0; z3<oneDimSizeL3; z3++) {
          dummy.matrix[0][0].x = dummy.matrix[1][1].x = dummy.matrix[2][2].x = dummy.matrix[3][3].x = phiField[countZ][0];
          dummy.matrix[0][0].y = dummy.matrix[1][1].y = phiField[countZ][3];
	  dummy.matrix[2][2].y = dummy.matrix[3][3].y =-phiField[countZ][3];
	  Btilde.insertMatrix(dummy,0,0);
	  
          dummy.matrix[0][0].y = dummy.matrix[1][1].y =-phiField[countZ][3];
	  dummy.matrix[2][2].y = dummy.matrix[3][3].y = phiField[countZ][3];
	  Btilde.insertMatrix(dummy,4,4);
	  
          dummy.matrix[0][0].x = dummy.matrix[1][1].x = phiField[countZ][2];
	  dummy.matrix[2][2].x = dummy.matrix[3][3].x =-phiField[countZ][2];
          dummy.matrix[0][0].y = dummy.matrix[1][1].y = phiField[countZ][1];
	  dummy.matrix[2][2].y = dummy.matrix[3][3].y =-phiField[countZ][1];
	  Btilde.insertMatrix(dummy,0,4);

          dummy.matrix[0][0].x = dummy.matrix[1][1].x =-phiField[countZ][2];
	  dummy.matrix[2][2].x = dummy.matrix[3][3].x = phiField[countZ][2];
	  Btilde.insertMatrix(dummy,4,0);
	  
	  countS = 0;
          for (s0=0; s0<oneDimSizeL0; s0++) {
            for (s1=0; s1<oneDimSizeL1; s1++) {
              for (s2=0; s2<oneDimSizeL2; s2++) {
                for (s3=0; s3<oneDimSizeL3; s3++) {
                  Dtilde = WilsonDiracOp->getSubMatrix(countZ*8,countS*8,8);
		  dummy8 = Parameter_yN*Btilde*Dtilde - (2*Parameter_rho)*Dtilde;
		  if (countZ == countS) {
 		    dummy8 = dummy8 - (2*Parameter_rho*Parameter_yN)*Btilde;
                    dummy8.matrix[0][0].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    dummy8.matrix[1][1].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    dummy8.matrix[2][2].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    dummy8.matrix[3][3].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    dummy8.matrix[4][4].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    dummy8.matrix[5][5].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    dummy8.matrix[6][6].x -= 2*Parameter_rho*Parameter_ExplicitMass;
                    dummy8.matrix[7][7].x -= 2*Parameter_rho*Parameter_ExplicitMass;
		  }
		  
		  FermionMatrix.insertMatrix(dummy8,countZ*8,countS*8);

                  countS++;		
		}
              }
	    }
	  }

          countZ++;
	}
      }
    }
  }
  if (LogLevel>2) printf("ready.\n");
}


void FermionMatrixOperations::executeDiracUnityMinusDiracQuasiInverseMatrixMultiplication(Complex* input, Complex* output, bool inFourierSpace) {
  Complex* source;
  Complex* target;
  Complex* EWData;
  if (inFourierSpace) {
    source = input;
    target = output;  
    EWData = NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWDataRescaled;
  } else {
    performFFT(input, InputVectorFourierTransform, ExtremeFFT4D_Forward); 
    source = InputVectorFourierTransform;
    target = InputVectorFourierTransform;  
    EWData = NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWData;
  }
 
  SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, source, target, 
   NeubergerDiracOpApplicationSinPStdData, EWData);

  if (!inFourierSpace) {
    performFFT(InputVectorFourierTransform, output, ExtremeFFT4D_Backward);
  }
}


void FermionMatrixOperations::executeDistributedDiracUnityMinusDiracQuasiInverseMatrixMultiplication(DistributedMemoryObject* input, DistributedMemoryObject* output, bool inFourierSpace) {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::executeDistributedDiracUnityMinusDiracQuasiInverseMatrixMultiplication: threadedOps==NULL\n");
    exit(0);
  }
  DistributedMemoryObject* source;
  DistributedMemoryObject* target;
  DistributedMemoryObject* EWData;
  if (inFourierSpace) {
    source = input;
    target = output;  
    EWData = Distributed_NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWDataRescaled;
  } else {
    performDistributedFFT(input, DistributedInputVectorFourierTransform, ExtremeFFT4D_Forward); 
    source = DistributedInputVectorFourierTransform;
    target = DistributedInputVectorFourierTransform;  
    EWData = Distributed_NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWData;
  }
 
  threadedOps->perform_DiracTypeOperatorMultiplicationInFourierSpace(source, target, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, Distributed_NeubergerDiracOpApplicationSinPStdData, EWData);

  if (!inFourierSpace) {
    performDistributedFFT(DistributedInputVectorFourierTransform, output, ExtremeFFT4D_Backward);
  }
}


void FermionMatrixOperations::executeDiracUnityMinusDiracQuasiInverseDaggerMatrixMultiplication(Complex* input, Complex* output, bool inFourierSpace) {
  Complex* source;
  Complex* target;
  Complex* EWData;
  if (inFourierSpace) {
    source = input;
    target = output;  
    EWData = NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWDataRescaled;
  } else {
    performFFT(input, InputVectorFourierTransform, ExtremeFFT4D_Forward); 
    source = InputVectorFourierTransform;
    target = InputVectorFourierTransform;  
    EWData = NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWData;
  }
 
  SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, source, target, 
   NeubergerDiracOpApplicationSinPDagData, EWData);

  if (!inFourierSpace) {
    performFFT(InputVectorFourierTransform, output, ExtremeFFT4D_Backward);
  }
}


void FermionMatrixOperations::executeDistributedDiracUnityMinusDiracQuasiInverseDaggerMatrixMultiplication(DistributedMemoryObject* input, DistributedMemoryObject* output, bool inFourierSpace) {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::executeDistributedDiracUnityMinusDiracQuasiInverseDaggerMatrixMultiplication: threadedOps==NULL\n");
    exit(0);
  }
  DistributedMemoryObject* source;
  DistributedMemoryObject* target;
  DistributedMemoryObject* EWData;
  if (inFourierSpace) {
    source = input;
    target = output;  
    EWData = Distributed_NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWDataRescaled;
  } else {
    performDistributedFFT(input, DistributedInputVectorFourierTransform, ExtremeFFT4D_Forward); 
    source = DistributedInputVectorFourierTransform;
    target = DistributedInputVectorFourierTransform;  
    EWData = Distributed_NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWData;
  }
 
  threadedOps->perform_DiracTypeOperatorMultiplicationInFourierSpace(source, target, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, Distributed_NeubergerDiracOpApplicationSinPDagData, EWData);
 
  if (!inFourierSpace) {
    performDistributedFFT(DistributedInputVectorFourierTransform, output, ExtremeFFT4D_Backward);
  }
}


void FermionMatrixOperations::executeDistributedMultiCombinedOperator_VA1_RorCopy_VAc_RorCopy_D_MatrixMultiplication(DistributedMemoryObject* v1, DistributedMemoryObject* v2, DistributedMemoryObject* v3, DistributedMemoryObject* v4, Complex alpha2, bool Ddag, bool inFourierSpace) {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::executeDistributedMultiCombinedOperator_VA1_RorCopy_VAc_RorCopy_D_InFourierSpace: threadedOps==NULL\n");
    exit(0);
  }
  if (!inFourierSpace) {
    printf("ERROR in FermionMatrixOperations::executeDistributedMultiCombinedOperator_VA1_RorCopy_VAc_RorCopy_D_InFourierSpace: only implemented for Fourier-Space\n");
    exit(0);
  }

  DistributedMemoryObject* EWData;
  DistributedMemoryObject* sinP;
  
  if (Ddag) {
    EWData = Distributed_NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWDataRescaled;
    sinP = Distributed_NeubergerDiracOpApplicationSinPDagData;
  } else {
    EWData = Distributed_NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWDataRescaled;
    sinP = Distributed_NeubergerDiracOpApplicationSinPStdData;
  }
 
  threadedOps->perform_MultiCombinedOperator_VA1_RorCopy_VAc_RorCopy_D_InFourierSpace(v1, v2, v3, v4, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, RPreconditioner_Usage, alpha2, sinP, EWData, Distributed_RPreconditionerDiagonalData, Distributed_RPreconditionerDiagonalData);
}


//Berechne Projector \hat P_\pm angewendet auf input  ==> OHNE Faktor 0.5
void FermionMatrixOperations::executeProjectorHatMultiplication(Complex* input, Complex* output, bool projPlus, bool daggered) {  
  Complex alphaPsign(-1, 0);
  if (projPlus) alphaPsign.x = 1;
  Complex alpha(-1.0/Parameter_rho, 0);
  
  if (daggered) {
    executeGamma5(input, Inverse_p);
    executeDiracDaggerMatrixMultiplication(Inverse_p, Inverse_s, false);
    SSE_ComplexVectorAddition(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, alpha, Inverse_s, Inverse_p);

    if (input != output) {
      SSE_ZCopy(vectorLengthXtrSize, input, 1, output, 1);      
    }
    SSE_ComplexVectorAddition(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, alphaPsign, Inverse_p, output);	
  } else {
    executeDiracMatrixMultiplication(input, Inverse_p, false);
    SSE_ZCopy(vectorLengthXtrSize, input, 1, Inverse_s, 1);      
    SSE_ComplexVectorAddition(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, alpha, Inverse_p, Inverse_s);	
    executeGamma5(Inverse_s, Inverse_s);
    
    if (input != output) {
      SSE_ZCopy(vectorLengthXtrSize, input, 1, output, 1);      
    }
    SSE_ComplexVectorAddition(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, alphaPsign, Inverse_s, output);	
  }
}


//Berechne Projector P_\pm angewendet auf input  ==> OHNE Faktor 0.5
void FermionMatrixOperations::executeProjectorMultiplication(Complex* input, Complex* output, bool projPlus) {
  Complex alphaPsign(-1, 0);
  if (projPlus) alphaPsign.x = 1;

  executeGamma5(input, Inverse_p);

  if (input != output) {
    SSE_ZCopy(vectorLengthXtrSize, input, 1, output, 1);
  }
  SSE_ComplexVectorAddition(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, alphaPsign, Inverse_p, output);
}


void FermionMatrixOperations::executeDiracMatrixMultiplication(Complex* input, Complex* output, bool inFourierSpace) {
  Complex* source;
  Complex* target;
  Complex* EWData;
  if (inFourierSpace) {
    source = input;
    target = output;  
    EWData = NeubergerDiracOpApplicationEWDataRescaled;
  } else {
    performFFT(input, InputVectorFourierTransform, ExtremeFFT4D_Forward); 
    source = InputVectorFourierTransform;
    target = InputVectorFourierTransform;  
    EWData = NeubergerDiracOpApplicationEWData;
  }
 
  SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, source, target, 
   NeubergerDiracOpApplicationSinPStdData, EWData);

  if (!inFourierSpace) {
    performFFT(InputVectorFourierTransform, output, ExtremeFFT4D_Backward);
  }
}


void FermionMatrixOperations::executeDiracDaggerMatrixMultiplication(Complex* input, Complex* output, bool inFourierSpace) {
  Complex* source;
  Complex* target;
  Complex* EWData;
  if (inFourierSpace) {
    source = input;
    target = output;  
    EWData = NeubergerDiracOpApplicationEWDataRescaled;
  } else {
    performFFT(input, InputVectorFourierTransform, ExtremeFFT4D_Forward); 
    source = InputVectorFourierTransform;
    target = InputVectorFourierTransform;  
    EWData = NeubergerDiracOpApplicationEWData;
  }
 
  SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, source, target, 
   NeubergerDiracOpApplicationSinPDagData, EWData);

  if (!inFourierSpace) {
    performFFT(InputVectorFourierTransform, output, ExtremeFFT4D_Backward);
  }
}


void FermionMatrixOperations::executeDiracDaggerPlusExplicitMassMatrixMultiplication(Complex* input, Complex* output, bool inFourierSpace) {
  Complex* source;
  Complex* target;
  Complex* EWData;
  if (inFourierSpace) {
    source = input;
    target = output;  
    EWData = NeubergerDiracPlusExplicitMassOpApplicationEWDataRescaled;
  } else {
    performFFT(input, InputVectorFourierTransform, ExtremeFFT4D_Forward); 
    source = InputVectorFourierTransform;
    target = InputVectorFourierTransform;  
    EWData = NeubergerDiracPlusExplicitMassOpApplicationEWData;
  }
 
  SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, source, target, 
   NeubergerDiracOpApplicationSinPDagData, EWData);

  if (!inFourierSpace) {
    performFFT(InputVectorFourierTransform, output, ExtremeFFT4D_Backward);
  }
}


void FermionMatrixOperations::executeFermionQuasiHermiteanMatrixMultiplication(Complex* input, Complex* output, double* phiField, bool useRPrec, bool daggered, bool inFourierSpace) {
  Complex* p1 = NULL;
  Complex* p2 = NULL;
  Complex* p3 = NULL;
  Complex* p4 = NULL;  
  Complex* p5 = NULL;  
  double vol = oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3;
  double fac = 1.0;  
  useRPrec = useRPrec && RPreconditioner_Usage;
  
  if ((xFFTusage) && (inFourierSpace) && (useRPrec)) {
    p1 = input;
    p2 = InputVectorFourierTransform;
    p3 = InputVectorFourierTransform;
    p4 = output;
    p5 = InterimVectorFourierTransform;
  } 
  if ((xFFTusage) && (inFourierSpace) && (!useRPrec)) {
    p2 = input;
    p3 = InputVectorFourierTransform;
    p4 = output;
    p5 = InterimVectorFourierTransform;
  } 
  if ((xFFTusage) && (!inFourierSpace) && (useRPrec)) {
    p1 = InputVectorFourierTransform;
    p2 = InputVectorFourierTransform;
    p3 = output;
    p4 = InputVectorFourierTransform;
    p5 = InterimVectorFourierTransform;
    fac = 1.0 / vol;
  } 
  if ((xFFTusage) && (!inFourierSpace) && (!useRPrec)) {
    p2 = InputVectorFourierTransform;
    p3 = output;
    p4 = InputVectorFourierTransform;
    p5 = InterimVectorFourierTransform;
    fac = 1.0 / vol;
  } 
  if ((!xFFTusage) && (inFourierSpace) && (useRPrec)) {
    p1 = input;
    p2 = InputVectorFourierTransform;
    p3 = InputVectorFourierTransform;
    p4 = output;
    p5 = InputVectorFourierTransform;
  } 
  if ((!xFFTusage) && (inFourierSpace) && (!useRPrec)) {
    p2 = input;
    p3 = InputVectorFourierTransform;
    p4 = output;
    p5 = InputVectorFourierTransform;
  } 
  if ((!xFFTusage) && (!inFourierSpace) && (useRPrec)) {
    p1 = InputVectorFourierTransform;
    p2 = InputVectorFourierTransform;
    p3 = output;
    p4 = InputVectorFourierTransform;
    p5 = InputVectorFourierTransform;
    fac = 1.0 / vol;
  } 
  if ((!xFFTusage) && (!inFourierSpace) && (!useRPrec)) {
    p2 = InputVectorFourierTransform;
    p3 = output;
    p4 = InputVectorFourierTransform;
    p5 = InputVectorFourierTransform;
    fac = 1.0 / vol;
  } 

  if (!inFourierSpace) performFFT(input, InputVectorFourierTransform, ExtremeFFT4D_Forward);  
    
  if (useRPrec) executeRPreconditionerMatrixMultiplication(p1, InputVectorFourierTransform, false, false, true);  
    
  if (p2 == input) {
    executePiModeRemoverOperator(p2, InputVectorFourierTransform, true);  
    if (daggered) {
      executeDiracUnityMinusDiracQuasiInverseDaggerMatrixMultiplication(p2, output, true);    
    } else {
      executeDiracUnityMinusDiracQuasiInverseMatrixMultiplication(p2, output, true);    
    }
  } else {
    if (daggered) {
      executeDiracUnityMinusDiracQuasiInverseDaggerMatrixMultiplication(p2, output, true);
    } else {
      executeDiracUnityMinusDiracQuasiInverseMatrixMultiplication(p2, output, true);
    }
    executePiModeRemoverOperator(p2, InputVectorFourierTransform, true);
  }
    
  performFFT(InputVectorFourierTransform, p5, ExtremeFFT4D_Backward);         

  if (daggered) {
    if (Parameter_MassSplit==1) {
      perform_yBDagger(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, fac*Parameter_yN/vol, phiField, p5, p5);
    } else {
      perform_yBsplitDagger(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, fac*Parameter_yN/vol, Parameter_MassSplit, phiField, p5, p5);    
    }
  } else {
    if (Parameter_MassSplit==1) {  
      perform_yB(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, fac*Parameter_yN/vol, phiField, p5, p5);
    } else {
      perform_yBsplit(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, fac*Parameter_yN/vol, Parameter_MassSplit, phiField, p5, p5);    
    }
  }

  performFFT(p5, InputVectorFourierTransform, ExtremeFFT4D_Forward);       

  executePiModeRemoverOperator(InputVectorFourierTransform, InputVectorFourierTransform, true);

  Complex alpha(fac,0);
  SSE_ComplexVectorAddition(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, alpha, p3, p4);

  if (useRPrec) executeRPreconditionerMatrixMultiplication(p4, p4, false, false, true);  

  if (!inFourierSpace) performFFT(p4, output, ExtremeFFT4D_Backward);  
}


void FermionMatrixOperations::executeFermionMatrixMultiplication(Complex* input, Complex* output, double* phiField, bool ScalarsForInvSolver, Complex* s1, Complex* s2, int preconLevel, int QPrec) {
  performFFT(input, InputVectorFourierTransform, ExtremeFFT4D_Forward);

  Complex* source;
  if ((Preconditioner_Usage) && (preconLevel>0)) {
    if (preconLevel==1) {
      executeFermionMatrixStaticInverseMultiplication(InputVectorFourierTransform, InputVectorFourierTransform, false, true);
    }
    if (preconLevel==2) {
      executeFermionMatrixStaticInverseMultiplication(InputVectorFourierTransform, InputVectorFourierTransform, true, true);
    }

    performFFT(InputVectorFourierTransform, output, ExtremeFFT4D_Backward);
    
    source = output;
    SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, InputVectorFourierTransform, InputVectorFourierTransform, 
     NeubergerDiracOpApplicationSinPStdData, NeubergerDiracOpApplicationEWDataRescaled);
  } else {
    source = input;
    SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, InputVectorFourierTransform, InputVectorFourierTransform, 
     NeubergerDiracOpApplicationSinPStdData, NeubergerDiracOpApplicationEWData);
  }
    
  performFFT(InputVectorFourierTransform, InterimVectorFourierTransform, ExtremeFFT4D_Backward);

  if ((ScalarsForInvSolver) && ((!QPreconditioner_Usage) || (QPrec!=1))) {
    if (Parameter_MassSplit == 1) { 
      performf_YBD_2rhoD_2rhoD_AndScalarProducts(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3,
       Parameter_yN, 2.0*Parameter_rho, Parameter_ExplicitMass, source, InterimVectorFourierTransform, phiField, output, Inverse_rest,
       Inverse_p, *s1, *s2);
     } else {
      performf_YBsplitD_2rhoD_2rhoD_AndScalarProducts(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3,
       Parameter_yN, Parameter_MassSplit, 2.0*Parameter_rho, Parameter_ExplicitMass, source, InterimVectorFourierTransform, phiField, output, Inverse_rest,
       Inverse_p, *s1, *s2);     
     }
  } else {
    if (Parameter_MassSplit == 1) {    
      performf_YBD_2rhoD_2rhoD(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3,
       Parameter_yN, 2.0*Parameter_rho, Parameter_ExplicitMass, source, InterimVectorFourierTransform, phiField, output);
    } else {
       performf_YBsplitD_2rhoD_2rhoD(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3,
       Parameter_yN, Parameter_MassSplit, 2.0*Parameter_rho, Parameter_ExplicitMass, source, InterimVectorFourierTransform, phiField, output);    
    }
  }

  if ((QPreconditioner_Usage) && (QPrec==1)) {
    executeQPreconditionerMatrixMultiplication(output, output, false, false);
 
    if (ScalarsForInvSolver) {
      SSE_ComplexScalarProduct(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, Inverse_p, Inverse_rest, s1[0]); 
      SSE_ComplexScalarProduct(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, Inverse_p, output, s2[0]); 
    }
  }
}


void FermionMatrixOperations::executeFermionDaggerMatrixMultiplication(Complex* input, Complex* output, double* phiField, int QPrec) {
  Complex* source = input;
  if ((QPreconditioner_Usage) && (QPrec == 1)) {
    executeQPreconditionerMatrixMultiplication(input, InterimVectorFourierTransform, false, false);
    source = InterimVectorFourierTransform;    
  }

  if (Parameter_ExplicitMass==0) {
    if (Parameter_MassSplit == 1) {      
      SSE_Performf_YB_2rho_AndCopyToOutput(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3,
       Parameter_yN, 2.0*Parameter_rho, source, phiField, InputVectorFourierTransform, output);
    } else {
      performf_YBsplit_2rho_AndCopyToOutput(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3,
       Parameter_yN, Parameter_MassSplit, 2.0*Parameter_rho, Parameter_ExplicitMass, source, phiField, InputVectorFourierTransform, output);  
    }
  } else {
    performf_YBsplit_2rho_AndCopyToOutput(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3,
     Parameter_yN, Parameter_MassSplit, 2.0*Parameter_rho, Parameter_ExplicitMass, source, phiField, InputVectorFourierTransform, output);  
  } 

  performFFT(InputVectorFourierTransform, InterimVectorFourierTransform, ExtremeFFT4D_Forward);
    
  SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, InterimVectorFourierTransform, InterimVectorFourierTransform, 
   NeubergerDiracOpApplicationSinPDagData, NeubergerDiracOpApplicationEWData);

  performFFT(InterimVectorFourierTransform, InputVectorFourierTransform, ExtremeFFT4D_Backward);

  SSE_ComplexVectorAdditionSPECIAL2(oneDimSizeL0,oneDimSizeL1,oneDimSizeL2,oneDimSizeL3,xtraSize1,xtraSize2,xtraSize3, 2*Parameter_rho, InputVectorFourierTransform, output);
}


void FermionMatrixOperations::executeFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(Complex* input, Complex* output, double* phiField, bool useRPrec, bool inFourierSpace) {
  Complex* p1 = NULL;
  Complex* p2 = NULL;
  Complex* p3 = NULL;
  Complex* p4 = NULL;  
  Complex* p5 = NULL;  
  double vol = oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3;
  double fac = 1.0;
  useRPrec = useRPrec && RPreconditioner_Usage;
  
  if ((xFFTusage) && (inFourierSpace) && (useRPrec)) {
    p1 = input;
    p2 = InputVectorFourierTransform;
    p3 = InputVectorFourierTransform;
    p4 = output;
    p5 = InterimVectorFourierTransform;
  } 
  if ((xFFTusage) && (inFourierSpace) && (!useRPrec)) {
    p2 = input;
    p3 = InputVectorFourierTransform;
    p4 = output;
    p5 = InterimVectorFourierTransform;
  } 
  if ((xFFTusage) && (!inFourierSpace) && (useRPrec)) {
    p1 = InputVectorFourierTransform;
    p2 = InputVectorFourierTransform;
    p3 = output;
    p4 = InputVectorFourierTransform;
    p5 = InterimVectorFourierTransform;
    fac = 1.0 / vol;
  } 
  if ((xFFTusage) && (!inFourierSpace) && (!useRPrec)) {
    p2 = InputVectorFourierTransform;
    p3 = output;
    p4 = InputVectorFourierTransform;
    p5 = InterimVectorFourierTransform;
    fac = 1.0 / vol;
  } 
  if ((!xFFTusage) && (inFourierSpace) && (useRPrec)) {
    p1 = input;
    p2 = InputVectorFourierTransform;
    p3 = InputVectorFourierTransform;
    p4 = output;
    p5 = InputVectorFourierTransform;
  } 
  if ((!xFFTusage) && (inFourierSpace) && (!useRPrec)) {
    p2 = input;
    p3 = InputVectorFourierTransform;
    p4 = output;
    p5 = InputVectorFourierTransform;
  } 
  if ((!xFFTusage) && (!inFourierSpace) && (useRPrec)) {
    p1 = InputVectorFourierTransform;
    p2 = InputVectorFourierTransform;
    p3 = InputVectorFourierTransform;
    p4 = output;
    p5 = InputVectorFourierTransform;
    fac = 1.0 / vol;
  } 
  if ((!xFFTusage) && (!inFourierSpace) && (!useRPrec)) {
    p2 = InputVectorFourierTransform;
    p3 = InputVectorFourierTransform;
    p4 = output;
    p5 = InputVectorFourierTransform;
    fac = 1.0 / vol;
  } 
  
  if (!inFourierSpace) performFFT(input, InputVectorFourierTransform, ExtremeFFT4D_Forward);  
    
  if (useRPrec) executeRPreconditionerMatrixMultiplication(p1, InputVectorFourierTransform, false, false, true);  
    
  if (p2 == input) {
    executePiModeRemoverOperator(p2, InputVectorFourierTransform, true);  
    executeDiracUnityMinusDiracQuasiInverseDaggerMatrixMultiplication(p2, output, true);    
  } else {
    executeDiracUnityMinusDiracQuasiInverseDaggerMatrixMultiplication(p2, output, true);
    executePiModeRemoverOperator(p2, InputVectorFourierTransform, true);
  }
    
  performFFT(InputVectorFourierTransform, p5, ExtremeFFT4D_Backward);         

  if (Parameter_MassSplit == 1) {
    perform_yBDagger(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, fac*Parameter_yN/vol, phiField, p5, p5);
  } else {
    perform_yBsplitDagger(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, fac*Parameter_yN/vol, Parameter_MassSplit, phiField, p5, p5);  
  }

  performFFT(p5, InputVectorFourierTransform, ExtremeFFT4D_Forward);       

  executePiModeRemoverOperator(InputVectorFourierTransform, InputVectorFourierTransform, true);

  Complex alpha(fac,0);
  SSE_ComplexVectorAddition(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, alpha, output, InputVectorFourierTransform);

  if (useRPrec) executeRPreconditionerMatrixMultiplication(InputVectorFourierTransform, InputVectorFourierTransform, true, false, true);  

  executeDiracUnityMinusDiracQuasiInverseMatrixMultiplication(InputVectorFourierTransform, output, true);    

  executePiModeRemoverOperator(InputVectorFourierTransform, InputVectorFourierTransform, true);
  
  performFFT(InputVectorFourierTransform, p5, ExtremeFFT4D_Backward);         
  if (Parameter_MassSplit == 1) {  
    perform_yB(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, Parameter_yN/vol, phiField, p5, p5);
  } else {
    perform_yBsplit(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, Parameter_yN/vol, Parameter_MassSplit, phiField, p5, p5);  
  }
  performFFT(p5, InputVectorFourierTransform, ExtremeFFT4D_Forward);       

  executePiModeRemoverOperator(InputVectorFourierTransform, InputVectorFourierTransform, true);

  alpha.x = 1.0;
  SSE_ComplexVectorAddition(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, alpha, p3, p4);

  if (useRPrec) executeRPreconditionerMatrixMultiplication(p4, p4, false, false, true);  

  if (!inFourierSpace) performFFT(p4, output, ExtremeFFT4D_Backward);  
}


void FermionMatrixOperations::executeDistributedFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(DistributedMemoryObject* input, DistributedMemoryObject* output, DistributedMemoryObject* phiField, bool useRPrec, bool inFourierSpace) {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::executeDistributedFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication: threadedOps==NULL\n");
    exit(0);
  }
  DistributedMemoryObject* p1 = NULL;
  DistributedMemoryObject* p2 = NULL;
  DistributedMemoryObject* p3 = NULL;
  DistributedMemoryObject* p4 = NULL;  
  DistributedMemoryObject* p5 = NULL;  
  double vol = oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3;
  double fac = 1.0;
  useRPrec = useRPrec && RPreconditioner_Usage;
  
  if ((xFFTusage) && (inFourierSpace) && (useRPrec)) {
    p1 = input;
    p2 = DistributedInputVectorFourierTransform;
    p3 = DistributedInputVectorFourierTransform;
    p4 = output;
    p5 = DistributedInterimVectorFourierTransform;
  } 
  if ((xFFTusage) && (inFourierSpace) && (!useRPrec)) {
    p2 = input;
    p3 = DistributedInputVectorFourierTransform;
    p4 = output;
    p5 = DistributedInterimVectorFourierTransform;
  } 
  if ((xFFTusage) && (!inFourierSpace) && (useRPrec)) {
    p1 = DistributedInputVectorFourierTransform;
    p2 = DistributedInputVectorFourierTransform;
    p3 = output;
    p4 = DistributedInputVectorFourierTransform;
    p5 = DistributedInterimVectorFourierTransform;
    fac = 1.0 / vol;
  } 
  if ((xFFTusage) && (!inFourierSpace) && (!useRPrec)) {
    p2 = DistributedInputVectorFourierTransform;
    p3 = output;
    p4 = DistributedInputVectorFourierTransform;
    p5 = DistributedInterimVectorFourierTransform;
    fac = 1.0 / vol;
  } 
  if ((!xFFTusage) && (inFourierSpace) && (useRPrec)) {
    p1 = input;
    p2 = DistributedInputVectorFourierTransform;
    p3 = DistributedInputVectorFourierTransform;
    p4 = output;
    p5 = DistributedInputVectorFourierTransform;
  } 
  if ((!xFFTusage) && (inFourierSpace) && (!useRPrec)) {
    p2 = input;
    p3 = DistributedInputVectorFourierTransform;
    p4 = output;
    p5 = DistributedInputVectorFourierTransform;
  } 
  if ((!xFFTusage) && (!inFourierSpace) && (useRPrec)) {
    p1 = DistributedInputVectorFourierTransform;
    p2 = DistributedInputVectorFourierTransform;
    p3 = DistributedInputVectorFourierTransform;
    p4 = output;
    p5 = DistributedInputVectorFourierTransform;
    fac = 1.0 / vol;
  } 
  if ((!xFFTusage) && (!inFourierSpace) && (!useRPrec)) {
    p2 = DistributedInputVectorFourierTransform;
    p3 = DistributedInputVectorFourierTransform;
    p4 = output;
    p5 = DistributedInputVectorFourierTransform;
    fac = 1.0 / vol;
  } 
  
  if (!inFourierSpace) performDistributedFFT(input, DistributedInputVectorFourierTransform, ExtremeFFT4D_Forward);
    
  if (useRPrec) executeDistributedRPreconditionerMatrixMultiplication(p1, DistributedInputVectorFourierTransform, false, false, true);  
    
  if (p2 == input) {
    executeDistributedPiModeRemoverOperator(p2, DistributedInputVectorFourierTransform, true);  
    executeDistributedDiracUnityMinusDiracQuasiInverseDaggerMatrixMultiplication(p2, output, true);    
  } else {
    executeDistributedDiracUnityMinusDiracQuasiInverseDaggerMatrixMultiplication(p2, output, true);
    executeDistributedPiModeRemoverOperator(p2, DistributedInputVectorFourierTransform, true);
  }
    
  performDistributedFFT(DistributedInputVectorFourierTransform, p5, ExtremeFFT4D_Backward);         

  threadedOps->perform_yBDaggerOperatorMultiplication(p5, p5, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, fac*Parameter_yN/vol, Parameter_MassSplit, phiField);

  performDistributedFFT(p5, DistributedInputVectorFourierTransform, ExtremeFFT4D_Forward);       

  executeDistributedPiModeRemoverOperator(DistributedInputVectorFourierTransform, DistributedInputVectorFourierTransform, true);
  
  Complex alpha(fac,0);
  threadedOps->vectorAdditionOfFermionVectors(output, DistributedInputVectorFourierTransform, alpha, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);

  if (useRPrec) executeDistributedRPreconditionerMatrixMultiplication(DistributedInputVectorFourierTransform, DistributedInputVectorFourierTransform, true, false, true);  

  executeDistributedDiracUnityMinusDiracQuasiInverseMatrixMultiplication(DistributedInputVectorFourierTransform, output, true);    

  executeDistributedPiModeRemoverOperator(DistributedInputVectorFourierTransform, DistributedInputVectorFourierTransform, true);
  
  performDistributedFFT(DistributedInputVectorFourierTransform, p5, ExtremeFFT4D_Backward);         
  threadedOps->perform_yBOperatorMultiplication(p5, p5, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, Parameter_yN/vol, Parameter_MassSplit, phiField);
  performDistributedFFT(p5, DistributedInputVectorFourierTransform, ExtremeFFT4D_Forward);       

  executeDistributedPiModeRemoverOperator(DistributedInputVectorFourierTransform, DistributedInputVectorFourierTransform, true);

  alpha.x = 1.0;
  threadedOps->vectorAdditionOfFermionVectors(p3, p4, alpha, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);

  if (useRPrec) executeDistributedRPreconditionerMatrixMultiplication(p4, p4, false, false, true);  

  if (!inFourierSpace) performDistributedFFT(p4, output, ExtremeFFT4D_Backward); 
}

/**
* Input must be different from output.
*/
void FermionMatrixOperations::executeFermionMatrixFermionDaggerMatrixMultiplication(Complex* input, Complex* output, double* phiField, int preconLevel, int QPrec, bool inFourierSpace) {
  double VolumeNorm = 1.0 / (oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3);

  Complex* source = input;
  Complex* p2 = InputVectorFourierTransform;
  Complex* p3 = InterimVectorFourierTransform;

  if (!inFourierSpace) {
    performFFT(source, p2, ExtremeFFT4D_Forward);
    source = p2;
  }

  if ((QPreconditioner_Usage) && (QPrec == 1)) {
    executeQPreconditionerMatrixMultiplication(source, p2, false, true);
    source = p2;
  } 

  if ((!inFourierSpace) && ((!QPreconditioner_Usage) || (QPrec != 1))) {
    if (Parameter_MassSplit == 1) { 
      perform_yBDagger(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, Parameter_yN/VolumeNorm, phiField, input, p3);
    } else {
      perform_yBsplitDagger(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, Parameter_yN/VolumeNorm, Parameter_MassSplit, phiField, input, p3);    
    }
  } else {
    performFFT(source, p3, ExtremeFFT4D_Backward);
    if (Parameter_MassSplit == 1) {     
      perform_yBDagger(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, Parameter_yN, phiField, p3, p3);
    } else {
      perform_yBsplitDagger(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, Parameter_yN, Parameter_MassSplit, phiField, p3, p3);    
    }
  }
  performFFT(p3, output, ExtremeFFT4D_Forward);

  executeDiracDaggerMatrixMultiplication(output, p3, true);
  if (Parameter_ExplicitMass==0) {
    executeDiracDaggerMatrixMultiplication(source, p2, true);
  } else {
    executeDiracDaggerPlusExplicitMassMatrixMultiplication(source, p2, true);  //plus ExplicitMass
  }
    
  double f0 = VolumeNorm;
  double f1 = 2*Parameter_rho;
  if ((Preconditioner_Usage) && (preconLevel>0)) {
    double g1 = f1 / VolumeNorm;
    for (int I2=0; I2<vectorLengthXtrSize; I2++) {
      output[I2].x = (p3[I2].x - f1*output[I2].x) - g1*p2[I2].x;
      output[I2].y = (p3[I2].y - f1*output[I2].y) - g1*p2[I2].y;
    }
    executeFermionMatrixStaticInverseMultiplication(output, output, true, true);
  } else {
    for (int I2=0; I2<vectorLengthXtrSize; I2++) {
      output[I2].x = f0*(p3[I2].x - f1*output[I2].x) - f1*p2[I2].x;
      output[I2].y = f0*(p3[I2].y - f1*output[I2].y) - f1*p2[I2].y;
    }
  }
  executeDiracMatrixMultiplication(output, p3, true);
  if (Parameter_ExplicitMass==0) {
    for (int I2=0; I2<vectorLengthXtrSize; I2++) {
      output[I2].x = p3[I2].x - f1*output[I2].x;
      output[I2].y = p3[I2].y - f1*output[I2].y;
    }
  } else {
    double ox,oy;
    for (int I2=0; I2<vectorLengthXtrSize; I2++) {
      ox = output[I2].x;
      oy = output[I2].y;
      output[I2].x = p3[I2].x - f1*output[I2].x;
      output[I2].y = p3[I2].y - f1*output[I2].y;
      p3[I2].x += Parameter_ExplicitMass*ox; 
      p3[I2].y += Parameter_ExplicitMass*oy;
    }
  }
  performFFT(output, p2, ExtremeFFT4D_Backward);
  if (Parameter_MassSplit == 1) {       
    perform_yB(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, Parameter_yN, phiField, p2, p2);
  } else {
    perform_yBsplit(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, Parameter_yN, Parameter_MassSplit, phiField, p2, p2);  
  }
  performFFT(p2, output, ExtremeFFT4D_Forward);

  if (!inFourierSpace) {
    f0 *= VolumeNorm;
    f1 *= VolumeNorm;
    for (int I2=0; I2<vectorLengthXtrSize; I2++) {
      p3[I2].x = f0*output[I2].x - f1*p3[I2].x;
      p3[I2].y = f0*output[I2].y - f1*p3[I2].y;
    }
    if ((QPreconditioner_Usage) && (QPrec == 1)) {
      executeQPreconditionerMatrixMultiplication(p3, p3, false, true);   
    }
    performFFT(p3, output, ExtremeFFT4D_Backward);
  } else {
    for (int I2=0; I2<vectorLengthXtrSize; I2++) {
      output[I2].x = f0*output[I2].x - f1*p3[I2].x;
      output[I2].y = f0*output[I2].y - f1*p3[I2].y;
    }
    if ((QPreconditioner_Usage) && (QPrec == 1)) {
      executeQPreconditionerMatrixMultiplication(output, output, false, true);   
    }
  }
}


void FermionMatrixOperations::executeQPreconditionerMatrixMultiplication(Complex* input, Complex* output, bool inverse, bool inFourierSpace) {
  Complex* target;
  Complex* source;
  Complex* mulData;
  if (inFourierSpace) {
    source = input;
    target = output;
    if (inverse) {
      mulData = inverseQPreconditionerDiagonalData;
    } else {
      mulData = QPreconditionerDiagonalData;
    }
  } else {
    performFFT(input, InputVectorFourierTransform, ExtremeFFT4D_Forward);
  
    source = InputVectorFourierTransform;
    target = InputVectorFourierTransform;
    if (inverse) {
      mulData = inverseQPreconditionerDiagonalDataNormalized;    
    } else {
      mulData = QPreconditionerDiagonalDataNormalized;
    }
  }
  
       
  int I0,I1,I2,I3, i;
  int count = 0;
  int countXtr = 0;
  int xtrAdd1 = xtraSize1*8*(oneDimSizeL3+xtraSize3)*(oneDimSizeL2+xtraSize2);
  int xtrAdd2 = xtraSize2*8*(oneDimSizeL3+xtraSize3);
  int xtrAdd3 = xtraSize3*8;
  for (I0=0; I0<oneDimSizeL0; I0++) {
    for (I1=0; I1<oneDimSizeL1; I1++) {
      for (I2=0; I2<oneDimSizeL2; I2++) {
        for (I3=0; I3<oneDimSizeL3; I3++) {
          for (i=0; i<8; i++) {
	    //Imagin�rteil wird hier ignoriert, weil er sowieso Null ist.
            target[countXtr].x = source[countXtr].x * mulData[count].x;
            target[countXtr].y = source[countXtr].y * mulData[count].x;

            countXtr++;
          }
          count++;
	}
        countXtr += xtrAdd3;
      }
      countXtr += xtrAdd2;
    }
    countXtr += xtrAdd1;
  }

  if (!inFourierSpace) {
    performFFT(InputVectorFourierTransform, output, ExtremeFFT4D_Backward);
  }
}


void FermionMatrixOperations::executeRPreconditionerMatrixMultiplication(Complex* input, Complex* output, bool dob, bool inverse, bool inFourierSpace) {
  Complex* target;
  Complex* source;
  Complex* mulData;
  if (inFourierSpace) {
    source = input;
    target = output;
    if (inverse) {
      if (dob) {
        mulData = inverseRSQRPreconditionerDiagonalData;
      } else {
        mulData = inverseRPreconditionerDiagonalData;      
      }
    } else {
      if (dob) {
        mulData = RSQRPreconditionerDiagonalData;
      } else {
        mulData = RPreconditionerDiagonalData;      
      }
    }
  } else {
    performFFT(input, InputVectorFourierTransform, ExtremeFFT4D_Forward);
  
    source = InputVectorFourierTransform;
    target = InputVectorFourierTransform;
    if (inverse) {
      if (dob) {
        mulData = inverseRSQRPreconditionerDiagonalDataNormalized;    
      } else {
        mulData = inverseRPreconditionerDiagonalDataNormalized;          
      }
    } else {
      if (dob) {
        mulData = RSQRPreconditionerDiagonalDataNormalized;
      } else {
        mulData = RPreconditionerDiagonalDataNormalized;      
      }
    }
  }
  
       
  int I0,I1,I2,I3, i;
  int count = 0;
  int countXtr = 0;
  int xtrAdd1 = xtraSize1*8*(oneDimSizeL3+xtraSize3)*(oneDimSizeL2+xtraSize2);
  int xtrAdd2 = xtraSize2*8*(oneDimSizeL3+xtraSize3);
  int xtrAdd3 = xtraSize3*8;
  for (I0=0; I0<oneDimSizeL0; I0++) {
    for (I1=0; I1<oneDimSizeL1; I1++) {
      for (I2=0; I2<oneDimSizeL2; I2++) {
        for (I3=0; I3<oneDimSizeL3; I3++) {
          for (i=0; i<4; i++) {
	    //Imagin�rteil von mulData enth�t Faktor fr bottom-quark.
            target[countXtr].x = source[countXtr].x * mulData[count].x;
            target[countXtr].y = source[countXtr].y * mulData[count].x;
            countXtr++;
          }
          for (i=0; i<4; i++) {
	    //Imagin�rteil von mulData enth�t Faktor fr bottom-quark.
            target[countXtr].x = source[countXtr].x * mulData[count].y;
            target[countXtr].y = source[countXtr].y * mulData[count].y;
            countXtr++;
          }
	  
          count++;
	}
        countXtr += xtrAdd3;
      }
      countXtr += xtrAdd2;
    }
    countXtr += xtrAdd1;
  }

  if (!inFourierSpace) {
    performFFT(InputVectorFourierTransform, output, ExtremeFFT4D_Backward);
  }
}


void FermionMatrixOperations::executeDistributedRPreconditionerMatrixMultiplication(DistributedMemoryObject* input, DistributedMemoryObject* output, bool dob, bool inverse, bool inFourierSpace) {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::executeDistributedRPreconditionerMatrixMultiplication: threadedOps==NULL\n");
    exit(0);
  }
  DistributedMemoryObject* target;
  DistributedMemoryObject* source;
  DistributedMemoryObject* mulData;
  if (inFourierSpace) {
    source = input;
    target = output;
    if (inverse) {
      if (dob) {
        mulData = Distributed_inverseRSQRPreconditionerDiagonalData;
      } else {
        mulData = Distributed_inverseRPreconditionerDiagonalData;      
      }
    } else {
      if (dob) {
        mulData = Distributed_RSQRPreconditionerDiagonalData;
      } else {
        mulData = Distributed_RPreconditionerDiagonalData;      
      }
    }
  } else {
    performDistributedFFT(input, DistributedInputVectorFourierTransform, ExtremeFFT4D_Forward);  
  
    source = DistributedInputVectorFourierTransform;
    target = DistributedInputVectorFourierTransform;
    if (inverse) {
      if (dob) {
        mulData = Distributed_inverseRSQRPreconditionerDiagonalDataNormalized;    
      } else {
        mulData = Distributed_inverseRPreconditionerDiagonalDataNormalized;          
      }
    } else {
      if (dob) {
        mulData = Distributed_RSQRPreconditionerDiagonalDataNormalized;
      } else {
        mulData = Distributed_RPreconditionerDiagonalDataNormalized;      
      }
    }
  }

  threadedOps->perform_RPreconditionerTypeMatrixMultiplicationInFourierSpace(source, target, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, mulData);  
       
  if (!inFourierSpace) {
    performDistributedFFT(DistributedInputVectorFourierTransform, output, ExtremeFFT4D_Backward);  
  }
}


/**
* Auch im Fourier-Raum mit Normierungsfaktor 1/V multipliziert!!!
**/
void FermionMatrixOperations::executeFermionMatrixStaticInverseMultiplication(Complex* input, Complex* output, bool dob, bool inFourierSpace) {
  Complex* source;
  Complex* target;
  if (inFourierSpace) {
    source = input;
    target = output;
  } else {
    performFFT(input, InputVectorFourierTransform, ExtremeFFT4D_Forward);
  
    source = InputVectorFourierTransform;
    target = InputVectorFourierTransform;
  }
  
  
  if (dob) {
    SSE_ExecuteFermionMatrixStaticInverseMultiplicationInFourierSpace(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, source, target,
     NeubergerDiracOpApplicationSinPStdData, 
     NeubergerWithChiStaticMdagMFermionInverseOpApplicationSSEAuxData, 
     Index_PlusPiPiPiPiXtrSize);
  } else {
    Complex* c1 = createFermionVector();
    Complex* c2 = createFermionVector();

    SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, source, c1, 
     NeubergerDiracOpApplicationSinPStdData, NeubergerWithChiStaticFermionInverseOpApplicationMainDiagonalData);
    SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, source, c2, 
     NeubergerDiracOpApplicationSinPStdData, NeubergerWithChiStaticFermionInverseOpApplicationSubDiagonalData);
       
    int I0,I1,I2,I3, i;
    int count = 0;
    int countPi = 0;
    int N0half = oneDimSizeL0/2;
    int xtrAdd1 = xtraSize1*8*(oneDimSizeL3+xtraSize3)*(oneDimSizeL2+xtraSize2);
    int xtrAdd2 = xtraSize2*8*(oneDimSizeL3+xtraSize3);
    int xtrAdd3 = xtraSize3*8;
    for (I0=oneDimSizeL0-1; I0>=N0half; I0--) {
      countPi += Index_PlusPiPiPiPiXtrSize[4*(I0+0)];
      for (I1=oneDimSizeL1-1; I1>=0; I1--) {
        countPi += Index_PlusPiPiPiPiXtrSize[4*(I1+128)];
        for (I2=oneDimSizeL2-1; I2>=0; I2--) {
          countPi += Index_PlusPiPiPiPiXtrSize[4*(I2+256)];
          for (I3=oneDimSizeL3-1; I3>=0; I3--) {
            countPi += Index_PlusPiPiPiPiXtrSize[4*(I3+384)];
            countPi /= 16;
            for (i=0; i<8; i++) {
 	      target[count].x = c1[count].x + c2[countPi].x;
	      target[count].y = c1[count].y + c2[countPi].y;
              target[countPi].x = c1[countPi].x + c2[count].x;
              target[countPi].y = c1[countPi].y + c2[count].y;
              count++;
              countPi++;
            }
            countPi -= 8;
            countPi *= 16;
	  
            countPi -= Index_PlusPiPiPiPiXtrSize[4*(I3+384)];
   	  }
	  count += xtrAdd3;
          countPi -= Index_PlusPiPiPiPiXtrSize[4*(I2+256)];
        }
        count += xtrAdd2;
        countPi -= Index_PlusPiPiPiPiXtrSize[4*(I1+128)];
      }
      count += xtrAdd1;
      countPi -= Index_PlusPiPiPiPiXtrSize[4*(I0+0)];
    }
    destroyFermionVector(c1);
    destroyFermionVector(c2);
  }


  if (!inFourierSpace) {
    performFFT(InputVectorFourierTransform, output, ExtremeFFT4D_Backward);
  }
}


double FermionMatrixOperations::checkLGSsolutionAccuracy(Complex* solution, Complex* result, double* phiField) {
  executeFermionDaggerMatrixMultiplication(solution, Inverse_s, phiField, 0);
  executeFermionMatrixMultiplication(Inverse_s, Inverse_s, phiField, false, NULL, NULL, 2, 0);

  transformFromXtraSizeArray(Inverse_s, Inverse_s);
  transformFromXtraSizeArray(result, Inverse_p);
  
  int I;
  double diff = 0;
  for (I=0; I<vectorLength; I++) {
    diff += sqr(Inverse_s[I].x-Inverse_p[I].x);
    diff += sqr(Inverse_s[I].y-Inverse_p[I].y);
  }
  diff = sqrt(diff);
  return diff;
}


void FermionMatrixOperations::applyDistributedFermionMatrixMMDaggerChebyshevPolynomial(DistributedMemoryObject* input, DistributedMemoryObject* output, DistributedMemoryObject* phiField, GeneralChebyshevApproximation* chebyPoly, bool quasiHermiteanMode, bool inFourierSpace) {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::applyDistributedFermionMatrixMMDaggerChebyshevPolynomial: threadedOps==NULL\n");
    exit(0);
  }
  int coeffCount = chebyPoly->getCoeffCount();
  double alpha = chebyPoly->getAlpha();
  double beta = chebyPoly->getBeta();
  double* gamma = chebyPoly->getGammas();
  
  zeroDistributedFermionVector(output);

  if (coeffCount <= 0) {
    return;
  }
  if (coeffCount == 1) {
    double fac = gamma[0];
    Complex a1(fac,0);
    threadedOps->vectorAdditionOfFermionVectors(input, output, a1, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
    return;
  }
  if (coeffCount == 2) {
    if (quasiHermiteanMode) {
      executeDistributedFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(input, DistributedInverse_b, phiField, true, inFourierSpace);
    } else {
      printf("ERROR in FermionMatrixOperations::applyDistributedFermionMatrixMMDaggerChebyshevPolynomial: Distributed version of executeFermionMatrixFermionDaggerMatrixMultiplication not implemented yet.\n");
      exit(0);
    }
    double fac = gamma[0]+0.5*beta*gamma[1];
    double fac2 = 0.5*alpha*gamma[1];    
    Complex a1(fac,0);
    Complex a2(fac2,0);
    
    threadedOps->vectorAdditionOfFermionVectors(input, output, a1, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
    threadedOps->vectorAdditionOfFermionVectors(DistributedInverse_b, output, a2, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
    return;
  }

  //Mehr als 2 Koeffizienten
  DistributedMemoryObject* tim1 = DistributedInverse_rest;
  DistributedMemoryObject* ti = DistributedInverse_s;
  DistributedMemoryObject* tip1 = DistributedInverse_p;  
  
  threadedOps->copyFermionVector(input, tim1, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
  if (quasiHermiteanMode) {
    executeDistributedFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(input, DistributedInverse_b, phiField, true, inFourierSpace);
  } else {
    printf("ERROR in FermionMatrixOperations::applyDistributedFermionMatrixMMDaggerChebyshevPolynomial: Distributed version of executeFermionMatrixFermionDaggerMatrixMultiplication not implemented yet.\n");
    exit(0);
  }
  double fac = gamma[0]+0.5*beta*gamma[1];
  double fac2 = 0.5*alpha*gamma[1];      
  double fac3 = 0.5*beta;
  double fac4 = 0.5*alpha;      
  Complex a1(fac,0);
  Complex a2(fac2,0);
  Complex a3(fac3,0);
  Complex a4(fac4,0);
  Complex aAlpha(alpha,0);
  Complex aBeta(beta,0);
  Complex aM1(-1,0);
  
  
  threadedOps->multiplyFermionVectorWithComplexNumber(input, ti, a3, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
  threadedOps->vectorAdditionOfFermionVectors(DistributedInverse_b, ti, a4, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
  threadedOps->vectorAdditionOfFermionVectors(input, output, a1, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
  threadedOps->vectorAdditionOfFermionVectors(DistributedInverse_b, output, a2, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
  
  for (int I=2; I<coeffCount; I++) {
    if (quasiHermiteanMode) {
      executeDistributedFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(ti, tip1, phiField, true, inFourierSpace);
    } else {
      printf("ERROR in FermionMatrixOperations::applyDistributedFermionMatrixMMDaggerChebyshevPolynomial: Distributed version of executeFermionMatrixFermionDaggerMatrixMultiplication not implemented yet.\n");
      exit(0);
    }
    
    threadedOps->multiplyFermionVectorWithComplexNumber(tip1, tip1, aAlpha, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
    threadedOps->vectorAdditionOfFermionVectors(ti, tip1, aBeta, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
    threadedOps->vectorAdditionOfFermionVectors(tim1, tip1, aM1, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);

    DistributedMemoryObject* dummy = tim1;
    tim1 = ti;
    ti = tip1;
    tip1 = dummy;

    Complex compFac(gamma[I],0);
    threadedOps->vectorAdditionOfFermionVectors(ti, output, compFac, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
  }
}


void FermionMatrixOperations::applyFermionMatrixMMDaggerChebyshevPolynomial(Complex* input, Complex* output, double* phiField, GeneralChebyshevApproximation* chebyPoly, bool quasiHermiteanMode, bool inFourierSpace) {
  int coeffCount = chebyPoly->getCoeffCount();
  double alpha = chebyPoly->getAlpha();
  double beta = chebyPoly->getBeta();
  double* gamma = chebyPoly->getGammas();
  
  if (coeffCount <= 0) {
    for (int I=0; I<vectorLengthXtrSize; I++) {
      output[I].x = 0;
      output[I].y = 0;      
    }
    return;
  }
  if (coeffCount == 1) {
    double fac = gamma[0];
    for (int I=0; I<vectorLengthXtrSize; I++) {
      output[I].x = fac*input[I].x;
      output[I].y = fac*input[I].y;      
    }
    return;
  }
  if (coeffCount == 2) {
    if (quasiHermiteanMode) {
      executeFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(input, Inverse_b, phiField, true, inFourierSpace);
    } else {
      executeFermionMatrixFermionDaggerMatrixMultiplication(input, Inverse_b, phiField, 2, 1, inFourierSpace);    
    }
    double fac = gamma[0]+0.5*beta*gamma[1];
    double fac2 = 0.5*alpha*gamma[1];    
    for (int I=0; I<vectorLengthXtrSize; I++) {
      output[I].x = fac*input[I].x + fac2*Inverse_b[I].x;
      output[I].y = fac*input[I].y + fac2*Inverse_b[I].y;      
    }
    return;
  }

  //Mehr als 2 Koeffizienten
  Complex* tim1 = Inverse_rest;
  Complex* ti = Inverse_s;
  Complex* tip1 = Inverse_p;  
  
  SSE_ZCopy(vectorLengthXtrSize, input, 1, tim1, 1);
  if (quasiHermiteanMode) {
    executeFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(input, Inverse_b, phiField, true, inFourierSpace);
  } else {
    executeFermionMatrixFermionDaggerMatrixMultiplication(input, Inverse_b, phiField, 2, 1, inFourierSpace);    
  }
  double fac = gamma[0]+0.5*beta*gamma[1];
  double fac2 = 0.5*alpha*gamma[1];      
  double fac3 = 0.5*beta;
  double fac4 = 0.5*alpha;      
  for (int I=0; I<vectorLengthXtrSize; I++) {
    ti[I].x = fac3*input[I].x + fac4*Inverse_b[I].x;
    ti[I].y = fac3*input[I].y + fac4*Inverse_b[I].y;      
    output[I].x = fac*input[I].x + fac2*Inverse_b[I].x;
    output[I].y = fac*input[I].y + fac2*Inverse_b[I].y;      
  }
  
  for (int I=2; I<coeffCount; I++) {
    if (quasiHermiteanMode) {
      executeFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(ti, tip1, phiField, true, inFourierSpace);
    } else {
      executeFermionMatrixFermionDaggerMatrixMultiplication(ti, tip1, phiField, 2, 1, inFourierSpace);    
    }
    for (int I2=0; I2<vectorLengthXtrSize; I2++) {
      tip1[I2].x = alpha*tip1[I2].x + beta*ti[I2].x - tim1[I2].x;
      tip1[I2].y = alpha*tip1[I2].y + beta*ti[I2].y - tim1[I2].y;
    }
    Complex* dummy = tim1;
    tim1 = ti;
    ti = tip1;
    tip1 = dummy;

    Complex compFac(gamma[I],0);
    SSE_ComplexVectorAddition(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, compFac, ti, output); 
  }
}


bool FermionMatrixOperations::applyDistributedFermionMatrixMMDaggerFunction(DistributedMemoryObject* input, DistributedMemoryObject* output, DistributedMemoryObject* phiField, double (*func)(double x),double TOL, int maxIter, int& neededIter, int auxVeccount, DistributedMemoryObject** auxVecs, bool quasiHermiteanMode, bool inFourierSpace) {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::applyDistributedFermionMatrixMMDaggerFunction: threadedOps==NULL\n");
    exit(0);
  }
  if (LogLevel>3) printf("Solving Distributed-MMdag-Function linear set of equations with %d auxVectors...\n",auxVeccount);
  if (maxIter<=0) maxIter = 10000;

  double error = 0;
  DistributedMemoryObject* vec_x = DistributedInverse_b;  
  DistributedMemoryObject* vec_q = DistributedInverse_rest;
  if (auxVeccount>0) vec_q = auxVecs[0];
  DistributedMemoryObject* vec_qold = DistributedInverse_p;
  DistributedMemoryObject* vec_v = DistributedInverse_s;
  double* alpha = new double[maxIter+2];
  double* beta = new double[maxIter+2];
  double* r = new double[maxIter+2];
  double normFac = 1.0 / (oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3);

  beta[0] = 0;
  Complex dummyComp;
  r[0] = 0;

  threadedOps->zeroFermionVector(vec_qold, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);

  threadedOps->vectorNormOfFermionVector(input, r[1], oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);

  if (inFourierSpace) {
    r[1] = sqrt(1.0 / r[1]);
    Complex a1(r[1], 0);
    threadedOps->multiplyFermionVectorWithComplexNumber(input, vec_q, a1, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
  } else {
    r[1] = sqrt(1.0 / (r[1]*normFac));
    performDistributedFFT(input, vec_q, ExtremeFFT4D_Forward);
    Complex a1(normFac * r[1], 0);
    threadedOps->multiplyFermionVectorWithComplexNumber(vec_q, vec_q, a1, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
  }

  neededIter = 0;
  while (neededIter<maxIter) {
    neededIter++;
    if (neededIter<auxVeccount) {
      vec_v = auxVecs[neededIter];
    } else {
      if ((DistributedInverse_p!=vec_q) && (DistributedInverse_p!=vec_qold)) vec_v = DistributedInverse_p;
      if ((DistributedInverse_rest!=vec_q) && (DistributedInverse_rest!=vec_qold)) vec_v = DistributedInverse_rest;
      if ((DistributedInverse_s!=vec_q) && (DistributedInverse_s!=vec_qold)) vec_v = DistributedInverse_s;      
    }
    if (quasiHermiteanMode) {
      executeDistributedFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(vec_q, vec_v, phiField, true, true);
    } else {    
      printf("ERROR in FermionMatrixOperations::applyDistributedFermionMatrixMMDaggerFunction: Distributed version of executeFermionMatrixFermionDaggerMatrixMultiplication not implemented yet.\n");
      exit(0);
    }

    threadedOps->scalarProductOfFermionVectors(vec_q, vec_v, dummyComp, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
    alpha[neededIter] = dummyComp.x;

    Complex a1(-alpha[neededIter], 0);
    Complex a2(-beta[neededIter-1], 0);
    threadedOps->vectorAdditionOfFermionVectors(vec_q, vec_v, a1, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
    threadedOps->vectorAdditionOfFermionVectors(vec_qold, vec_v, a2, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);

    threadedOps->vectorNormOfFermionVector(vec_v, beta[neededIter], oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
    beta[neededIter] = sqrt(beta[neededIter]);

    vec_qold = vec_q;
    vec_q = vec_v;
    double fac = 1.0 / beta[neededIter];
    Complex a3(fac, 0);
    threadedOps->multiplyFermionVectorWithComplexNumber(vec_v, vec_q, a3, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);

    r[neededIter+1] = -(r[neededIter]*alpha[neededIter] + r[neededIter-1]*beta[neededIter-1]) / beta[neededIter];

    error = fabs(r[1]) / fabs(r[neededIter+1]);

    if (LogLevel>4) printf("...ERROR: %f\n", error);

    if (error < TOL) break;
  }

  if (error >= TOL) {
    printf("solveDistributedFermionMatrixMMDaggerFunctionLGS: Lanczos-iteration did not converge within %d steps!!! ERROR!!!\n",neededIter);
    exit(0);
  }

  char compz = 'I';
  long int n = neededIter;
  double* d = new double[neededIter];
  for (int I=0; I<neededIter; I++) d[I] = alpha[I+1];
  double* e = new double[neededIter];
  for (int I=0; I<neededIter-1; I++) e[I] = beta[I+1];
  double* z = new double[neededIter*neededIter];
  long int ldz = neededIter;
  double* work = new double[2*neededIter-2];
  long int info = 0;

  dsteqr_(&compz, &n, d, e, z, &ldz, work, &info);

  if (info != 0) {
    printf("solveDistributedFermionMatrixMMDaggerFunctionLGS: LAPACL-routine dsteqr did not work!!! ERRROR %ld!!!\n", info);
    exit(0);
  }

  int pos = 0;
  double fac = 1.0 / r[1];
  for (int I=0; I<neededIter; I++) {
    if (func != NULL) {
      e[I] = fac*z[pos] * ((*func)(d[I]));
    } else {
      e[I] = fac*z[pos] / sqrt(d[I]);    
    }
    pos += neededIter;
  }
  double* y = new double[neededIter];
  for (int I=0; I<neededIter; I++) {
    y[I] = 0;
    pos = I;
    for (int I2=0; I2<neededIter; I2++) {
      y[I] += z[pos] * e[I2];
      pos += neededIter;
    }
  }
  

  //Build result vector vec_x from saved vec_q if available
  int iterStart = 1;
  if (auxVeccount>0) {
    if (inFourierSpace) {
      vec_x = output;
    }  
    threadedOps->zeroFermionVector(vec_x, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
  
    for (int iter=1; iter<=neededIter; iter++) {
      if (iter>=auxVeccount) break;
      Complex compFac(y[iter-1], 0);
      threadedOps->vectorAdditionOfFermionVectors(auxVecs[iter-1], vec_x, compFac, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
    }
    iterStart = auxVeccount;

    vec_q = auxVecs[auxVeccount-1];
    if (auxVeccount>1) {
      vec_qold = auxVecs[auxVeccount-2];
    } else {
      vec_qold = DistributedInverse_p;
      threadedOps->zeroFermionVector(vec_qold, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
    }
  } else {
    threadedOps->zeroFermionVector(vec_qold, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
    if (inFourierSpace) {
      Complex a4(r[1], 0);
      threadedOps->multiplyFermionVectorWithComplexNumber(input, vec_q, a4, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
    } else {
      performDistributedFFT(input, vec_q, ExtremeFFT4D_Forward);
      Complex a4(normFac * r[1], 0);
      threadedOps->multiplyFermionVectorWithComplexNumber(vec_q, vec_q, a4, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
    }
    if (inFourierSpace) {
      vec_x = output;
    }  
    threadedOps->zeroFermionVector(vec_x, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
  }
  

  for (int iter=iterStart; iter<=neededIter; iter++) {
    Complex compFac(y[iter-1], 0);
    threadedOps->vectorAdditionOfFermionVectors(vec_q, vec_x, compFac, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);

    if ((DistributedInverse_p!=vec_q) && (DistributedInverse_p!=vec_qold)) vec_v = DistributedInverse_p;
    if ((DistributedInverse_rest!=vec_q) && (DistributedInverse_rest!=vec_qold)) vec_v = DistributedInverse_rest;
    if ((DistributedInverse_s!=vec_q) && (DistributedInverse_s!=vec_qold)) vec_v = DistributedInverse_s;      

    if (quasiHermiteanMode) {
      executeDistributedFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(vec_q, vec_v, phiField, true, true);
    } else {    
      printf("ERROR in FermionMatrixOperations::applyDistributedFermionMatrixMMDaggerFunction: Distributed version of executeFermionMatrixFermionDaggerMatrixMultiplication not implemented yet.\n");
      exit(0);
    }

    Complex a5(-alpha[iter], 0);
    Complex a6(-beta[iter-1], 0);
    threadedOps->vectorAdditionOfFermionVectors(vec_q, vec_v, a5, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
    threadedOps->vectorAdditionOfFermionVectors(vec_qold, vec_v, a6, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);

    vec_qold = vec_q;
    vec_q = vec_v;
    double fac = 1.0 / beta[iter];
    Complex a7(fac, 0);
    threadedOps->multiplyFermionVectorWithComplexNumber(vec_v, vec_q, a7, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
  }

  if (!inFourierSpace) {
    performDistributedFFT(vec_x, output, ExtremeFFT4D_Backward);
  }

  delete[] y;
  delete[] d;
  delete[] e;
  delete[] z;
  delete[] work;
  if (LogLevel>3) printf("... ready after %d recursions with error %1.15f.\n",neededIter,error);
  if ((LogLevel>0) && (neededIter>vectorLength)) printf("ALARM: More iterations than matrix size for MMDagger-g inverse calculation (%d)!!!\n",neededIter);
  delete[] alpha;
  delete[] beta;
  delete[] r;
  return ((error <= TOL) && (info==0));
}


bool FermionMatrixOperations::applyFermionMatrixMMDaggerFunction(Complex* input, Complex* output, double* phiField, double (*func)(double x),double TOL, int maxIter, int& neededIter, int auxVeccount, Complex** auxVecs, bool quasiHermiteanMode, bool inFourierSpace) {
  if (LogLevel>3) printf("Solving MMdag-Function linear set of equations with %d auxVectors...\n",auxVeccount);
  if (maxIter<=0) maxIter = 10000;

  double error = 0;
  Complex* vec_x = Inverse_b;  
  Complex* vec_q = Inverse_rest;
  if (auxVeccount>0) vec_q = auxVecs[0];
  Complex* vec_qold = Inverse_p;
  Complex* vec_v = Inverse_s;
  double* alpha = new double[maxIter+2];
  double* beta = new double[maxIter+2];
  double* r = new double[maxIter+2];
  double normFac = 1.0 / (oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3);

  beta[0] = 0;
  Complex dummyComp;
  r[0] = 0;

  for (int I=0; I<vectorLengthXtrSize; I++) {
    vec_qold[I].x = 0;
    vec_qold[I].y = 0;
  }

  SSE_ComplexSquareNorm(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, input, r[1]);

  if (inFourierSpace) {
    r[1] = sqrt(1.0 / r[1]);
    for (int I=0; I<vectorLengthXtrSize; I++) {
      vec_q[I].x = input[I].x * r[1];
      vec_q[I].y = input[I].y * r[1];
    }  
  } else {
    r[1] = sqrt(1.0 / (r[1]*normFac));
    performFFT(input, vec_q, ExtremeFFT4D_Forward);
    for (int I=0; I<vectorLengthXtrSize; I++) {
      vec_q[I].x *= normFac * r[1];
      vec_q[I].y *= normFac * r[1];
    }
  }

  neededIter = 0;
  while (neededIter<maxIter) {
    neededIter++;
    if (neededIter<auxVeccount) {
      vec_v = auxVecs[neededIter];
    } else {
      if ((Inverse_p!=vec_q) && (Inverse_p!=vec_qold)) vec_v = Inverse_p;
      if ((Inverse_rest!=vec_q) && (Inverse_rest!=vec_qold)) vec_v = Inverse_rest;
      if ((Inverse_s!=vec_q) && (Inverse_s!=vec_qold)) vec_v = Inverse_s;      
    }
    if (quasiHermiteanMode) {
      executeFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(vec_q, vec_v, phiField, true, true);
    } else {    
      executeFermionMatrixFermionDaggerMatrixMultiplication(vec_q, vec_v, phiField, 2, 1, true);
    }

    SSE_ComplexScalarProduct(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, vec_q, vec_v, dummyComp);
    alpha[neededIter] = dummyComp.x;

    for (int I=0; I<vectorLengthXtrSize; I++) {
      vec_v[I].x = vec_v[I].x - alpha[neededIter]*vec_q[I].x - beta[neededIter-1]*vec_qold[I].x;
      vec_v[I].y = vec_v[I].y - alpha[neededIter]*vec_q[I].y - beta[neededIter-1]*vec_qold[I].y;
    }
    SSE_ComplexSquareNorm(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, vec_v, beta[neededIter]);
    beta[neededIter] = sqrt(beta[neededIter]);

    vec_qold = vec_q;
    vec_q = vec_v;
    double fac = 1.0 / beta[neededIter];
    for (int I=0; I<vectorLengthXtrSize; I++) {
      vec_q[I].x = fac * vec_v[I].x;
      vec_q[I].y = fac * vec_v[I].y;
    }

    r[neededIter+1] = -(r[neededIter]*alpha[neededIter] + r[neededIter-1]*beta[neededIter-1]) / beta[neededIter];

    error = fabs(r[1]) / fabs(r[neededIter+1]);

    if (LogLevel>4) printf("...ERROR: %f\n", error);

    if (error < TOL) break;
  }
#ifndef KRYLOVTESTMODUS
  if (error >= TOL) {
    printf("solveFermionMatrixMMDaggerFunctionLGS: Lanczos-iteration did not converge within %d steps!!! ERROR!!!\n",neededIter);
    exit(0);
  }
#else  
  KRYLOVERRROR = error;
#endif  

  char compz = 'I';
  long int n = neededIter;
  double* d = new double[neededIter];
  for (int I=0; I<neededIter; I++) d[I] = alpha[I+1];
  double* e = new double[neededIter];
  for (int I=0; I<neededIter-1; I++) e[I] = beta[I+1];
  double* z = new double[neededIter*neededIter];
  long int ldz = neededIter;
  double* work = new double[2*neededIter-2];
  long int info = 0;

  dsteqr_(&compz, &n, d, e, z, &ldz, work, &info);

  if (info != 0) {
    printf("solveFermionMatrixMMDaggerFunctionLGS: LAPACL-routine dsteqr did not work!!! ERRROR %ld!!!\n", info);
    exit(0);
  }

  int pos = 0;
  double fac = 1.0 / r[1];
  for (int I=0; I<neededIter; I++) {
    if (func != NULL) {
      e[I] = fac*z[pos] * ((*func)(d[I]));
    } else {
      e[I] = fac*z[pos] / sqrt(d[I]);    
    }
    pos += neededIter;
  }
  double* y = new double[neededIter];
  for (int I=0; I<neededIter; I++) {
    y[I] = 0;
    pos = I;
    for (int I2=0; I2<neededIter; I2++) {
      y[I] += z[pos] * e[I2];
      pos += neededIter;
    }
  }
  

  //Build result vector vec_x from saved vec_q if available
  int iterStart = 1;
  if (auxVeccount>0) {
    if (inFourierSpace) {
      vec_x = output;
    }  
    for (int I=0; I<vectorLengthXtrSize; I++) {
      vec_x[I].x = 0;
      vec_x[I].y = 0;
    }
  
    for (int iter=1; iter<=neededIter; iter++) {
      if (iter>=auxVeccount) break;
      Complex compFac(y[iter-1], 0);
      SSE_ComplexVectorAddition(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, compFac, auxVecs[iter-1], vec_x);
    }
    iterStart = auxVeccount;

    vec_q = auxVecs[auxVeccount-1];
    if (auxVeccount>1) {
      vec_qold = auxVecs[auxVeccount-2];
    } else {
      vec_qold = Inverse_p;
      for (int I=0; I<vectorLengthXtrSize; I++) {
        vec_qold[I].x = 0;
        vec_qold[I].y = 0;
      }
    }
  } else {
    for (int I=0; I<vectorLengthXtrSize; I++) {
      vec_qold[I].x = 0;
      vec_qold[I].y = 0;
    }
    if (inFourierSpace) {
      for (int I=0; I<vectorLengthXtrSize; I++) {
        vec_q[I].x = input[I].x * r[1];
        vec_q[I].y = input[I].y * r[1];
      }  
    } else {
      performFFT(input, vec_q, ExtremeFFT4D_Forward);
      for (int I=0; I<vectorLengthXtrSize; I++) { 
        vec_q[I].x *= normFac * r[1];
        vec_q[I].y *= normFac * r[1];
      }
    }
    if (inFourierSpace) {
      vec_x = output;
    }  
    for (int I=0; I<vectorLengthXtrSize; I++) {
      vec_x[I].x = 0;
      vec_x[I].y = 0;
    }
  }
  

  for (int iter=iterStart; iter<=neededIter; iter++) {
    Complex compFac(y[iter-1], 0);
    SSE_ComplexVectorAddition(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, compFac, vec_q, vec_x);

    if ((Inverse_p!=vec_q) && (Inverse_p!=vec_qold)) vec_v = Inverse_p;
    if ((Inverse_rest!=vec_q) && (Inverse_rest!=vec_qold)) vec_v = Inverse_rest;
    if ((Inverse_s!=vec_q) && (Inverse_s!=vec_qold)) vec_v = Inverse_s;      

    if (quasiHermiteanMode) {
      executeFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(vec_q, vec_v, phiField, true, true);
    } else {    
      executeFermionMatrixFermionDaggerMatrixMultiplication(vec_q, vec_v, phiField, 2, 1, true);
    }

    for (int I=0; I<vectorLengthXtrSize; I++) {
      vec_v[I].x = vec_v[I].x - alpha[iter]*vec_q[I].x - beta[iter-1]*vec_qold[I].x;
      vec_v[I].y = vec_v[I].y - alpha[iter]*vec_q[I].y - beta[iter-1]*vec_qold[I].y;
    }

    vec_qold = vec_q;
    vec_q = vec_v;
    double fac = 1.0 / beta[iter];
    for (int I=0; I<vectorLengthXtrSize; I++) {
      vec_q[I].x = fac * vec_v[I].x;
      vec_q[I].y = fac * vec_v[I].y;
    }
  }

  if (!inFourierSpace) {
    performFFT(vec_x, output, ExtremeFFT4D_Backward);
  }

  delete[] y;
  delete[] d;
  delete[] e;
  delete[] z;
  delete[] work;
  if (LogLevel>3) printf("... ready after %d recursions with error %1.15f.\n",neededIter,error);
  if ((LogLevel>0) && (neededIter>vectorLength)) printf("ALARM: More iterations than matrix size for MMDagger-g inverse calculation (%d)!!!\n",neededIter);
  delete[] alpha;
  delete[] beta;
  delete[] r;
  return ((error <= TOL) && (info==0));
}



bool FermionMatrixOperations::solveFermionMatrixLGS(Complex* input, Complex* output, double* phiField, double TOL, bool doubleFermionMatrix, bool quasiHermiteanMode, int maxIter, int& neededIter) {
  int iterCount = 0;
  Complex errorOld;
  Complex errorNew;
  Complex alpha;
  Complex alphaZ;
  Complex alphaN;
  double beta;

  if (LogLevel>3) printf("Solving linear set of equations (QHMode=%d)...\n", quasiHermiteanMode);

  if (doubleFermionMatrix) {
    performFFT(input, Inverse_b, ExtremeFFT4D_Forward);
    double normFac = 1.0 / (oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3);
    for (int I=0; I<vectorLengthXtrSize; I++) {
      Inverse_b[I].x *= normFac;
      Inverse_b[I].y *= normFac;
    }
  } else {  
    if (quasiHermiteanMode) {
      executeFermionQuasiHermiteanMatrixMultiplication(input, Inverse_b, phiField, true, true, false);
    } else {
      executeFermionDaggerMatrixMultiplication(input, Inverse_b, phiField, 0);
    }
  }

  //Startvector
  zeroFermionVector(output);
  SSE_ZCopy(vectorLengthXtrSize, Inverse_b, 1, Inverse_rest, 1);

  //Weitere Initialisierung
  SSE_ZCopy(vectorLengthXtrSize, Inverse_rest, 1, Inverse_p, 1);
  SSE_ComplexSquareNorm(oneDimSizeL0,oneDimSizeL1,oneDimSizeL2,oneDimSizeL3, xtraSize1,xtraSize2,xtraSize3, Inverse_rest, errorOld.x);
  if (LogLevel>5) errorOld.print();

  //Hauptschleife
  neededIter = 0;
  while ((sqrt(errorOld.x) > TOL) && ((neededIter<maxIter) || (maxIter<=0))) {
    neededIter++;
    iterCount++;
    
    if (doubleFermionMatrix) {
      if (quasiHermiteanMode) {
        executeFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(Inverse_p, Inverse_s, phiField, true, true);
      } else {      
        executeFermionMatrixFermionDaggerMatrixMultiplication(Inverse_p, Inverse_s, phiField, 2, 1, true);
      }
    } else {
      if (quasiHermiteanMode) {
        executeFermionQuasiHermiteanMatrixMultiplication(Inverse_p, Inverse_b, phiField, true, false, false);
        executeFermionQuasiHermiteanMatrixMultiplication(Inverse_b, Inverse_s, phiField, true, true, false);    
      } else {
        executeFermionMatrixMultiplication(Inverse_p, Inverse_b, phiField, false, NULL, NULL, 0, 0);
        executeFermionDaggerMatrixMultiplication(Inverse_b, Inverse_s, phiField, 0);    
        if ((Preconditioner_Usage) || (QPreconditioner_Usage)) {
          printf("Hier fehlt noch die Anwendung der P- bzw. Q- Matrizen!!!\n");      
          exit(0);     
        }
      }
    }
    SSE_ComplexScalarProduct(oneDimSizeL0,oneDimSizeL1,oneDimSizeL2,oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, Inverse_p, Inverse_rest, alphaZ);
    SSE_ComplexScalarProduct(oneDimSizeL0,oneDimSizeL1,oneDimSizeL2,oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, Inverse_p, Inverse_s, alphaN);
    
    alpha = alphaZ / alphaN;
    
//    cblas_zaxpy(vectorLength, &alpha, Inverse_p, 1, output, 1);
    SSE_ComplexVectorAddition(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1,xtraSize2,xtraSize3, alpha, Inverse_p, output);

    errorNew.y=0;
/*    errorNew.x=0;
    double dummy;
    for (I=0; I<vectorLength; I++) {
      dummy = alpha.x*Inverse_s[I].x - alpha.y*Inverse_s[I].y;
      Inverse_rest[I].y -= alpha.x*Inverse_s[I].y + alpha.y*Inverse_s[I].x;
      Inverse_rest[I].x -= dummy;
      errorNew.x += Inverse_rest[I].x*Inverse_rest[I].x + Inverse_rest[I].y*Inverse_rest[I].y;    
    }*/
    alpha.x *= -1.0;
    alpha.y *= -1.0;
    SSE_ComplexVectorAdditionWithSquaredNorm(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1,xtraSize2,xtraSize3, alpha, Inverse_s, Inverse_rest, errorNew.x);

    beta = errorNew.x / errorOld.x;    
/*    for (I=0; I<vectorLength; I++) {
      Inverse_p[I].x = beta*Inverse_p[I].x + Inverse_rest[I].x;
      Inverse_p[I].y = beta*Inverse_p[I].y + Inverse_rest[I].y;
    }*/
    SSE_ComplexVectorAdditionSPECIAL1(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1,xtraSize2,xtraSize3, beta, Inverse_rest, Inverse_p);

    errorOld.x = errorNew.x;
    if (LogLevel>5) printf("%1.15f + %1.15fi\n",errorOld.x,errorOld.y);
  }

  if (doubleFermionMatrix) {
    performFFT(output, Inverse_p, ExtremeFFT4D_Backward);
    SSE_ZCopy(vectorLengthXtrSize, Inverse_p, 1, output, 1);
  }

  if (LogLevel>3) printf("... ready after %d recursions with error %1.15f.\n",iterCount,sqrt(errorOld.x));
  if ((LogLevel>0) && (iterCount>vectorLength)) printf("ALARM: More iterations than matrix size for inverse calculation (%d)!!!\n",iterCount);
  return (sqrt(errorOld.x) <= TOL);
}


bool FermionMatrixOperations::executeFermionMatrixDaggerInverseMatrixMultiplication(Complex* input, Complex* output, double* phiField, double TOL, bool quasiHermiteanMode, int maxIter, int& neededIter) {
  bool success = false;
  
  if (quasiHermiteanMode) {
    if (RPreconditioner_Usage) {
      executeRPreconditionerMatrixMultiplication(input, output, false, false, false);
      executeFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(output, output, phiField, true, false);
    } else {
      executeFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(input, output, phiField, true, false);    
    }
    success = solveFermionMatrixLGS(output, output, (double*) phiField, TOL, true, quasiHermiteanMode, maxIter, neededIter);    
    if (RPreconditioner_Usage) {
      executeRPreconditionerMatrixMultiplication(output, output, false, false, false);
    }  
  } else {
    executeFermionMatrixMultiplication(input, output, (double*) phiField, false, NULL, NULL, 2, 1);
    success = solveFermionMatrixLGS(output, output, (double*) phiField, TOL, true, quasiHermiteanMode, maxIter, neededIter);    
    if (QPreconditioner_Usage) {
      executeQPreconditionerMatrixMultiplication(output, output, false, false);
    }  
  }
  
  return success;
}


void FermionMatrixOperations::executeDistributedYBOperatorMultiplication(DistributedMemoryObject* input, DistributedMemoryObject* output, double xtrFac, DistributedMemoryObject* phiObj) {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::executeDistributedYBOperatorMultiplication: threadedOps==NULL\n");
    exit(0);
  }
  threadedOps->perform_yBOperatorMultiplication(input, output, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtrFac*Parameter_yN, Parameter_MassSplit, phiObj);
}


void FermionMatrixOperations::executeDistributedYBDaggerOperatorMultiplication(DistributedMemoryObject* input, DistributedMemoryObject* output, double xtrFac, DistributedMemoryObject* phiObj) {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::executeDistributedYBDaggerOperatorMultiplication: threadedOps==NULL\n");
    exit(0);
  }
  threadedOps->perform_yBDaggerOperatorMultiplication(input, output, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtrFac*Parameter_yN, Parameter_MassSplit, phiObj);
}


void FermionMatrixOperations::executeDistributedMultiplicationVectorWithDerivativesOfB(DistributedMemoryObject* leftInput, DistributedMemoryObject* rightInput, DistributedMemoryObject* output) {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::executeDistributedMultiplicationVectorWithDerivativesOfB: threadedOps==NULL\n");
    exit(0);
  }
  threadedOps->perform_MultiplicationVectorWithDerivativesOfB(leftInput, rightInput, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, Parameter_MassSplit, output);
}


void FermionMatrixOperations::executeMultiplicationVectorWithDerivativesOfB(Complex* leftInput, Complex* rightInput, Complex* output) {
  int I0,I1,I2,I3;
  int inIndex = 0;
  int outIndex = 0;
  Complex g1,g2,g3,g4, h1, h2, h3, h4;
  int xtrAdd1 = xtraSize1*8*(oneDimSizeL3+xtraSize3)*(oneDimSizeL2+xtraSize2);
  int xtrAdd2 = xtraSize2*8*(oneDimSizeL3+xtraSize3);
  int xtrAdd3 = xtraSize3*8;

  if (Parameter_MassSplit == 1) {
    for (I0=oneDimSizeL0-1; I0>=0; I0--) {
      for (I1=oneDimSizeL1-1; I1>=0; I1--) {
        for (I2=oneDimSizeL2-1; I2>=0; I2--) {
          for (I3=oneDimSizeL3-1; I3>=0; I3--) {
            //Derivative with respect to Phi_0
            SSE_ComplexScalarProduct(1,1,1,1,0,0,0, &(leftInput[inIndex]), &(rightInput[inIndex]), output[outIndex]);

            calcGamma5VectorProduct(&(leftInput[inIndex]), &(rightInput[inIndex]), g1);
            calcGamma5VectorProduct(&(leftInput[inIndex]), &(rightInput[inIndex+4]), g2);
            calcGamma5VectorProduct(&(leftInput[inIndex+4]), &(rightInput[inIndex]), g3);
            calcGamma5VectorProduct(&(leftInput[inIndex+4]), &(rightInput[inIndex+4]), g4);
    
            output[outIndex+1] = ComplexI*(g2+g3);
            output[outIndex+2] = g2-g3;
            output[outIndex+3] = ComplexI*(g1-g4);
  
            outIndex += 4;
            inIndex += 8; 
 	  }
	  inIndex += xtrAdd3;
        }
        inIndex += xtrAdd2;
      } 
      inIndex += xtrAdd1;
    }
  } else {  
    double fp = 0.5*(1+Parameter_MassSplit);
    double fm = 0.5*(1-Parameter_MassSplit);
    for (I0=oneDimSizeL0-1; I0>=0; I0--) {
      for (I1=oneDimSizeL1-1; I1>=0; I1--) {
        for (I2=oneDimSizeL2-1; I2>=0; I2--) {
          for (I3=oneDimSizeL3-1; I3>=0; I3--) {
            //Derivative with respect to Phi_0
            cblas_zdotc_sub(4, &(leftInput[inIndex]), 1, &(rightInput[inIndex]), 1, &h1);
            cblas_zdotc_sub(4, &(leftInput[inIndex]), 1, &(rightInput[inIndex+4]), 1, &h2);
            cblas_zdotc_sub(4, &(leftInput[inIndex+4]), 1, &(rightInput[inIndex]), 1, &h3);
            cblas_zdotc_sub(4, &(leftInput[inIndex+4]), 1, &(rightInput[inIndex+4]), 1, &h4);
	    	    
            calcGamma5VectorProduct(&(leftInput[inIndex]), &(rightInput[inIndex]), g1);
            calcGamma5VectorProduct(&(leftInput[inIndex]), &(rightInput[inIndex+4]), g2);
            calcGamma5VectorProduct(&(leftInput[inIndex+4]), &(rightInput[inIndex]), g3);
            calcGamma5VectorProduct(&(leftInput[inIndex+4]), &(rightInput[inIndex+4]), g4);

            output[outIndex+0] = h1 + Parameter_MassSplit * h4;
            output[outIndex+1] = fp*ComplexI*(g2+g3) + fm*ComplexI*(h2-h3);
            output[outIndex+2] = fp*(g2-g3) + fm*(h2+h3);
            output[outIndex+3] = ComplexI*(g1- Parameter_MassSplit * g4);
  
            outIndex += 4;
            inIndex += 8; 
 	  }
	  inIndex += xtrAdd3;
        }
        inIndex += xtrAdd2;
      } 
      inIndex += xtrAdd1;
    }
  }
}


void FermionMatrixOperations::executeMultiplicationVectorWithDerivativesOfMMdaggerInverse(Complex* omega, Complex* outMMdaggerInverseOmega, double* output, double* phiField, double TOL) {
  int I;

  if (LogLevel>3) printf("Calculating Derivatives of MMdaggerInverse...\n");

  if (Parameter_yN == 0.0) {
    for (I=0; I<vectorLength/2; I++) {
      output[I] = 0;
    }
    return;
  }

  int neededIter = 0;
  solveFermionMatrixLGS(omega, outMMdaggerInverseOmega, phiField, TOL, true, false, -1, neededIter);

  Complex* interim = createFermionVector();
  Complex* interim2 = createFermionVector();
  
  executeFermionDaggerMatrixMultiplication(outMMdaggerInverseOmega, interim, phiField, 0);
  if (Preconditioner_Usage) {
    executeFermionMatrixStaticInverseMultiplication(interim, interim, true, false);
  }
  executeDiracMatrixMultiplication(interim, interim2, false);
  
  Complex fac(-2.0*Parameter_rho,0);
  SSE_ComplexVectorAdditionSPECIAL1(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1,xtraSize2,xtraSize3, fac.x, interim2, interim);
    
  executeMultiplicationVectorWithDerivativesOfB(outMMdaggerInverseOmega, interim, interim2);

  double fac2 = -2 * Parameter_yN;
  for (I=0; I<vectorLength/2; I++) {
    output[I] = fac2 * interim2[I].x;
  }
  
  destroyFermionVector(interim);
  destroyFermionVector(interim2);

  if (LogLevel>3) printf("... ready.\n");
}


void FermionMatrixOperations::exactFermionMatrixConditionNumber(double* phiField, double &eigMin, double &eigMax, double &invCond, bool doubleM, int probeCount, bool quasiHermiteanMode) {
  double TOL = 1E-2;
  
  if (Parameter_yN == 0) {
    eigMin = 0;
    if (quasiHermiteanMode) {
      eigMax = 2*Parameter_rho;    
    } else {
      eigMax = -4*Parameter_rho*Parameter_rho;
    }
    invCond = 0;
    if (doubleM) eigMax *= eigMax;
    return;
  }
  
  Complex* lowEW = calcFermionMatrixARPACKEigenValues(5, probeCount,  phiField, TOL, doubleM, NULL, quasiHermiteanMode, true);
  Complex* highEW = calcFermionMatrixARPACKEigenValues(4, probeCount,  phiField, TOL, doubleM, NULL, quasiHermiteanMode, true);
  
  int I;
  double lowNorm = sqrt(sqr(lowEW[0].x) + sqr(lowEW[0].y));
  double highNorm = sqrt(sqr(highEW[0].x) + sqr(highEW[0].y));
  for (I=0; I<probeCount; I++) {
    double lowN = sqrt(sqr(lowEW[I].x) + sqr(lowEW[I].y));
    double highN = sqrt(sqr(highEW[I].x) + sqr(highEW[I].y));
    
    if (lowN<lowNorm) lowNorm = lowN;
    if (highN>highNorm) highNorm = highN;
  }
  eigMax = highNorm;
  eigMin = lowNorm;
  invCond = eigMin/eigMax;
  delete[] lowEW;
  delete[] highEW;
}


void FermionMatrixOperations::getPreconditionerParameter(bool& use, double& m, double& s) {
  m = Preconditioner_M;
  s = Preconditioner_S;
  use = Preconditioner_Usage;
}


void FermionMatrixOperations::getQPreconditionerParameter(bool& use, double& mu, double& beta) {
  mu = QPreconditioner_mu;
  beta = QPreconditioner_beta;
  use = QPreconditioner_Usage;
}


void FermionMatrixOperations::getRPreconditionerParameter(bool& use, double& m, double& f) {
  m = RPreconditioner_m;
  f = RPreconditioner_f;
  use = RPreconditioner_Usage;
}


void FermionMatrixOperations::printPreconditionerParameter() {
  if (LogLevel>1) {
    if (Preconditioner_Usage) {
      printf("Preconditioning is turned on with PrecM = %1.3f, PrecS = %1.3f.\n", Preconditioner_M, Preconditioner_S);
    } else {
      printf("Preconditioning is turned off.\n");
    }
    if (QPreconditioner_Usage) {
      printf("Q - Preconditioning is turned on with PrecMU = %1.3f, PrecBETA = %1.3f.\n", QPreconditioner_mu, QPreconditioner_beta);
    } else {
      printf("Q - Preconditioning is turned off.\n");
    }
    if (RPreconditioner_Usage) {
      printf("R - Preconditioning is turned on with PrecM = %1.3f, PrecF = %1.3f.\n", RPreconditioner_m, RPreconditioner_f);
    } else {
      printf("R - Preconditioning is turned off.\n");
    }
  }
}


void FermionMatrixOperations::transformToXtraSizeArray(Complex* input, Complex* output) {
  int I3, I2, I1, I0, i;
  int p1, p2;
  int o0 = 8*oneDimSizeL3*oneDimSizeL2*oneDimSizeL1;
  int o1 = 8*oneDimSizeL3*oneDimSizeL2;
  int o2 = 8*oneDimSizeL3;
  int o3 = 8;
  int d0 = 8*(oneDimSizeL3 + xtraSize3)*(oneDimSizeL2 + xtraSize2)*(oneDimSizeL1 + xtraSize1);
  int d1 = 8*(oneDimSizeL3 + xtraSize3)*(oneDimSizeL2 + xtraSize2);
  int d2 = 8*(oneDimSizeL3 + xtraSize3);
  int d3 = 8;
  
  for (I0=oneDimSizeL0-1; I0>=0; I0--) {
    for (I1=oneDimSizeL1-1; I1>=0; I1--) {
      for (I2=oneDimSizeL2-1; I2>=0; I2--) {
        for (I3=oneDimSizeL3-1; I3>=0; I3--) {
	  p1 = I0*o0 + I1*o1 + I2*o2 + I3*o3;
	  p2 = I0*d0 + I1*d1 + I2*d2 + I3*d3;
	  
	  for (i=0; i<8; i++) {
  	    output[p2+i] = input[p1+i];
	  }
        }
      }
    }
  }
}


void FermionMatrixOperations::transformFromXtraSizeArray(Complex* input, Complex* output) {
  int I3, I2, I1, I0, i;
  int p1, p2;
  int o0 = 8*oneDimSizeL3*oneDimSizeL2*oneDimSizeL1 - oneDimSizeL3*8*oneDimSizeL2*oneDimSizeL1;
  int o1 = 8*oneDimSizeL3*oneDimSizeL2 - oneDimSizeL3*8*oneDimSizeL2;
  int o2 = 8*oneDimSizeL3 - oneDimSizeL3*8;
  int o3 = 8;
  int d0 = 8*(oneDimSizeL3 + xtraSize3)*(oneDimSizeL2 + xtraSize2)*(oneDimSizeL1 + xtraSize1) - oneDimSizeL1*8*(oneDimSizeL3 + xtraSize3)*(oneDimSizeL2 + xtraSize2);
  int d1 = 8*(oneDimSizeL3 + xtraSize3)*(oneDimSizeL2 + xtraSize2) - oneDimSizeL2*8*(oneDimSizeL3 + xtraSize3);
  int d2 = 8*(oneDimSizeL3 + xtraSize3) - oneDimSizeL3*8;
  int d3 = 8;
  
  p1 = 0;
  p2 = 0;
  for (I0=0; I0<oneDimSizeL0; I0++) {
    for (I1=0; I1<oneDimSizeL1; I1++) {
      for (I2=0; I2<oneDimSizeL2; I2++) {
        for (I3=0; I3<oneDimSizeL3; I3++) {
	  for (i=0; i<8; i++) {
  	    output[p1+i] = input[p2+i];
	  }
          p1 += o3;
          p2 += d3;
        }
        p1 += o2;
        p2 += d2;
      }
      p1 += o1;
      p2 += d1;
    }
    p1 += o0;
    p2 += d0;
  }
}


void FermionMatrixOperations::setxFFTusage(bool xFFTu) {
  xFFTusage = xFFTu;
  if (xFFTu) {
    if (LogLevel>1) printf("Usage of xFFT activated!!!\n");
  } else {
    if (LogLevel>1) printf("Usage of xFFT deactivated!!!\n");
  }
}


bool FermionMatrixOperations::getxFFTusage() {
  return xFFTusage;
}


void FermionMatrixOperations::setxFFT_DistributedFFT_ThreadCount(int xFFTThreadCnt) {
  if (xFFTThreadCnt != xFFT_DistributedFFT_ThreadCount) {
    if (LogLevel>1) printf("Set number of threads for distributed xFFT to %d.\n", xFFTThreadCnt);  
  }
  xFFT_DistributedFFT_ThreadCount = xFFTThreadCnt;
}


int FermionMatrixOperations::getxFFT_DistributedFFT_ThreadCount() {
  return xFFT_DistributedFFT_ThreadCount;
}


void FermionMatrixOperations::testFourierTrafo(bool forw) {
  if (!xFFTusage) {
    if (LogLevel>0) printf("Used FFT is: FFTW\n");
    return;
  }
  if (LogLevel>0) printf("Used FFT is: xFFT\n");
    
  if (LogLevel>0) printf("\nTesting xFFT-Fast-Fourier-Transformation with forward=%d...\n", forw);
  int VL = vectorLength;
  int I;
  for (I=0; I<VL; I++) {
    Inverse_b[I] = Complex(2*(zufall()-0.5),2*(zufall()-0.5));
  }
  transformToXtraSizeArray(Inverse_b, Inverse_b);

  setxFFTusage(false);
  performFFT(Inverse_b,Inverse_p, forw);
  setxFFTusage(true);
  performFFT(Inverse_b,Inverse_s, forw);
  
  transformFromXtraSizeArray(Inverse_p,Inverse_p);
  transformFromXtraSizeArray(Inverse_s,Inverse_s);
  
  double diff = 0;
  int correct = 0;
  for (I=0; I<VL; I++) {
    if ((I>=0) && (I<8)) {
      if (LogLevel>0) { 
        printf("xFFT: ");
	Inverse_s[I].print();
        printf("FFTW: ");
        Inverse_p[I].print();
      }
    }
    double dsqr = sqr(Inverse_p[I].x-Inverse_s[I].x)+sqr(Inverse_p[I].y-Inverse_s[I].y);
    if (sqrt(dsqr)<1E-14*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3) correct++;
    diff += dsqr;
  }
  diff = sqrt(diff);
  if (LogLevel>0) printf("Difference %1.15f:\n",diff);
  if (LogLevel>0) printf("Number of correct entries %d:\n",correct);
  if ((diff>1E-14*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3) || (isNaN(diff)) || (correct != 8*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3)) {
    printf("Error in xFFT!!!   ==> EXITING\n");
    exit(0);
  }
}


void FermionMatrixOperations::testDistributedFourierTrafo(bool forw) {
  if (!xFFTusage) {
    if (LogLevel>0) printf("Used FFT is: FFTW\n");
    return;
  }
  if (LogLevel>0) printf("Used FFT is: xFFT\n");
    
  if (LogLevel>0) printf("\nTesting Distributed xFFT-Fast-Fourier-Transformation with forward=%d...\n", forw);
  int VL = vectorLength;
  int I;
  for (I=0; I<VL; I++) {
    Inverse_b[I] = Complex(2*(zufall()-0.5),2*(zufall()-0.5));
  }
  transformToXtraSizeArray(Inverse_b, Inverse_b);
  copyFermionVectorToDistributedFermionVector(Inverse_b, DistributedInverse_b);
  
  setxFFTusage(false);
  performDistributedFFT(DistributedInverse_b,DistributedInverse_p, forw);  
  setxFFTusage(true);
  performDistributedFFT(DistributedInverse_b,DistributedInverse_s, forw);  

  copyDistributedFermionVectorToFermionVector(DistributedInverse_p, Inverse_p);
  copyDistributedFermionVectorToFermionVector(DistributedInverse_s, Inverse_s);
  transformFromXtraSizeArray(Inverse_p,Inverse_p);
  transformFromXtraSizeArray(Inverse_s,Inverse_s);
  
  double diff = 0;
  int correct = 0;
  for (I=0; I<VL; I++) {
    if ((I>=0) && (I<8)) {
      if (LogLevel>0) { 
        printf("xFFT: ");
	Inverse_s[I].print();
        printf("FFTW: ");
        Inverse_p[I].print();
      }
    }
    double dsqr = sqr(Inverse_p[I].x-Inverse_s[I].x)+sqr(Inverse_p[I].y-Inverse_s[I].y);
    if (sqrt(dsqr)<1E-14*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3) correct++;
    diff += dsqr;
  }
  diff = sqrt(diff);
  if (LogLevel>0) printf("Difference %1.15f:\n",diff);
  if (LogLevel>0) printf("Number of correct entries %d:\n",correct);
  if ((diff>1E-14*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3) || (isNaN(diff)) || (correct != 8*oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3)) {
    printf("Error in xFFT!!!   ==> EXITING\n");
    exit(0);
  }
}


void FermionMatrixOperations::executeStarTransformationT(Complex* input, Complex* output, bool adjoining, bool inverse) {
  Complex dummy[8];
  double facX = 1.0;
  double facY = 1.0;
  if (adjoining) {
    facY = -1.0;  
  }
  if (inverse) {
    facX *= -1.0;
    facY *= -1.0;    
  }

  int I0,I1,I2,I3, i;
  int countXtr = 0;
  int xtrAdd1 = xtraSize1*8*(oneDimSizeL3+xtraSize3)*(oneDimSizeL2+xtraSize2);
  int xtrAdd2 = xtraSize2*8*(oneDimSizeL3+xtraSize3);
  int xtrAdd3 = xtraSize3*8;
  for (I0=0; I0<oneDimSizeL0; I0++) {
    for (I1=0; I1<oneDimSizeL1; I1++) {
      for (I2=0; I2<oneDimSizeL2; I2++) {
        for (I3=0; I3<oneDimSizeL3; I3++) {
          for (i=0; i<8; i++) {
            dummy[i].x = facX*input[countXtr+i].x;
            dummy[i].y = facY*input[countXtr+i].y;
	  }
	  
	  if (adjoining && inverse) {
  	    output[countXtr+0].x = -dummy[5].y;   // +i
	    output[countXtr+0].y =  dummy[5].x;   // +i
	    output[countXtr+1].x =  dummy[4].y;   // -i
	    output[countXtr+1].y = -dummy[4].x;   // -i
	    output[countXtr+2].x = -dummy[7].y; 
	    output[countXtr+2].y =  dummy[7].x; 
	    output[countXtr+3].x =  dummy[6].y; 
	    output[countXtr+3].y = -dummy[6].x; 
	    output[countXtr+4].x =  dummy[1].y; 
	    output[countXtr+4].y = -dummy[1].x; 
	    output[countXtr+5].x = -dummy[0].y; 
	    output[countXtr+5].y =  dummy[0].x; 
	    output[countXtr+6].x =  dummy[3].y; 
	    output[countXtr+6].y = -dummy[3].x; 
	    output[countXtr+7].x = -dummy[2].y; 
	    output[countXtr+7].y =  dummy[2].x; 
	  } else {
  	    output[countXtr+0].x =  dummy[5].y;   // -i
	    output[countXtr+0].y = -dummy[5].x;   // -i
	    output[countXtr+1].x = -dummy[4].y;   // +i
	    output[countXtr+1].y =  dummy[4].x;   // +i
	    output[countXtr+2].x =  dummy[7].y; 
	    output[countXtr+2].y = -dummy[7].x; 
	    output[countXtr+3].x = -dummy[6].y; 
	    output[countXtr+3].y =  dummy[6].x; 
	    output[countXtr+4].x = -dummy[1].y; 
	    output[countXtr+4].y =  dummy[1].x; 
	    output[countXtr+5].x =  dummy[0].y; 
	    output[countXtr+5].y = -dummy[0].x; 
	    output[countXtr+6].x = -dummy[3].y; 
	    output[countXtr+6].y =  dummy[3].x; 
	    output[countXtr+7].x =  dummy[2].y; 
	    output[countXtr+7].y = -dummy[2].x; 
	  }
	  
          countXtr+=8;
	}
        countXtr += xtrAdd3;
      }
      countXtr += xtrAdd2;
    }
    countXtr += xtrAdd1;
  }
}


void FermionMatrixOperations::executePiModeRemoverOperator(Complex* input, Complex* output, bool inFourierSpace) {
  Complex* source;
  Complex* target;
  if (inFourierSpace) {
    source = input;
    target = output;  
    if (source != target) {
      SSE_ZCopy(vectorLengthXtrSize, source, 1, target, 1);
    }
  } else {
    performFFT(input, InputVectorFourierTransform, ExtremeFFT4D_Forward); 
    source = InputVectorFourierTransform;
    target = InputVectorFourierTransform;  
    double fac = 1.0/(oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3);
    for (int I=0; I<vectorLengthXtrSize; I++) {
      target[I].x *= fac;
      target[I].y *= fac;
    }
  }
 
  for (int I=0; I<15; I++) {
    int p = (int) Index_PiModes[I];
    if (p>=0) {
      for (int I2=0; I2<8; I2++) {
        target[p+I2].x = 0.0;
        target[p+I2].y = 0.0;
      }
    }
  }

  if (!inFourierSpace) {
    performFFT(target, output, ExtremeFFT4D_Backward);
  }
}


void FermionMatrixOperations::executeDistributedPiModeRemoverOperator(DistributedMemoryObject* input, DistributedMemoryObject* output, bool inFourierSpace) {
  if (threadedOps==NULL) {
    printf("ERROR in FermionMatrixOperations::executeDistributedPiModeRemoverOperator: threadedOps==NULL\n");
    exit(0);
  }
  if (inFourierSpace) {
    threadedOps->perform_PiModeRemoverOperatorMultiplicationInFourierSpace(input, output, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, Distributed_Index_PiModes); 
  } else {
    performDistributedFFT(input, DistributedInputVectorFourierTransform, ExtremeFFT4D_Forward); 
    double fac = 1.0/(oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3);
    Complex alpha(fac, 0);
    threadedOps->multiplyFermionVectorWithComplexNumber(DistributedInputVectorFourierTransform, DistributedInputVectorFourierTransform, alpha, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
    threadedOps->perform_PiModeRemoverOperatorMultiplicationInFourierSpace(DistributedInputVectorFourierTransform, DistributedInputVectorFourierTransform, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, Distributed_Index_PiModes);     
    performDistributedFFT(DistributedInputVectorFourierTransform, output, ExtremeFFT4D_Backward); 
  } 
}


void FermionMatrixOperations::executeTau2(Complex* input, Complex* output) {
  Complex dummy[8];

  int I0,I1,I2,I3, i;
  int countXtr = 0;
  int xtrAdd1 = xtraSize1*8*(oneDimSizeL3+xtraSize3)*(oneDimSizeL2+xtraSize2);
  int xtrAdd2 = xtraSize2*8*(oneDimSizeL3+xtraSize3);
  int xtrAdd3 = xtraSize3*8;
  for (I0=0; I0<oneDimSizeL0; I0++) {
    for (I1=0; I1<oneDimSizeL1; I1++) {
      for (I2=0; I2<oneDimSizeL2; I2++) {
        for (I3=0; I3<oneDimSizeL3; I3++) {
          for (i=0; i<8; i++) {
            dummy[i].x = input[countXtr+i].x;
            dummy[i].y = input[countXtr+i].y;
	  }
	  
	  output[countXtr+0].x =  dummy[4].y;   
	  output[countXtr+0].y = -dummy[4].x;   
	  output[countXtr+1].x =  dummy[5].y;   
	  output[countXtr+1].y = -dummy[5].x;   
	  output[countXtr+2].x =  dummy[6].y; 
	  output[countXtr+2].y = -dummy[6].x; 
	  output[countXtr+3].x =  dummy[7].y; 
	  output[countXtr+3].y = -dummy[7].x; 
	  output[countXtr+4].x = -dummy[0].y; 
	  output[countXtr+4].y =  dummy[0].x; 
	  output[countXtr+5].x = -dummy[1].y; 
	  output[countXtr+5].y =  dummy[1].x; 
	  output[countXtr+6].x = -dummy[2].y; 
	  output[countXtr+6].y =  dummy[2].x; 
	  output[countXtr+7].x = -dummy[3].y; 
	  output[countXtr+7].y =  dummy[3].x; 
	  
          countXtr+=8;
	}
        countXtr += xtrAdd3;
      }
      countXtr += xtrAdd2;
    }
    countXtr += xtrAdd1;
  }
}


void FermionMatrixOperations::executeGamma0(Complex* input, Complex* output) {
  Complex dummy[8];

  int I0,I1,I2,I3, i;
  int countXtr = 0;
  int xtrAdd1 = xtraSize1*8*(oneDimSizeL3+xtraSize3)*(oneDimSizeL2+xtraSize2);
  int xtrAdd2 = xtraSize2*8*(oneDimSizeL3+xtraSize3);
  int xtrAdd3 = xtraSize3*8;
  for (I0=0; I0<oneDimSizeL0; I0++) {
    for (I1=0; I1<oneDimSizeL1; I1++) {
      for (I2=0; I2<oneDimSizeL2; I2++) {
        for (I3=0; I3<oneDimSizeL3; I3++) {
          for (i=0; i<8; i++) {
            dummy[i].x = input[countXtr+i].x;
            dummy[i].y = input[countXtr+i].y;
	  }
	  
	  output[countXtr+0].x =  dummy[2].x;   
	  output[countXtr+0].y =  dummy[2].y;   
	  output[countXtr+1].x =  dummy[3].x;   
	  output[countXtr+1].y =  dummy[3].y;   
	  output[countXtr+2].x =  dummy[0].x; 
	  output[countXtr+2].y =  dummy[0].y; 
	  output[countXtr+3].x =  dummy[1].x; 
	  output[countXtr+3].y =  dummy[1].y; 
	  output[countXtr+4].x =  dummy[6].x; 
	  output[countXtr+4].y =  dummy[6].y; 
	  output[countXtr+5].x =  dummy[7].x; 
	  output[countXtr+5].y =  dummy[7].y; 
	  output[countXtr+6].x =  dummy[4].x; 
	  output[countXtr+6].y =  dummy[4].y; 
	  output[countXtr+7].x =  dummy[5].x; 
	  output[countXtr+7].y =  dummy[5].y; 
	  
          countXtr+=8;
	}
        countXtr += xtrAdd3;
      }
      countXtr += xtrAdd2;
    }
    countXtr += xtrAdd1;
  }
}


void FermionMatrixOperations::executeGamma5(Complex* input, Complex* output) {
  Complex dummy[8];

  int I0,I1,I2,I3, i;
  int countXtr = 0;
  int xtrAdd1 = xtraSize1*8*(oneDimSizeL3+xtraSize3)*(oneDimSizeL2+xtraSize2);
  int xtrAdd2 = xtraSize2*8*(oneDimSizeL3+xtraSize3);
  int xtrAdd3 = xtraSize3*8;
  for (I0=0; I0<oneDimSizeL0; I0++) {
    for (I1=0; I1<oneDimSizeL1; I1++) {
      for (I2=0; I2<oneDimSizeL2; I2++) {
        for (I3=0; I3<oneDimSizeL3; I3++) {
          for (i=0; i<8; i++) {
            dummy[i].x = input[countXtr+i].x;
            dummy[i].y = input[countXtr+i].y;
	  }
	  
	  output[countXtr+0].x =  dummy[0].x;   
	  output[countXtr+0].y =  dummy[0].y;   
	  output[countXtr+1].x =  dummy[1].x;   
	  output[countXtr+1].y =  dummy[1].y;   
	  output[countXtr+2].x = -dummy[2].x; 
	  output[countXtr+2].y = -dummy[2].y; 
	  output[countXtr+3].x = -dummy[3].x; 
	  output[countXtr+3].y = -dummy[3].y; 
	  output[countXtr+4].x =  dummy[4].x; 
	  output[countXtr+4].y =  dummy[4].y; 
	  output[countXtr+5].x =  dummy[5].x; 
	  output[countXtr+5].y =  dummy[5].y; 
	  output[countXtr+6].x = -dummy[6].x; 
	  output[countXtr+6].y = -dummy[6].y; 
	  output[countXtr+7].x = -dummy[7].x; 
	  output[countXtr+7].y = -dummy[7].y; 
	  
          countXtr+=8;
	}
        countXtr += xtrAdd3;
      }
      countXtr += xtrAdd2;
    }
    countXtr += xtrAdd1;
  }
}


/*
c     %-----------------------------------------------%
c     |                                               | 
c     | Specifications for ARPACK usage are set       | 
c     | below:                                        |
c     |                                               |
c     |    1) NEV = 4  asks for 4 eigenvalues to be   |  
c     |       computed.                               | 
c     |                                               |
c     |    2) NCV = 20 sets the length of the Arnoldi |
c     |       factorization                           |
c     |                                               |
c     |    3) This is a standard problem              |
c     |         (indicated by bmat  = 'I')            |
c     |                                               |
c     |    4) Ask for the NEV eigenvalues of          |
c     |       largest magnitude                       |
c     |         (indicated by which = 'LM')           |
c     |       See documentation in ZNAUPD  for the     |
c     |       other options SM, LR, SR, LI, SI.       | 
c     |                                               |
c     | Note: NEV and NCV must satisfy the following  |
c     | conditions:                                   |
c     |              NEV <= MAXNEV                    |
c     |          NEV + 2 <= NCV <= MAXNCV             |
c     |                                               |
c     |                                                     |
c     | Specification of stopping rules and initial         |
c     | conditions before calling ZNAUPD                     |
c     |                                                     |
c     | TOL  determines the stopping criterion.             |
c     |                                                     |
c     |      Expect                                         |
c     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
c     |               computed   true                       |
c     |                                                     |
c     |      If TOL .le. 0,  then TOL <- macheps            |
c     |           (machine precision) is used.              |
c     |                                                     |
c     | IDO  is the REVERSE COMMUNICATION parameter         |
c     |      used to specify actions to be taken on return  |
c     |      from ZNAUPD . (see usage below)                 |
c     |                                                     |
c     |      It MUST initially be set to 0 before the first |
c     |      call to ZNAUPD .                                | 
c     |                                                     |
c     | INFO on entry specifies starting vector information |
c     |      and on return indicates error codes            |
c     |                                                     |
c     |      Initially, setting INFO=0 indicates that a     | 
c     |      random starting vector is requested to         |
c     |      start the ARNOLDI iteration.  Setting INFO to  |
c     |      a nonzero value on the initial call is used    |
c     |      if you want to specify your own starting       |
c     |      vector (This vector must be placed in RESID).  | 
c     |                                                     |
c     | The work array WORKL is used in ZNAUPD  as           | 
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.                                  |
c     |                                                     |
c     | Specification of Algorithm Mode:                  |
c     |                                                   |
c     | This program uses the exact shift strategy        |
c     | (indicated by setting IPARAM(1) = 1).             |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of ZNAUPD  is used     |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | ZNAUPD .                                           |
c        | Eigenvalues are returned in the one           |
c        | dimensional array D and the corresponding     |
c        | eigenvectors are returned in the first        |
c        | NCONV (=IPARAM(5)) columns of the two         |
c        | dimensional array V if requested.  Otherwise, |
c        | an orthogonal basis for the invariant         |
c        | subspace corresponding to the eigenvalues in  |
c        | D is returned in V.                           |
c     
c     what:     0: largest real part
c               1: smallest real value
c               2: largest imaginary part
c               3: smallest imaginary value
c               4: largest absolute value
c               5: smallest abolute value
c
c     %-----------------------------------------------------%
*/
Complex* FermionMatrixOperations::calcFermionMatrixARPACKEigenValues(int what, int nev, double* phiField, double TOL, bool doubleFermionMatrix, Complex** rightEigenVectors, bool quasiHermiteanMode, bool inFourierSpace) {
#ifdef useARPACK
  int n = vectorLength;
  int ncv = 2*nev+2;
  if (ncv<8) ncv = 8;
  int ldv=n;
  
  int iparam[11];
  int ipntr[14];
    
  //Allocate Memory
  int* select = new int[ncv];   
  Complex* d = new Complex[ncv];
  Complex* v = new Complex[ldv*ncv];
  Complex* workd = new Complex[3*n];
  Complex* workev = new Complex[2*ncv];
  Complex* resid = new Complex[n];
  Complex* workl = new Complex[3*ncv*ncv+5*ncv];
  double* rwork = new double[ncv];

  char bmat[1];
  char which[2];
  int  ido, lworkl, info, ierr, ishfts, maxitr, mode1;
       
  Complex sigma;
  int rvec;
  
  bmat[0]  = 'I';
  
  if (what==0) {
    which[0] = 'L';
    which[1] = 'R';
  } else if (what==1) {
    which[0] = 'S';
    which[1] = 'R';
  } else if (what==2) {
    which[0] = 'L';
    which[1] = 'I';
  } else if (what==3) {
    which[0] = 'S';
    which[1] = 'I';
  } else if (what==4) {
    which[0] = 'L';
    which[1] = 'M';
  } else if (what==5) {
    which[0] = 'S';
    which[1] = 'M';
  } else {
    printf("calcFermionMatrixARPACKEigenValues: What-Mode not supported!!!\n");
    exit(0);
  }
     
  lworkl  = 3*ncv*ncv+5*ncv; 
  double tol = TOL;
  ido    = 0;
  info   = 0;

  ishfts = 1;
  maxitr = 10000;           //Maximale Anzahl erlaubter Iterationen !!!
  mode1 = 1;

  iparam[0] = ishfts;                
  iparam[2] = maxitr;   
  iparam[3] = 1;          //Blocksize. works only for = 1         
  iparam[6] = mode1;


  while(true) {  
    znaupd_(&ido, &(bmat[0]), &n, &(which[0]), &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, 
            workd, workl, &lworkl, rwork, &info);

    int xInd = ipntr[0]-1;
    int yInd = ipntr[1]-1;
    int zInd = ipntr[2]-1;

    if (ido == 1) {
      for (int I=0; I<n; I++) {
        workd[I+zInd].x=workd[I+xInd].x;
        workd[I+zInd].y=workd[I+xInd].y;
      }
    }       
    
    if ((ido==1) || (ido==-1)) {
      //Matrix application: Operand in xInd, result in yInd
      for (int I=0; I<n; I++) {
        workd[I+yInd]= Complex(1.3,-1.2)*workd[I+xInd];
      } 
      Complex* ARinput = ((Complex*) &(workd[xInd].x));
      Complex* ARoutput = ((Complex*) &(workd[yInd].x));
      transformToXtraSizeArray(ARinput,Inverse_p);
            
      if (doubleFermionMatrix) {
        if (quasiHermiteanMode) {
          executeFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(Inverse_p, Inverse_s, phiField, true, inFourierSpace);
	} else {
          executeFermionMatrixFermionDaggerMatrixMultiplication(Inverse_p, Inverse_s, phiField, 2, 1, inFourierSpace);	
	} 	
      } else {
        if (quasiHermiteanMode) {
          executeFermionQuasiHermiteanMatrixMultiplication(Inverse_p, Inverse_s, phiField, true, false, inFourierSpace);
	} else {
          executeFermionMatrixMultiplication(Inverse_p, Inverse_s, phiField, false, NULL, NULL, 1, 1);
	} 	
      }

      transformFromXtraSizeArray(Inverse_s, ARoutput);
    }
   
    if ((ido == 3) || (ido == 4)) {
      printf("calcFermionMatrixARPACKEigenValues: Not supplied mode!!!\n");
      exit(0);
    }
    
    if (ido == 99) {
      break;
    }    
  }

  if (info<0) {
    printf("calcFermionMatrixARPACKEigenValues: Error info code: %d\n",info);
    exit(0);
  } 
  
  if (info==1) {
    printf("calcFermionMatrixARPACKEigenValues: Maximum number of iterations (%d) reached. Number can be changed in routine.\n", maxitr);
    exit(0);
  } 
  if (info==3) {
    printf("calcFermionMatrixARPACKEigenValues: No shifts could be applied during implicit Arnoldi update, try increasing NCV or derease Accuracy.\n");
    exit(0);
  }  
  
  rvec = false;
  char howMany = 'P';
  if (rightEigenVectors != NULL) {
    rvec = true;  
    howMany = 'A';
  }
  

  zneupd_(&rvec, &howMany, select, d, v, &ldv, &sigma,
          workev, &(bmat[0]), &n, &(which[0]), &nev, &tol, resid, &ncv,
          v, &ldv, iparam, ipntr, workd, workl, &lworkl,
          rwork, &ierr);

  if (ierr!=0) {
    printf("calcFermionMatrixARPACKEigenValues: Error in zneupd: %d\n",ierr);
    if (ierr==-14) {
      printf("ZNAUPD did not find any eigenvalues to sufficient accuracy.\n");
    }
    exit(0);
  }

  //Get Result
  Complex* eigenV = new Complex[nev];
  for (int I=0; I<nev; I++) {
   eigenV[I] = d[I];
  }
  
  //Get Right-Eigenvectors
  if (rightEigenVectors != NULL) {
    for (int I=0; I<nev; I++) {
      Complex* vec = rightEigenVectors[I];
      int vpos = I*vectorLength;
      for (int I2=0; I2<vectorLength; I2++) {
        vec[I2].x = v[vpos+I2].x;
        vec[I2].y = v[vpos+I2].y;     
      }
      transformToXtraSizeArray(vec,vec);
    }
  }

  //Free Memory
  delete[] select;   
  delete[] d;
  delete[] v;
  delete[] workd;
  delete[] workev;
  delete[] resid;
  delete[] workl;
  delete[] rwork;

  return eigenV;

#else
  printf("ARPACK is not included!!!\n");
  exit(0);
#endif
}
