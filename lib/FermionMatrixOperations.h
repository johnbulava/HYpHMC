#ifndef FermionMatrixOperations_included
#define FermionMatrixOperations_included

#include "Global.h"
#ifdef useBLAS
  #include CBLASIncludeFile
#endif
#include FFTWIncludeFile
#include "Complex.h"
#include "Quat.h"
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include "NeubergerMatrix.h"
#include "WilsonMatrix.h"
#include "SSEroutines.h"
#include "ExtremeFFT4D.h"
#include "GeneralChebyshevApproximation.h"
#include "MultiThreadedOperations.h"


#define FFTPlanDBMax 1000


class FermionMatrixOperations {
protected:
  void ini();
  void desini();

  Complex* InputVectorFourierTransform;  
  Complex* InterimVectorFourierTransform;    
  Complex* Inverse_rest;
  Complex* Inverse_b;
  Complex* Inverse_p;
  Complex* Inverse_s;
  DistributedMemoryObject* DistributedInputVectorFourierTransform;  
  DistributedMemoryObject* DistributedInterimVectorFourierTransform;    
  DistributedMemoryObject* DistributedInverse_rest;
  DistributedMemoryObject* DistributedInverse_b;
  DistributedMemoryObject* DistributedInverse_p;
  DistributedMemoryObject* DistributedInverse_s;
  fftw_plan FFTPlanDataBase[FFTPlanDBMax];
  Complex* InputVectorDataBase[FFTPlanDBMax];
  Complex* OutputVectorDataBase[FFTPlanDBMax];
  bool BackwardForwardDataBase[FFTPlanDBMax];
  int FFTDataBaseCounter;
  
  ExtremeFFT4D* xFFT;
  MultiThreadedOperations* threadedOps;

  double Parameter_rho;
  double Parameter_r;
  double Parameter_yN;
  double Parameter_MassSplit;
  double Parameter_ExplicitMass;
  bool Parameter_AntiPeriodicBoundaryConditionInTime;
  NeubergerMatrix* NeubergerDiracOp;
  WilsonMatrix* WilsonDiracOp;
  int oneDimSizeL0;
  int oneDimSizeL1;
  int oneDimSizeL2;
  int oneDimSizeL3;  
  int oneDimSizeLargest;
  int timeDirection;
  int vectorLength;
  int vectorLengthXtrSize;
  
  void generateOpApplicationData();
  void generateDistributedOpApplicationData();
  void setSizeDistributed(int sizeL0, int sizeL1, int sizeL2, int sizeL3);  

  double Preconditioner_M;
  double Preconditioner_S;
  bool Preconditioner_Usage;
  double QPreconditioner_mu;
  double QPreconditioner_beta;
  bool QPreconditioner_Usage;
  double RPreconditioner_m;
  double RPreconditioner_f;
  bool RPreconditioner_Usage;  
  bool xFFTusage;
  int xFFT_DistributedFFT_ThreadCount;
  
public:
  Complex* NeubergerDiracOpApplicationEWData;
  Complex* NeubergerDiracOpApplicationEWDataRescaled;
  Complex* NeubergerDiracPlusExplicitMassOpApplicationEWData;
  Complex* NeubergerDiracPlusExplicitMassOpApplicationEWDataRescaled;
  Complex* NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWData;
  Complex* NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWDataRescaled;
  Complex* NeubergerDiracOpApplicationSinPStdData;
  Complex* NeubergerDiracOpApplicationSinPDagData;  
  Complex* NeubergerWithChiStaticFermionInverseOpApplicationMainDiagonalData;
  Complex* NeubergerWithChiStaticFermionInverseOpApplicationSubDiagonalData;
  Complex* NeubergerWithChiStaticMdagMFermionInverseOpApplicationMainDiagonalData;
  Complex* NeubergerWithChiStaticMdagMFermionInverseOpApplicationSubDiagonalData;
  Complex* NeubergerWithChiStaticMdagMFermionInverseOpApplicationSSEAuxData;
  Complex* QPreconditionerDiagonalData;
  Complex* QPreconditionerDiagonalDataNormalized;
  Complex* inverseQPreconditionerDiagonalData;
  Complex* inverseQPreconditionerDiagonalDataNormalized;
  Complex* RPreconditionerDiagonalData;
  Complex* RPreconditionerDiagonalDataNormalized;
  Complex* inverseRPreconditionerDiagonalData;
  Complex* inverseRPreconditionerDiagonalDataNormalized;  
  Complex* RSQRPreconditionerDiagonalData;
  Complex* RSQRPreconditionerDiagonalDataNormalized;
  Complex* inverseRSQRPreconditionerDiagonalData;
  Complex* inverseRSQRPreconditionerDiagonalDataNormalized;
  Complex* halfMomentumForwardFFTfactors;
  Complex* halfMomentumBackwardFFTfactors;
  

  long int* Index_PlusPiPiPiPi;
  long int* Index_PlusPiPiPiPiXtrSize;
  long int* Index_PiModes;

  DistributedMemoryObject* Distributed_NeubergerDiracOpApplicationEWData;
  DistributedMemoryObject* Distributed_NeubergerDiracOpApplicationEWDataRescaled;
  DistributedMemoryObject* Distributed_NeubergerDiracPlusExplicitMassOpApplicationEWData;
  DistributedMemoryObject* Distributed_NeubergerDiracPlusExplicitMassOpApplicationEWDataRescaled;
  DistributedMemoryObject* Distributed_NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWData;
  DistributedMemoryObject* Distributed_NeubergerDiracUnityMinusDiracQuasiInverseOpApplicationEWDataRescaled;
  DistributedMemoryObject* Distributed_NeubergerDiracOpApplicationSinPStdData;
  DistributedMemoryObject* Distributed_NeubergerDiracOpApplicationSinPDagData;  
  DistributedMemoryObject* Distributed_NeubergerWithChiStaticFermionInverseOpApplicationMainDiagonalData;
  DistributedMemoryObject* Distributed_NeubergerWithChiStaticFermionInverseOpApplicationSubDiagonalData;
  DistributedMemoryObject* Distributed_NeubergerWithChiStaticMdagMFermionInverseOpApplicationMainDiagonalData;
  DistributedMemoryObject* Distributed_NeubergerWithChiStaticMdagMFermionInverseOpApplicationSubDiagonalData;
  DistributedMemoryObject* Distributed_NeubergerWithChiStaticMdagMFermionInverseOpApplicationSSEAuxData;
  DistributedMemoryObject* Distributed_QPreconditionerDiagonalData;
  DistributedMemoryObject* Distributed_QPreconditionerDiagonalDataNormalized;
  DistributedMemoryObject* Distributed_inverseQPreconditionerDiagonalData;
  DistributedMemoryObject* Distributed_inverseQPreconditionerDiagonalDataNormalized;
  DistributedMemoryObject* Distributed_RPreconditionerDiagonalData;
  DistributedMemoryObject* Distributed_RPreconditionerDiagonalDataNormalized;
  DistributedMemoryObject* Distributed_inverseRPreconditionerDiagonalData;
  DistributedMemoryObject* Distributed_inverseRPreconditionerDiagonalDataNormalized;  
  DistributedMemoryObject* Distributed_RSQRPreconditionerDiagonalData;
  DistributedMemoryObject* Distributed_RSQRPreconditionerDiagonalDataNormalized;
  DistributedMemoryObject* Distributed_inverseRSQRPreconditionerDiagonalData;
  DistributedMemoryObject* Distributed_inverseRSQRPreconditionerDiagonalDataNormalized;  
  DistributedMemoryObject* Distributed_halfMomentumForwardFFTfactors;
  DistributedMemoryObject* Distributed_halfMomentumBackwardFFTfactors;

  DistributedMemoryObject* Distributed_Index_PlusPiPiPiPi;
  DistributedMemoryObject* Distributed_Index_PlusPiPiPiPiXtrSize;
  DistributedMemoryObject* Distributed_Index_PiModes;


  FermionMatrixOperations();
  FermionMatrixOperations(int sizeL0, int sizeL1, int sizeL2, int sizeL3, double rho, double r, double yN);
  ~FermionMatrixOperations();

  void setSize(int sizeL0, int sizeL1, int sizeL2, int sizeL3);
  void setDiracParameters(double rho, double r);
  void getDiracParameters(double& rho, double& r);
  void setYukawaCoupling(double yN);
  double getYukawaCoupling();
  void setPreconditioner(bool usePrecon, double m, double s);
  void setQPreconditioner(bool useqPrecon, double mu, double beta);
  void setRPreconditioner(bool userPrecon, double m, double f);
  void setMassSplitRatio(double split);
  void setExplicitMass(double mass);
  void setAntiPeriodicBoundaryConditionsInTime(bool aPerBCinTime);
  double getMassSplitRatio();
  double getExplicitMass();
  void activateMultiThreadedOps(int OpMode, bool enforceCoreTie);  
  void deactivateMultiThreadedOps();  
  bool isMultiThreadedOpsActivated();  
  MultiThreadedOperations* getMultiThreadedOps();
  
  
  Complex* createFermionVector();
  Complex* createFermionVector(int localSize); 
  void destroyFermionVector(Complex* &v);
  DistributedMemoryObject* createDistributedFermionVector();
  DistributedMemoryObject* createDistributedUniformLatticeComplexVector(int localSize); 
  void destroyDistributedComplexVector(DistributedMemoryObject* &v);
  void copyFermionVectorToDistributedFermionVector(Complex* vec, DistributedMemoryObject* memObj);
  void copyDistributedFermionVectorToFermionVector(DistributedMemoryObject* memObj, Complex* vec);
  void copyPhiFieldUniformlyToDistributedPhiField(double* phiField, DistributedMemoryObject* memObj);
  void zeroFermionVector(Complex* v);
  void zeroDistributedFermionVector(DistributedMemoryObject* v);
  
  void planeWaveFermionVector(Complex* v, double p0, double p1, double p2, double p3);
  void fillGaussRandomVector(Complex* v, int count);

  void executeFermionMatrixMultiplication(Complex* input, Complex* output, double* phiField, bool ScalarsForInvSolver, Complex* s1, Complex* s2, int preconLevel, int QPrec);
  void executeFermionDaggerMatrixMultiplication(Complex* input, Complex* output, double* phiField, int QPrec);

  void executeFermionQuasiHermiteanMatrixMultiplication(Complex* input, Complex* output, double* phiField, bool useRPrec, bool daggered, bool inFourierSpace);

  void multiplyFermionVectorWithTimeIndexedComplexScalars(Complex* input, Complex* output, Complex* scalars);

  void executeFermionMatrixFermionDaggerMatrixMultiplication(Complex* input, Complex* output, double* phiField, int preconLevel, int QPrec, bool inFourierSpace);
  void executeFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(Complex* input, Complex* output, double* phiField, bool useRPrec, bool inFourierSpace);
  void executeDistributedFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(DistributedMemoryObject* input, DistributedMemoryObject* output, DistributedMemoryObject* phiField, bool useRPrec, bool inFourierSpace);


  void executeDiracMatrixMultiplication(Complex* input, Complex* output, bool inFourierSpace);
  void executeDiracDaggerMatrixMultiplication(Complex* input, Complex* output, bool inFourierSpace);
  void executeDiracDaggerPlusExplicitMassMatrixMultiplication(Complex* input, Complex* output, bool inFourierSpace);
  void executeDiracUnityMinusDiracQuasiInverseMatrixMultiplication(Complex* input, Complex* output, bool inFourierSpace);
  void executeDistributedDiracUnityMinusDiracQuasiInverseMatrixMultiplication(DistributedMemoryObject* input, DistributedMemoryObject* output, bool inFourierSpace);
  void executeDiracUnityMinusDiracQuasiInverseDaggerMatrixMultiplication(Complex* input, Complex* output, bool inFourierSpace);
  void executeDistributedDiracUnityMinusDiracQuasiInverseDaggerMatrixMultiplication(DistributedMemoryObject* input, DistributedMemoryObject* output, bool inFourierSpace);
  void executeDistributedMultiCombinedOperator_VA1_RorCopy_VAc_RorCopy_D_MatrixMultiplication(DistributedMemoryObject* v1, DistributedMemoryObject* v2, DistributedMemoryObject* v3, DistributedMemoryObject* v4, Complex alpha2, bool Ddag, bool inFourierSpace);
  
  void executeProjectorHatMultiplication(Complex* input, Complex* output, bool projPlus, bool daggered); //OHNE Faktor 0.5
  void executeProjectorMultiplication(Complex* input, Complex* output, bool projPlus); //OHNE Faktor 0.5

  
  
  void executeFermionMatrixStaticInverseMultiplication(Complex* input, Complex* output, bool dob, bool inFourierSpace);
  void executeQPreconditionerMatrixMultiplication(Complex* input, Complex* output, bool inverse, bool inFourierSpace);
  void executeRPreconditionerMatrixMultiplication(Complex* input, Complex* output, bool dob, bool inverse, bool inFourierSpace);
  void executeDistributedRPreconditionerMatrixMultiplication(DistributedMemoryObject* input, DistributedMemoryObject* output, bool dob, bool inverse, bool inFourierSpace);

  
  double checkLGSsolutionAccuracy(Complex* solution, Complex* result, double* phiField);  
  bool solveFermionMatrixLGS(Complex* input, Complex* output, double* phiField, double TOL, bool doubleFermionMatrix, bool quasiHermiteanMode, int maxIter, int& neededIter);
  bool executeFermionMatrixDaggerInverseMatrixMultiplication(Complex* input, Complex* output, double* phiField, double TOL, bool quasiHermiteanMode, int maxIter, int& neededIter);


  bool applyFermionMatrixMMDaggerFunction(Complex* input, Complex* output, double* phiField, double (*func)(double x),double TOL, int maxIter, int& neededIter, int auxVeccount, Complex** auxVecs, bool quasiHermiteanMode, bool inFourierSpace);
  bool applyDistributedFermionMatrixMMDaggerFunction(DistributedMemoryObject* input, DistributedMemoryObject* output, DistributedMemoryObject* phiField, double (*func)(double x),double TOL, int maxIter, int& neededIter, int auxVeccount, DistributedMemoryObject** auxVecs, bool quasiHermiteanMode, bool inFourierSpace);
  void applyFermionMatrixMMDaggerChebyshevPolynomial(Complex* input, Complex* output, double* phiField, GeneralChebyshevApproximation* chebyPoly, bool quasiHermiteanMode, bool inFourierSpace);
  void applyDistributedFermionMatrixMMDaggerChebyshevPolynomial(DistributedMemoryObject* input, DistributedMemoryObject* output, DistributedMemoryObject* phiField, GeneralChebyshevApproximation* chebyPoly, bool quasiHermiteanMode, bool inFourierSpace);


  void executeStarTransformationT(Complex* input, Complex* output, bool adjoining, bool inverse);
  void executeGamma0(Complex* input, Complex* output);
  void executeGamma5(Complex* input, Complex* output);
  void executeTau2(Complex* input, Complex* output);
  void executePiModeRemoverOperator(Complex* input, Complex* output, bool inFourierSpace);
  void executeDistributedPiModeRemoverOperator(DistributedMemoryObject* input, DistributedMemoryObject* output, bool inFourierSpace);

  void executeDistributedYBOperatorMultiplication(DistributedMemoryObject* input, DistributedMemoryObject* output, double xtrFac, DistributedMemoryObject* phiObj);
  void executeDistributedYBDaggerOperatorMultiplication(DistributedMemoryObject* input, DistributedMemoryObject* output, double xtrFac, DistributedMemoryObject* phiObj);


  void executeMultiplicationVectorWithDerivativesOfB(Complex* leftInput, Complex* rightInput, Complex* output);
  void executeDistributedMultiplicationVectorWithDerivativesOfB(DistributedMemoryObject* leftInput, DistributedMemoryObject* rightInput, DistributedMemoryObject* output);
  void executeMultiplicationVectorWithDerivativesOfMMdaggerInverse(Complex* omega, Complex* outMMdaggerInverseOmega, double* output, double* phiField, double TOL);
  void constructNeubergerWithXiFermionMatrix(ComplexMatrix& FermionMatrix, vector4D* phiField);
  void constructNeubergerWithXiFermionMatrix(ComplexMatrix& FermionMatrix, vector4D* phiField, double massSplitF);
  void constructWilsonWithOutXiFermionMatrix(ComplexMatrix& FermionMatrix, vector4D* phiField);
  void constructNeubergerWithOutXiFermionMatrix(ComplexMatrix& FermionMatrix, vector4D* phiField, double massSplitF);  
  void constructNeubergerWithOutXiFermionMatrix(ComplexMatrix& FermionMatrix, vector4D* phiField);
  void constructWilsonWithXiFermionMatrix(ComplexMatrix& FermionMatrix, vector4D* phiField);
  void exactFermionMatrixConditionNumber(double* phiField, double &eigMin, double &eigMax, double &invCond, bool doubleM, int probeCount, bool quasiHermiteanMode);
  int getVectorLength();
  int getVectorLengthXtrSize();
  int get1DSizeL0();
  int get1DSizeL1();
  int get1DSizeL2();
  int get1DSizeL3();  
  int get1DSizeLargest(); 
  int getTimeDirection(); 
  double getYN();
  void clearTrafoPlans();
  fftw_plan getFFTPlan(Complex* input, Complex* output, bool forw);
  void performFFT(Complex* input, Complex* output, bool forw);
  void performDistributedFFT(DistributedMemoryObject* input, DistributedMemoryObject* output, bool forw);  
  void getPreconditionerParameter(bool& use, double& m, double& s);
  void getQPreconditionerParameter(bool& use, double& mu, double& beta);  
  void getRPreconditionerParameter(bool& use, double& m, double& f);    
  void printPreconditionerParameter();
  void transformToXtraSizeArray(Complex* input, Complex* output);  
  void transformFromXtraSizeArray(Complex* input, Complex* output);
  void setxFFTusage(bool xFFTu);
  bool getxFFTusage();
  void setxFFT_DistributedFFT_ThreadCount(int xFFTThreadCnt);
  int getxFFT_DistributedFFT_ThreadCount();
  void tuneDistributedFFT(int tuneLevel);
  void tuneDistributedFFT(DistributedMemoryObject* input, DistributedMemoryObject* output, int tuneLevel);
  void tuneDistributedFFT(DistributedMemoryObject* input, DistributedMemoryObject* output, int tuneLevel, char* &fftPlanDescriptor);
  bool setDistributedFFTPlan(char* fftPlanDescriptor);
  
  
  void testFourierTrafo(bool forw);  
  void testDistributedFourierTrafo(bool forw);  
  Complex* calcFermionMatrixARPACKEigenValues(int what, int nev,  double* phiField, double TOL, bool doubleFermionMatrix, Complex** rightEigenVectors, bool quasiHermiteanMode, bool inFourierSpace);
};

#endif
