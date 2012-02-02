#ifndef pHMCForce_included
#define pHMCForce_included

#include <math.h>
#include "Complex.h"
#include "Quat.h"
#include "Global.h"
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include "Tools.h"
#include "FermionMatrixOperations.h"
#include "Force.h"
#include "LatticeSiteBins.h"
#include "LatticeMomentumBins.h"
#include "GeneralChebyshevApproximation.h"


#define pHMCForcePolynomialsMAX 20
#define pHMCForceDirectOmegaSamplingRelAccuracy 1E-12


class pHMCForce : public Force{
protected:
  double* storedActOmegaAction;
  double* storedTraStartOmegaAction;
  double omegaMomentumAction;
  double exactOmegaMMdagInverseSQRTAction;
  Complex** approxpolyRoots;
  int* approxpolyDegree;
  double* approxpolyLambda;
  Complex** auxVectors;
  Complex** fastEvalVectors;
  int fastEvalVectorsCount;
  Complex* dSdOmega;
  bool* storedActForcesAvail;
  bool* storedTraStartForcesAvail;
  Complex** storedActdSdOmega;
  Complex** storedTraStartdSdOmega;
  Complex** storedActdSdPhi;
  Complex** storedTraStartdSdPhi;
  Complex* omegaMomenta;
  Complex* omegaOLD;
  Complex* omegaMomentaOLD;
  GeneralChebyshevApproximation** chebyApproxOfInverseSQRTofPolynomial;
  int AdditionalAuxVectorsCount;
  Complex** AdditionalAuxVectors;
  int VectorCollectionCount;
  Complex** VectorCollection;
  
  DistributedMemoryObject* Distributed_Phi;
  DistributedMemoryObject* Distributed_dSdPhiDummy;
  DistributedMemoryObject* Distributed_dSdPhiSum;
  DistributedMemoryObject* Distributed_omegaField;
  DistributedMemoryObject** Distributed_auxVectors;
  DistributedMemoryObject** Distributed_fastEvalVectors;
  DistributedMemoryObject** Distributed_AdditionalAuxVectors;
  DistributedMemoryObject** Distributed_VectorCollection;
  
  double theta;
  double* nu;

  bool omegaMassAdaptionMode;
  double* omegaForceStrengthSumPREC;
  double* omegaForceStrengthSqrSumPREC;
  int* omegaForceStrengthCountPREC;
  double* omegaForceStrengthSumGLOBAL;
  double* omegaForceStrengthSqrSumGLOBAL;
  int* omegaForceStrengthCountGLOBAL;    
  LatticeMomentumBins** omegaForceStrengthAnalysisPREC; 
  LatticeMomentumBins** omegaForceStrengthAnalysisGLOBAL; 
  double* MomentaMasses;
  double* avgOmegaForces;
  bool quasiHermiteanMode;
  bool bMatFactorizationMode;
  bool tuneMode;
  
  
  void iniFastEvalVectors(int degMax);
  void desiniFastEvalVectors();
  void iniAdditionalFields();
  void desiniAdditionalFields();
  Complex* calcMonomialApplication(int polynomSlot, Complex* input, Complex* output, double* phi, int startM, int endM, bool inFourierSpace);
  DistributedMemoryObject* calcDistributedMonomialApplication(int polynomSlot, DistributedMemoryObject* input, DistributedMemoryObject* output, DistributedMemoryObject* phi, int startM, int endM, bool inFourierSpace);
  Complex* calcAllForcesNotHermiteanMode(int polynomSlot, double* phi);
  Complex* calcAllForcesQuasiHermiteanMode(int polynomSlot, double* phi);
  Complex* calcDistributedAllForcesQuasiHermiteanMode(int polynomSlot, double* phi);  
  void setChebyApproxOfInverseSQRTofPolynomial(int polynomSlot);
  void setVectorCollection();
  
public:
  pHMCForce(FermionMatrixOperations* fOps, bool loc, double tht, int id, int addAuxVecsCount); 
  ~pHMCForce();

  Complex* calcAllForces(int polynomSlot, double* phi, double& S);
  double calcOmegaAction(int polynomSlot, double* phi);
  double calcOmegaMomentumAction();
  double calcExactMMdagInverseSQRTomegaAction(double* phi, double alpha, int &neededIter);
  void applyInverseSQRTPolynomialKRYLOV(Complex* input, Complex* output, double* phi, int polynomSlot);
  void applyInverseSQRTPolynomialCHEBYSHEV(Complex* input, Complex* output, double* phi, int polynomSlot);
  void applyDistributedInverseSQRTPolynomialCHEBYSHEV(DistributedMemoryObject* input, DistributedMemoryObject* output, DistributedMemoryObject* phi, int polynomSlot);
  void applyFermionDoubleMatrixMMPolPol(Complex* input, Complex* output, double* phi, bool inFourierSpace);
  bool applyInverseFermionDoubleMatrixMMPolPol(Complex* input, Complex* output, double* phi, double TOL, int maxIter, int& neededIter);

  void sampleOmegaMomenta();
  void sampleOmegaMomenta(Complex* data);  
  void sampleOmegaFieldsPurelyGaussian();    
  void sampleOmegaFields(double* phiField);  
  void sampleOmegaFields(double* phiField, Complex* data);
  void drawAndWasteOmegaMomentaRandomNumbers();
  void randomOmegaField();
  void saveOmegaField();
  void restoreOmegaField();
  void doOmegaStep(double eps);
  void doOmegaMomentumStep(double eps);
  double calcGaussianWeightFactor(double* phi, double TOL, int &Ncg, int &NmmdagApplications);
  
  double getActOmegaAction(int polySlot);
  void setActOmegaAction(int polySlot, double s);
  double getExactOmegaMMdagInverseSQRTAction();
  void setExactOmegaMMdagInverseSQRTAction(double s);
  void activateForceStoring(int count, int* activateIndices);
  void clearActStore();
  void clearTraStartStore();
  void copyTraStartStoreToActStore();
  void copyActStoreToTraStartStore();
  void storeCurrentTopLevelForcesToTraStartStore();
  double getOmegaMomentumAction();
  void setOmegaMomentumAction(double s);  
  void setApproxPolyRoots(int polynomSlot, Complex* roots, int deg, double n, double polLam);
  Complex* getApproxPolyRoots(int polynomSlot);
  int getApproxPolyDegree(int polynomSlot);
  double getApproxPolyLambda(int polynomSlot);
  GeneralChebyshevApproximation* getApproxChebyApproxOfInverseSQRTofPolynomial(int polynomSlot);
  void setNu(int polynomSlot, double n);
  double getNu(int polynomSlot);
  void setQuasiHermiteanMode(bool qHM);
  void setBMatFactorizationMode(bool bMatFac);
  double getTheta();
  void setTheta(double tht);
  void plotApproxPoly(int polynomSlot, double eps, double lam, char* filename);
  void plotChebyApproxOfInversePolynomialSQRT(int polynomSlot, double eps, double lam, char* filename);
  void plotApproxPolyRoots(int polynomSlot, char* filename);
  void multiplyOmegaMomenta(double fac);
  Complex* getOmegaMomenta();
  void setOmegaMassAdaptionMode(bool omMassAdapMode);
  void setTuneMode(bool tuneM);
  void calcOmegaMassAdaption();
  void writeOmegaForceStrengthToDiskPREC(char *fileName);
  void writeOmegaForceStrengthToDiskGLOBAL(char *fileName);
  void readOmegaForceStrengthFromDiskPREC(char *fileName);
  void readOmegaForceStrengthFromDiskGLOBAL(char *fileName);
  void resetOmegaForceStrengthPREC();
  void resetOmegaForceStrengthGLOBAL();
  void calcAverageOmegaForceStrengthPREC(int polynomSlot);  
  double* getAverageOmegaForceStrengthPREC();
};
 

#endif
