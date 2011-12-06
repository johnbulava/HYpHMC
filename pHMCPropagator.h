#ifndef pHMCPropagator_included
#define pHMCPropagator_included

#include <math.h>
#include "Complex.h"
#include "Quat.h"
#include "Global.h"
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include "Tools.h"
#include "FermionMatrixOperations.h"
#include "Propagator.h"
#include "pHMCForce.h"
#include "PolynomialApproximation.h"


#ifdef UseMPI
  #include <mpi.h>
#endif 

#define pHMCPropagator_UpperEWboundLogMAX 10000


// X-Werte duerfen negativ werden bis -N. Nach oben keine Schranke.
#define pHMCPropagator_SETPOLYROOTS 1
#define pHMCPropagator_RESTOREFORCE 3
#define pHMCPropagator_SAMPLEOMEGAMOMENTA 4
#define pHMCPropagator_PROPAGATEOMEGA 5
#define pHMCPropagator_PROPAGATEOMEGAMOMENTA 6
#define pHMCPropagator_CALCOMEGAACTION 7
#define pHMCPropagator_CALCOMEGAMOMENTUMACTION 8
#define pHMCPropagator_SAVEFORCE 9
#define pHMCPropagator_SETNU 10
#define pHMCPropagator_WRITEOMEGAFIELDTODISK 11
#define pHMCPropagator_READOMEGAFIELDFROMDISK 12
#define pHMCPropagator_MULTIPLYOMEGAMOMENTA 13
#define pHMCPropagator_SETTHETA 14
#define pHMCPropagator_WRITEOMEGAFORCESTRENGTHTODISK 15
#define pHMCPropagator_READOMEGAFORCESTRENGTHFROMDISK 16
#define pHMCPropagator_SETOMEGAMASSADAPTIONMODE 17
#define pHMCPropagator_CALCOMEGAMASSADAPTION 18
#define pHMCPropagator_CLEARACTSTORE 19
#define pHMCPropagator_CLEARTRASTARTSTORE 20
#define pHMCPropagator_COPYTRASTARTSTORETOACTSTORE 21
#define pHMCPropagator_COPYACTSTORETOTRASTARTSTORE 22
#define pHMCPropagator_STORECURRENTTOPLEVELFORCES 23
#define pHMCPropagator_ACTIVATEFORCESTORING 24
#define pHMCPropagator_CALCEXACTMMDAGINVERSESQRTOMEGAACTION 25
#define pHMCPropagator_RESETEXACTMMDAGINVERSESQRTOMEGAACTION 26
#define pHMCPropagator_GETAVERAGEPRECSUBPOL0OMEGAFORCESTRENGTH 27
#define pHMCPropagator_RESETAVERAGEOMEGAFORCESTRENGTH 28
#define pHMCPropagator_SAMPLEOMEGAFIELDS 29


#define pHMCPropagator_SETRANDOMSEED 50
#define pHMCPropagator_CHECKRANDOMSYNCHRON 51
#define pHMCPropagator_SETPHI 52
#define pHMCPropagator_RESTOREPROP 53
#define pHMCPropagator_SAMPLEPHIMOMENTA 54
#define pHMCPropagator_PROPAGATEPHI 55
#define pHMCPropagator_SAVEPROP 56
#define pHMCPropagator_DRAWRANDOMNUMBERSONSLAVENODES 57
#define pHMCPropagator_MULTIPLYPHIMOMENTA 58
#define pHMCPropagator_SYNCQPRECDATA 59
#define pHMCPropagator_SYNCRPRECDATA 60
#define pHMCPropagator_SETQUASIHERMITEANMODE 61
#define pHMCPropagator_SETMODELSELECTION 62
#define pHMCPropagator_SETTUNEMODE 63




class pHMCPropagator: public Propagator  {
private:

  char* fileNameIdentifier;
  
  pHMCForce** pHMCforces;
  Complex*** polyRoots;
  PolynomialApproximation** approxPoly;  
  
  int subPolyCount;
  double* polyEpsilon;
  double* polyLambda;  
  int* polyDegree;
  int precMassCount;
  double* precMasses;
  int polyDigit;
  double polyAlpha;
  int maxPolyDegreePerNode;
  double globalForceTheta;
  double syncRandom;
  bool nodesReady;
  bool quasiHermiteanMode;
  int AdditionalAuxVectorsCount;
  int upperEWboundLogCount;
  double** upperEWboundLog;
  
  double calcOmegaAction();
  double calcOmegaMomentumAction();
  void setPolyData();

  void ThreadController(int nr, int mode, int& para_int, double& para_double, double* data);
  void threadedExecute(int mode, double eps);
  void threadedExecutePROP(int mode, double eps);
  void threadedExecuteFORCE(int mode, double eps);
  void setNf(int nf);
  void iniAdditionalFields();
  void desiniAdditionalFields();
  int getFermionForceSubCategoryCount();
  int getFermionForceSubCategory(double x);  
  bool LeapOmelyanMarkovStep(int iterations, double epsilon, double propTOL, double finalTOL, double lambda, double rho, double theta, double mu);
  bool LeapOmelyanMultiTimeScaleMarkovStep(int level, int* iterations, double epsilon, int* subPolNr, double* propTOL, double finalTOL, double* lambda, double* rho, double* theta, double* mu);
  void LeapOmelyanMultiTimeScalePropagation(double* dSdPhi, int level, int* iterations, double epsilon, int* subPolNr, double* propTOL, double finalTOL, double* lambda, double* rho, double* theta, double* mu);
  void LeapOmelyanMultiTimeScalePropagation(int level, int* iterations, double epsilon, int* subPolNr, double* propTOL, double finalTOL, double* lambda, double* rho, double* theta, double* mu);
  void distributeMonomialsToForces();
  char* buildForceOmegaFieldSaveFileName(int nr);
  char* buildOmegaForceStrengthSaveFileName(bool PREC);
  void phiFieldWasChanged();
  void PreconditionerWasChanged();
  void clearForceStorage();
  void changeOfFACCtype();
public:
  
  pHMCPropagator(FermionMatrixOperations* fOps, double lam, double kap, double current, double c6, double c8, double c10, double lam6, double lam8, double lam10, int nf, double gam, bool sphMode, double sphZeta, double tht, int subPolCnt, double* polEps, double* polLam, int* polDeg, int precMCnt, double* precMss, int digit, double alpha, int maxPolDegPerNod, int addAuxVecsCount); 
  ~pHMCPropagator();

  void setFileNameIdentifier(char* identifier);
  void improvePreconditioningParameters(int iterGrob, int iterFein, double testTOL);
  void improveRPreconditioningParameters(int thermStepID, int ParameterAdaptionMode, double PolLambda, double upperEWsafetyFac);  
  double calcTotalAction(double finalTOL);
  void SlaveController();
  void setTheta(double tht);
  
  void getNodesReady();
  void saveALLfields();
  void restoreALLfields(bool clearStore);
  void sampleALLMomenta();
  void sampleOmegaFields();
  void writeAllOmegaFieldsToDisk();
  void readAllOmegaFieldsFromDisk();
  void negateAllMomenta();
  double getSyncSingleRandom();
  void writeOmegaForceStrengthsToDiskPREC();
  void writeOmegaForceStrengthsToDiskGLOBAL();
  void readOmegaForceStrengthsFromDiskPREC();
  void readOmegaForceStrengthsFromDiskGLOBAL();
  void setOmegaMassAdaptionMode(int omMassAdapMode);  
  void calcOmegaMassAdaption();
  void analyzeAllPolynomialForces();
  void smoothStartPhiConfiguration();
  void activateForceStoring(bool act);
  void synchronizedChangeOfQPreconsitionerData(bool useqPrecon, double mu, double beta);
  void synchronizedChangeOfRPreconsitionerData(bool userPrecon, double m, double f);
  void resetExactMMdagInverseSQRTOmegaAction();
  void calcExactMMdagInverseSQRTOmegaAction();
  double getExactReweighingFactorFromMMdagInverseSQRTOmegaAction();
  void synchronizedChangeOfQuasiHermiteanMode(bool qHM);
  void synchronizedChangeOfModelSelection(int modelSel);
  void synchronizedChangeOfTuneMode(bool tM);
  void saveUpperEWboundLogToDisk(char* fileName);
  void loadUpperEWboundLogFromDisk(char* fileName);
  pHMCForce* getForce(int forceNr);
};


#endif
