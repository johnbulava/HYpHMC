#ifndef Propagator_included
#define Propagator_included

#include <math.h>
#include "Complex.h"
#include "Quat.h"
#include "Global.h"
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include "Tools.h"
#include "FermionMatrixOperations.h"
#include "Force.h"
#include FFTWIncludeFile
#include "LatticeMomentumBins.h"

#ifdef UseMPI
  #include <mpi.h>
#endif 

//Types
struct NeighboursType {
  int neighbourIndex[8];
};

// X-Werte duerfen negativ werden bis -N. Nach oben keine Schranke.
#define LPOS(x0, x1, x2, x3) (((x0+L0) % L0)*L3*L2*L1 + ((x1+L1) % L1)*L3*L2 + ((x2+L2) % L2)*L3 + ((x3+L3) % L3))
#define Propagator_EXIT -1
#define Propagator_FORCE 2
#define Propagator_SETFOURIERACCELERATIONTYPE 100
#define Propagator_SETFOURIERACCELERATIONMASSES 101



class Propagator {
protected:
  FermionMatrixOperations* fermiOps;
  double lambda;
  double kappa;
  double explicitCurrent;
  double ModelExtension_c6;
  double ModelExtension_c8;
  double ModelExtension_c10;
  double ModelExtension_lambda6;
  double ModelExtension_lambda8;
  double ModelExtension_lambda10;  
  bool SphericalMode;
  double SphericalZeta;
  int Nf;
  int forceCount;
  NeighboursType* neighbours;
  vector4D* momentaField;
  vector4D* momentaFieldOLD;  
  vector4D* phiFieldOLD;
  vector4D* phiForceInternal;
  vector4D* vectorInterim; 
  double* MomentaMasses;
  double* MomentumSpaceActionCoefficients;
  fftw_plan phiForceInternalPlanForward;
  fftw_plan phiForceInternalPlanBackward;  
  fftw_plan vectorInterimPlanForward;
  fftw_plan vectorInterimPlanBackward;
  fftw_plan momentaToInterimPlanForward;
  fftw_plan interimToMomentaPlanBackward;  
  fftw_plan internalToInterimPlanForward;
  fftw_plan phiFieldToInterimPlanForward;
  
  
  
  int phiForceFourierType;
  double phiForceFourierPara;
  bool momentumMassesDetermined;
  bool tuneMode;
  LatticeMomentumBins* phiTotalForceAnalysisPREC;
  LatticeMomentumBins* phiTotalForceAnalysisGLOBAL;  
  LatticeMomentumBins* phiTotalSphericalProjectedForceAnalysisPREC;
  LatticeMomentumBins* phiTotalSphericalProjectedForceAnalysisGLOBAL;    
  LatticeMomentumBins** phiFermionForceAnalysisPREC;
  LatticeMomentumBins** phiFermionForceAnalysisGLOBAL;  
  LatticeMomentumBins* phiHiggsForceAnalysisPREC;
  LatticeMomentumBins* phiHiggsForceAnalysisGLOBAL;    
  LatticeMomentumBins* phiChangeAnalysisNoFACC;  
  LatticeMomentumBins* phiChangeAnalysisFACC;  
  Force** forces;
  
  void iniFields();
  void desiniFields();
  void desiniForceCalculators();
  void desini();
  void calcPhiDerivativesOfSPhi();
  void calcPhiDerivativesOfSPhi(double* dSdPhi);
  void calcFullPhiDerivatives(double propTol, int forceSelect);
  void calcFullPhiDerivatives(double* dSdPhi, double propTol, int forceSelect);
  void MiniPhiMomentumStep(double eps);
  void MiniPhiMomentumStep(double* dSdPhi, double eps);
  void MiniPhiStep(double eps);  
  void checkValidityOfSetting();
  double calcPhiAction();
  double calcPhiMomentumAction();
  double Sold;
  double SoldOLD;
  double gamma;
  void ini(FermionMatrixOperations* fOps, double lam, double kap, double current, double c6, double c8, double c10, double lam6, double lam8, double lam10, int nf, double gam, bool sphMode, double sphZeta);
  void iniPhiForceFourierData();
  double optimalMomentumMass(double avgForce, double traLength, double latMomsqr, int i0, int i1, int i2, int i3);
  void calcMomentumSpaceActionCoefficients();
  void multiplyPhiMomenta(double fac);
  void LeapOmelyanPhiPropagation(int iterations, double epsilon, double lambda, double rho, double theta, double mu);
  void LeapFrogPhiPropagation(int iterations, double epsilon);
  void OmelyanO2PhiPropagation(int iterations, double epsilon);
  void OmelyanO4PhiPropagation(int iterations, double epsilon);


  virtual void ThreadController(int nr, int mode, int& para_int, double& para_double, double* data) {};
  virtual void threadedExecute(int mode, double TOL) {};
  virtual void setNf(int nf) {};
  virtual void iniAdditionalFields() {};
  virtual void desiniAdditionalFields() {};
  virtual int getFermionForceSubCategoryCount() {return 1;};
  virtual int getFermionForceSubCategory(double x) {return 0;};  
  virtual bool LeapOmelyanMarkovStep(int iterations, double epsilon, double propTOL, double finalTOL, double lambda, double rho, double theta, double mu) { return false; };
  virtual void LeapOmelyanMultiTimeScalePropagation(int level, int* iterations, double epsilon, int* subPolNr, double* propTOL, double finalTOL, double* lambda, double* rho, double* theta, double* mu) { };
  virtual bool LeapOmelyanMultiTimeScaleMarkovStep(int level, int* iterations, double epsilon, int* subPolNr, double* propTOL, double finalTOL, double* lambda, double* rho, double* theta, double* mu) { return false; };
  virtual void phiFieldWasChanged() {};
  virtual void PreconditionerWasChanged() {};
  virtual void changeOfFACCtype() {};

public:
  vector4D* phiField;
  double Sact;
  double SbeforeProp;
  double SafterProp;

  Propagator(FermionMatrixOperations* fOps, double lam, double kap, double current, double c6, double c8, double c10, double lam6, double lam8, double lam10, int nf, double gam, bool sphMode, double sphZeta); 
  virtual ~Propagator() {};

  void samplePhiMomentumField();
  void setKappa(double kap);
  void setCurrent(double current);
  void setModelExtensionParameterC6(double c6);
  void setModelExtensionParameterC8(double c8);
  void setModelExtensionParameterC10(double c10);
  void setModelExtensionParameterLambda6(double lam6);
  void setModelExtensionParameterLambda8(double lam8);
  void setModelExtensionParameterLambda10(double lam10);
  void setLambda(double lam);
  void setSphericalMode(bool sphMode, double sphZeta);
  void setGamma(double gam);
  void measure(double& measurePhiNorm, double& measureStaggeredPhiNorm, double& avgNorm, double& sigmaNorm);
  void savePhiFields();
  void restorePhiFields();
  void killSlaves();  
  bool LeapFrogMarkovStep(int iterations, double epsilon, double propTOL, double finalTOL);
  bool OmelyanO2MarkovStep(int iterations, double epsilon, double propTOL, double finalTOL);
  bool OmelyanO4MarkovStep(int iterations, double epsilon, double propTOL, double finalTOL);
  void multiTimeScalePropagation(int LevelCount, int* iterations, double epsilon, int* integrators, int* subPolNr, double* propTOL, double finalTOL);
  bool multiTimeScaleMarkovStep(int LevelCount, int* iterations, double epsilon, int* integrators, int* subPolNr, double* propTOL, double finalTOL);


  
  void resetPhiTotalSphericalProjectedForceFourierComponentsPREC();
  void resetPhiTotalSphericalProjectedForceFourierComponentsGLOBAL();
  void resetPhiTotalForceFourierComponentsPREC();
  void resetPhiTotalForceFourierComponentsGLOBAL();
  void resetPhiFermionForceFourierComponentsPREC();
  void resetPhiFermionForceFourierComponentsGLOBAL();
  void resetPhiHiggsForceFourierComponentsPREC();
  void resetPhiHiggsForceFourierComponentsGLOBAL();
  void resetPhiChangeFourierComponentsNoFACC();
  void resetPhiChangeFourierComponentsFACC();
  void calcPhiMomentumMasses(double traLength);
  void calcPhiForceFourierComponents(int forceType, int fermionForceSubCat, bool takeOverFFTresult);
  void calcPhiForceFourierComponents(double* dSdPhi, int forceType, int fermionForceSubCat, bool takeOverFFTresult);
  void calcPhiChangeFourierComponents();
  void savePhiTotalSphericalProjectedForceFourierComponentsPREC(char* fileName);
  void savePhiTotalSphericalProjectedForceFourierComponentsGLOBAL(char* fileName);
  void savePhiTotalForceFourierComponentsPREC(char* fileName);
  void savePhiTotalForceFourierComponentsGLOBAL(char* fileName);  
  void savePhiHiggsForceFourierComponentsPREC(char* fileName);
  void savePhiHiggsForceFourierComponentsGLOBAL(char* fileName);  
  void savePhiFermionForceFourierComponentsPREC(char* fileName);
  void savePhiFermionForceFourierComponentsGLOBAL(char* fileName);  
  void savePhiChangeFourierComponentsNoFACC(char* fileName);
  void savePhiChangeFourierComponentsFACC(char* fileName);
  void writeMomentumMasses(char* filename);
  void loadPhiTotalSphericalProjectedForceFourierComponentsPREC(char* fileName);
  void loadPhiTotalSphericalProjectedForceFourierComponentsGLOBAL(char* fileName);
  void loadPhiTotalForceFourierComponentsPREC(char* fileName);
  void loadPhiTotalForceFourierComponentsGLOBAL(char* fileName);  
  void loadPhiHiggsForceFourierComponentsPREC(char* fileName);
  void loadPhiHiggsForceFourierComponentsGLOBAL(char* fileName);  
  void loadPhiFermionForceFourierComponentsPREC(char* fileName);
  void loadPhiFermionForceFourierComponentsGLOBAL(char* fileName);  
  void loadPhiChangeFourierComponentsNoFACC(char* fileName);
  void loadPhiChangeFourierComponentsFACC(char* fileName);
  void setPhiForceFourierType(int type, double para);
  void improvePreconditioningParametersFAST();
  void getCopyOfSavedPhi(double* copy);
  
  
  virtual void improvePreconditioningParameters(int iterGrob, int iterFein, double testTOL) {};
  virtual double calcTotalAction(double finalTOL) {return 0.0;};
  virtual void SlaveController() {};
  virtual void smoothStartPhiConfiguration() {};
};


#endif
