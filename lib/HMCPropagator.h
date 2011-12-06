#ifndef HMCPropagator_included
#define HMCPropagator_included

#include <math.h>
#include "Complex.h"
#include "Quat.h"
#include "Global.h"
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include "Tools.h"
#include "FermionMatrixOperations.h"
#include "Propagator.h"
#include "HMCForce.h"

#ifdef UseMPI
  #include <mpi.h>
#endif 


// X-Werte duerfen negativ werden bis -N. Nach oben keine Schranke.
#define HMCPropagator_SAMPLE 1
#define HMCPropagator_INVERSION 3



class HMCPropagator: public Propagator  {
private:
  Complex* GaussVector;
  HMCForce** HMCforces;
  
  double calcOmegaAction(bool outMMdaggerInverseOmega_READY, double finalTOL);
  bool LeapOmelyanMarkovStep(int iterations, double epsilon, double propTOL, double finalTOL, double lambda, double rho, double theta, double mu);

  void ThreadController(int nr, int mode, int& para_int, double& para_double, double* data);
  void threadedExecute(int mode, double TOL);
  void setNf(int nf);
  void sampleGaussVector();
  void iniAdditionalFields();
  void desiniAdditionalFields();
  int getFermionForceSubCategoryCount();
  int getFermionForceSubCategory(double x);  
  
public:
  HMCPropagator(FermionMatrixOperations* fOps, double lam, double kap, int nf, double gam); 
  ~HMCPropagator();

  void sampleOmegaFields();

  void improvePreconditioningParameters(int iterGrob, int iterFein, double testTOL);
  double calcTotalAction(double finalTOL);
  void SlaveController();
};

#endif
