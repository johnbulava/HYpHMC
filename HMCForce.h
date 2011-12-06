#ifndef HMCForce_included
#define HMCForce_included

#include <math.h>
#include "Complex.h"
#include "Quat.h"
#include "Global.h"
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include "Tools.h"
#include "FermionMatrixOperations.h"
#include "Force.h"



class HMCForce : public Force{
protected:
  Complex* outMMdaggerInverseOmega;
  double omegaMMdaggerInverseOmegaScalarProduct;
  
  void iniAdditionalFields();
  void desiniAdditionalFields();

public:
  HMCForce(FermionMatrixOperations* fOps, bool loc); 
  ~HMCForce();

  void sampleOmegaField(vector4D* phi, Complex* GaussVector);
  void calcPhiForce(vector4D* phi, double TOL);
  void calcInverse(vector4D* phi, double TOL);
  void calcOmegaMMdaggerInverseOmegaScalarProduct();
  Complex* getMMdaggerInverseOmega();
  double getOmegaMMdaggerInverseOmegaScalarProduct();
  void setOmegaMMdaggerInverseOmegaScalarProduct(double s);
};

#endif
