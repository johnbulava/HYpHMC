#ifndef PiModeBoosterMatrix_included
#define PiModeBoosterMatrix_included

#include "Global.h"
#ifdef useBLAS
  #include CBLASIncludeFile
#endif
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include "DiracMatrix.h"
#include "Tools.h"


class PiModeBoosterMatrix: public DiracMatrix {
protected:
  double boost;
  void iniB(double b, int OneDimLSize, int NCopies);  
  
public:
  PiModeBoosterMatrix();
  PiModeBoosterMatrix(double b, int OneDimLSize, int NCopies);  
  PiModeBoosterMatrix(const PiModeBoosterMatrix& w);
  ~PiModeBoosterMatrix();
  Complex analyticalEigenvalue(vector4D p);
  PiModeBoosterMatrix& operator = (PiModeBoosterMatrix);
  
  double getBoost();
};

#endif
