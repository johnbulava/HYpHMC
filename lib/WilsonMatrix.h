#ifndef WilsonMatrix_included
#define WilsonMatrix_included

#include "Global.h"
#ifdef useBLAS
  #include CBLASIncludeFile
#endif
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include "DiracMatrix.h"
#include "Tools.h"


class WilsonMatrix: public DiracMatrix {
protected:
  double r;
  void iniW(double rValue, int OneDimLSize0, int OneDimLSize1, int OneDimLSize2, int OneDimLSize3, int NCopies);  
  
public:
  WilsonMatrix();
  WilsonMatrix(double rValue, int OneDimLSizeL0, int OneDimLSizeL1, int OneDimLSizeL2, int OneDimLSizeL3, int NCopies);  
  WilsonMatrix(const WilsonMatrix& w);
  ~WilsonMatrix();
  Complex analyticalEigenvalue(vector4D p);
  WilsonMatrix& operator = (WilsonMatrix);
  
  double getR();
};

#endif
