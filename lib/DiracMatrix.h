#ifndef DiracMatrix_included
#define DiracMatrix_included

#include "Global.h"
#ifdef useBLAS
  #include CBLASIncludeFile
#endif
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include "Tools.h"


class DiracMatrix: public ComplexMatrix {
protected:
  int NestedCopies;
  int OneDimLatticeSizeL0;
  int OneDimLatticeSizeL1;
  int OneDimLatticeSizeL2;
  int OneDimLatticeSizeL3;
  void constructOperator();
  void iniD(int OneDimLSizeL0, int OneDimLSizeL1, int OneDimLSizeL2, int OneDimLSizeL3, int NCopies);  
  
public:
  DiracMatrix();
  DiracMatrix(int N);
  DiracMatrix(const DiracMatrix& d);
  virtual ~DiracMatrix() {}
  virtual Complex analyticalEigenvalue(vector4D p) { return Complex(0,0);}
  void analyticalEigenvectors(vector4D p, ComplexVector v[4]);
  
};

#endif
