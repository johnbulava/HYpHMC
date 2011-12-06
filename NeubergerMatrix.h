#ifndef NeubergerMatrix_included
#define NeubergerMatrix_included

#include "Global.h"
#ifdef useBLAS
  #include CBLASIncludeFile
#endif
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include "DiracMatrix.h"
#include "Tools.h"


class NeubergerMatrix: public DiracMatrix {
protected:
  double r;
  double rho;
  void iniN(double rhoValue, double rValue, int OneDimLSizeL0, int OneDimLSizeL1, int OneDimLSizeL2, int OneDimLSizeL3, int NCopies);  
  
public:
  NeubergerMatrix();
  NeubergerMatrix(double rhoValue, double rValue, int OneDimLSizeL0, int OneDimLSizeL1, int OneDimLSizeL2, int OneDimLSizeL3, int NCopies);  
  NeubergerMatrix(const NeubergerMatrix& n);
  ~NeubergerMatrix();
  Complex analyticalEigenvalue(vector4D p);
  NeubergerMatrix& operator = (NeubergerMatrix);
  
  double getRho();
  double getR();
};

#endif
