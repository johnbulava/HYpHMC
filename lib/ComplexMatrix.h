#ifndef ComplexMatrix_included
#define ComplexMatrix_included

#include "Complex.h"
#include "Quat.h"
#include "Global.h"
#ifdef useBLAS
  #include CBLASIncludeFile
#endif
#include "ComplexVector.h"


class ComplexMatrix {
protected:
  void ini(int size);
  void ini(const ComplexMatrix& m);
  void desini();
  bool calcEigenvaluesAndEigenvectorsHelper(bool calcVr);
  
  
  
public:
  ComplexMatrix(int size);
  ComplexMatrix();
  ComplexMatrix(Quat q);
  ComplexMatrix(ComplexVector v);  
  ComplexMatrix(const ComplexMatrix& m);  
  ComplexMatrix(ComplexMatrix outer, ComplexMatrix inner);    
  ~ComplexMatrix();
  void resize(int size);
  Complex* matrixElements;
  Complex** matrix;
  Complex* eigenvalues;
  ComplexVector** rightEigenVectors;
  Complex determinant;
  int matrixSize;
  bool calcEigenvalues();
  bool calcEigenvaluesAndEigenvectors();
  bool calcDeterminant(); 
  void print();
  bool insertMatrix(ComplexMatrix& m, int zeile, int spalte);
	bool addMatrix(ComplexMatrix& m, int zeile, int spalte);
  ComplexMatrix getSubMatrix(int zeile,int spalte, int size);
  void dagger();
  void setZero();
  bool invert();
  Complex tres();
  
  
  ComplexMatrix operator + (ComplexMatrix);
  ComplexMatrix operator - (ComplexMatrix);
  ComplexMatrix operator * (ComplexMatrix);
  
  ComplexVector operator * (ComplexVector);
  ComplexMatrix& operator = (ComplexMatrix);

  
};

ComplexMatrix operator * (double A, ComplexMatrix B);
ComplexMatrix operator * (Complex A, ComplexMatrix B);

#endif
