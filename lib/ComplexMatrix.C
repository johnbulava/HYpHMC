#include "ComplexMatrix.h"
#include <stdio.h>

ComplexMatrix::ComplexMatrix(const ComplexMatrix& m) {
  ini(m);
}


ComplexMatrix::ComplexMatrix() {
  ini(1);
}


ComplexMatrix& ComplexMatrix::operator = (ComplexMatrix m) {
  desini();
  ini(m);

  return *this;
}


void ComplexMatrix::ini(int size) {
  if (size <= 0) size = 0;
//  printf("ComplexMatrix: Creating matrix with size %d...",size);
  matrixSize = size;
  matrixElements = new Complex[matrixSize*matrixSize];
  int I;
  for (I=0; I<matrixSize*matrixSize; I++) {
    matrixElements[I].x = 0;
    matrixElements[I].y = 0;    
  }
  matrix = new Complex*[matrixSize];
  for (I=0; I<matrixSize; I++) {
    matrix[I] = (Complex*)&(matrixElements[I*matrixSize]);
  }  
  eigenvalues = new Complex[matrixSize];
  for (I=0; I<matrixSize; I++) {
    eigenvalues[I].x = 0;
    eigenvalues[I].y = 0;
  }  
  determinant.x = 0;
  determinant.y = 0;
  rightEigenVectors = new ComplexVector*[matrixSize];
  for (I=0; I<matrixSize; I++) {
    rightEigenVectors[I] = NULL;
  }
//  printf("sucessfully.\n");
}


void ComplexMatrix::ini(const ComplexMatrix& m) {
  ini(m.matrixSize);
  determinant = m.determinant;
  matrixSize = m.matrixSize;
 
  int I;
  for (I=0; I<matrixSize*matrixSize; I++) {
    matrixElements[I] = m.matrixElements[I];
  }
  for (I=0; I<matrixSize; I++) {
    eigenvalues[I] = m.eigenvalues[I];
  }
  for (I=0; I<matrixSize; I++) {
    if (m.rightEigenVectors[I] != NULL) {
      rightEigenVectors[I] = new ComplexVector((*(m.rightEigenVectors[I])));
    }
  }
}


void ComplexMatrix::desini() {
  delete [] matrixElements;
  delete [] eigenvalues;
  delete [] matrix;
  for (int I=0; I<matrixSize; I++) {
    if (rightEigenVectors[I] != NULL) {
      delete rightEigenVectors[I];
    }
  }
  delete[] rightEigenVectors;
  matrixSize = 0;
}


ComplexMatrix::ComplexMatrix(int size) {
  ini(size);
}


void ComplexMatrix::resize(int size) {
  desini();
  ini(size);
}


void ComplexMatrix::dagger() {
  Complex dummy;
  int I,I2;
  for (I=0; I<matrixSize; I++) {
    for (I2=I; I2<matrixSize; I2++) {
      dummy = matrix[I][I2];
      matrix[I][I2] = adj(matrix[I2][I]);
      matrix[I2][I] = adj(dummy);
    }
  }
}


void ComplexMatrix::setZero() {
  int I;
  for (I=0; I<matrixSize*matrixSize; I++) {
    matrixElements[I].x = 0.0;
    matrixElements[I].y = 0.0;
  }
}


/**
* Constructs a 2x2 matrix qTheta. (not Theta BAR)
**/
ComplexMatrix::ComplexMatrix(Quat q) {
  ini(2);

  matrixElements[0].x = q.x0;
  matrixElements[0].y = -q.x3;
  matrixElements[3].x = q.x0;
  matrixElements[3].y = +q.x3;
 
  matrixElements[1].x = -q.x2;
  matrixElements[1].y = -q.x1;
  matrixElements[2].x = q.x2;
  matrixElements[2].y = -q.x1;
}


ComplexMatrix::ComplexMatrix(ComplexVector v) {
  ini(v.vectorSize);
  int I,I2;
  int count = 0;
  
  for (I=0; I<matrixSize; I++) {
    for (I2=0; I2<matrixSize; I2++) {
      matrixElements[count] = v.vectorElements[I] * adj(v.vectorElements[I2]);
      count++;
    }
  }
}


ComplexMatrix::ComplexMatrix(ComplexMatrix outer, ComplexMatrix inner) {  
  ini(outer.matrixSize * inner.matrixSize);
  
  int index = 0;
  for (int o0=0; o0<outer.matrixSize; o0++) {
    for (int i0=0; i0<inner.matrixSize; i0++) {
      for (int o1=0; o1<outer.matrixSize; o1++) {
        for (int i1=0; i1<inner.matrixSize; i1++) {
          matrixElements[index] = outer.matrix[o0][o1] * inner.matrix[i0][i1];
	  index++;
	}
      }
    }
  }
}


ComplexMatrix::~ComplexMatrix() {
  desini();
}


bool ComplexMatrix::calcEigenvaluesAndEigenvectorsHelper(bool calcVr) {
#ifdef useLAPACK
  if (matrixSize <=0) return false;
    
  long int N = matrixSize;
  long int lda = N;
  Complex* A = new Complex[N*N];  
  long int lwork = -1;
  Complex* work = new Complex[2*N];
  
  long int info = 0;
  char jobvl = 'N';
  char jobvr = 'N';
  Complex* vl = NULL;
  Complex* vr = NULL;
  long int ldvl = N;
  long int ldvr = N;
  
  if (calcVr) {
    jobvl = 'V';        //left and right are mixed up in LAPACK
    vl = new Complex[N*N];
  }
  double* rwork = new double[2*N];
  
  cblas_zcopy(matrixSize*matrixSize, matrixElements, 1, A, 1);
  zgeev_(&jobvl, &jobvr, &N, A, &lda, eigenvalues, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);
  
  lwork = (int) work[0].x;
  delete [] work;
  work = new Complex[lwork];

  zgeev_(&jobvl, &jobvr, &N, A, &lda, eigenvalues, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);

  if (calcVr) {
    for (int I=0; I<matrixSize; I++) {
      if (rightEigenVectors[I] != NULL) {
        delete rightEigenVectors[I];
      }
      rightEigenVectors[I] = new ComplexVector(matrixSize);
    }
    int count = 0;
    for (int I=0; I<matrixSize; I++) {
      for (int I2=0; I2<matrixSize; I2++) {
        rightEigenVectors[I]->vectorElements[I2] = vl[count];
        count++;
      }
    }
  }

  delete [] work;
  delete [] A;
  delete [] vl;
  delete [] vr;
  delete [] rwork;  

  if (info != 0) {
    printf("Error in calcEigenvaluesAndEigenvectorsHelper: Could not calculate eigenvalues and eigenvectors!\n");
    return false;
  }
  return true;
  #else
    printf("Error in calcEigenvaluesAndEigenvectorsHelper: Programm wurde ohne LAPACK-Bibliotheken kompiliert!!!\n");
    exit(0);
  #endif
}


bool ComplexMatrix::calcEigenvalues() {
  return calcEigenvaluesAndEigenvectorsHelper(false);
}


bool ComplexMatrix::calcEigenvaluesAndEigenvectors() {
  return calcEigenvaluesAndEigenvectorsHelper(true);
}


bool ComplexMatrix::calcDeterminant() {
  int I;
  if (matrixSize <=0) return false;
  
  bool b = calcEigenvalues();
  if (!b) return false;
  
  determinant.x = 1;
  determinant.y = 0;
  
  for (I=0; I<matrixSize; I++) {
    determinant = determinant * eigenvalues[I];
  }
  return true;
}


void ComplexMatrix::print() {
  int I,I2;
  
  for (I2=0; I2<matrixSize; I2++) {
    for (I=0; I<matrixSize; I++) {
      double r = matrix[I2][I].x;
      double i = matrix[I2][I].y;
      if (r>=0) printf(" ");
      printf("%1.3f",r);
      if (i>=0) printf("+");
      printf("%1.3fi ",i);
    }
    printf("\n");
  }    
}


ComplexMatrix ComplexMatrix::operator + (ComplexMatrix param) {
  if (param.matrixSize!=matrixSize) {
    printf("ComplexMatrix::operator +: ERROR: Matrix sizes do not match!\n");
    return ComplexMatrix(0);
  }
  int I;
  ComplexMatrix result(matrixSize);
  for (I=0; I<matrixSize*matrixSize; I++) {
    result.matrixElements[I] = matrixElements[I] + param.matrixElements[I];
  }
  return result;
}


ComplexMatrix ComplexMatrix::operator - (ComplexMatrix param) {
  if (param.matrixSize!=matrixSize) {
    printf("ComplexMatrix::operator -: ERROR: Matrix sizes do not match!\n");
    return ComplexMatrix(0);
  }
  int I;
  ComplexMatrix result(matrixSize);
  for (I=0; I<matrixSize*matrixSize; I++) {
    result.matrixElements[I] = matrixElements[I] - param.matrixElements[I];
  }
  return result;
}


ComplexMatrix ComplexMatrix::operator * (ComplexMatrix param) {
  if (param.matrixSize!=matrixSize) {
    printf("ComplexMatrix::operator *: ERROR: Matrix sizes do not match!\n");
    return ComplexMatrix(0);
  }
  ComplexMatrix result(matrixSize);
  
  #ifdef useBLAS
    cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, matrixSize, matrixSize,
                matrixSize, &ComplexUnity, matrixElements, matrixSize, param.matrixElements, matrixSize,
                &ComplexZero, result.matrixElements, matrixSize);
  #else
    printf("Programm wurde ohne BLAS-Bibliotheken kompiliert!!!\n");
    exit(0);
  #endif

  return result;
}


ComplexVector ComplexMatrix::operator * (ComplexVector param) {
  if (param.vectorSize!=matrixSize) {
    printf("ComplexMatrix::operator *: ERROR: Matrix sizes do not match!\n");
    return ComplexVector(0);
  }
  ComplexVector result(matrixSize);

  #ifdef useBLAS
    cblas_zgemv(CblasRowMajor,CblasNoTrans, matrixSize, matrixSize,
                &ComplexUnity, matrixElements, matrixSize, param.vectorElements, 1,
                &ComplexZero, result.vectorElements, 1);
  #else
    printf("Programm wurde ohne BLAS-Bibliotheken kompiliert!!!\n");
    exit(0);
  #endif

  return result;
}


ComplexMatrix ComplexMatrix::getSubMatrix(int zeile,int spalte, int size) {
  ComplexMatrix res(size);
  
  if ((zeile<0) || (spalte<0) || (zeile+size>matrixSize) || (spalte+size>matrixSize)) {
    return res;
  }
  int pos = zeile*matrixSize+spalte;
  int resPos = 0;
  int s,z;
  for (z=0; z<size; z++) {
    for (s=0; s<size; s++) {
      res.matrixElements[resPos] = matrixElements[pos];
      pos++;
      resPos++;
    }
    pos += matrixSize-size;
  }
  
  return res;
}


bool ComplexMatrix::insertMatrix(ComplexMatrix& m, int zeile, int spalte) {
//inline bool ComplexMatrix::insertMatrix(ComplexMatrix& m, int zeile, int spalte) {
  if ((spalte<0) || (zeile<0) || (spalte+m.matrixSize>matrixSize) || (zeile+m.matrixSize>matrixSize)) {
    return false;
  }
  int pos = zeile * matrixSize+spalte;
  int mpos = 0;
  int I,I2;
  
  for (I=0; I<m.matrixSize; I++) {
    for (I2=0; I2<m.matrixSize; I2++) {
      matrixElements[pos] = m.matrixElements[mpos];
      mpos++;
      pos++;
    }
    pos += matrixSize - m.matrixSize;    
  }
  return true;
}

//inline
bool ComplexMatrix::addMatrix(ComplexMatrix& m, int zeile, int spalte) {
  if ((spalte<0) || (zeile<0) || (spalte+m.matrixSize>matrixSize) || (zeile+m.matrixSize>matrixSize)) {
    return false;
  }
  int pos = zeile * matrixSize+spalte;
  int mpos = 0;
  int I,I2;
  
  for (I=0; I<m.matrixSize; I++) {
    for (I2=0; I2<m.matrixSize; I2++) {
      matrixElements[pos] = matrixElements[pos] + m.matrixElements[mpos];
      mpos++;
      pos++;
    }
    pos += matrixSize - m.matrixSize;    
  }
  return true;
}


ComplexMatrix operator * (double A, ComplexMatrix B) {
  ComplexMatrix res(B.matrixSize);
  int I;
  
  for (I=0; I<B.matrixSize*B.matrixSize; I++) {
    res.matrixElements[I] = A*B.matrixElements[I];  
  }

  return res;
}


ComplexMatrix operator * (Complex A, ComplexMatrix B) {
  ComplexMatrix res(B.matrixSize);
  int I;
  
  for (I=0; I<B.matrixSize*B.matrixSize; I++) {
    res.matrixElements[I] = A*B.matrixElements[I];  
  }

  return res;
}


Complex ComplexMatrix::tres() {
  Complex res(0,0);
  for (int I=0; I<matrixSize; I++) {
    res = res + matrix[I][I];
  }
  return res;
}


bool ComplexMatrix::invert() {
#ifdef useLAPACK
  if (matrixSize <=0) return false;
    
  long int N = matrixSize;  
  long int lda = N;
  Complex* A = matrixElements;  
  long int* ipiv = new long int[N];
  long int info = 0;  
  
  zgetrf_(&N, &N, A, &lda, ipiv, &info);
 
  if (info != 0) {
    printf("Error in ComplexMatrix::invert: Could not do LU factorization (Error = %ld)!\n",info);
    delete[] ipiv;
    return false;  
  }  
  
  long int lwork = -1;
  Complex* work = new Complex[2*N];
  
  zgetri_(&N, A, &lda, ipiv, work, &lwork, &info);    

  if (info != 0) {
    printf("Error in ComplexMatrix::invert: Could not determine optimal blocksize!\n");  
    delete[] ipiv;
    delete[] work;
    return false;
  }
  
  lwork = (int) work[0].x;
  delete[] work;
  work = new Complex[lwork];

  zgetri_(&N, A, &lda, ipiv, work, &lwork, &info);    

  delete [] work;
  delete [] ipiv;  

  if (info != 0) {
    printf("Error in calcEigenvaluesAndEigenvectorsHelper: Could not calculate eigenvalues and eigenvectors!\n");
    return false;
  }
  return true;
#else
  printf("Error in calcEigenvaluesAndEigenvectorsHelper: Programm wurde ohne LAPACK-Bibliotheken kompiliert!!!\n");
  exit(0);
#endif
}
