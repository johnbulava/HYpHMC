#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include <pthread.h>

#include "Complex.h"
#include "Quat.h"
#include "Global.h"
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include "Tools.h"
#include "FermionMatrixOperations.h"
#include "ComplexPolynom.h"
#include "HighPrecisionComplex.h"
#include "HighPrecisionComplexPolynom.h"
#include "PolynomialApproximation.h"
#include "GeneralChebyshevApproximation.h"


PolynomialApproximation* mainPoly;
int mainPolyDegree = 220;
double mainPolyEpsilon = 0.0005;
//int mainPolyDegree = 18;
//double mainPolyEpsilon = 0.1;
double mainPolLambda = 1.0;
Complex* mainPolyRoots;

int coeffCount;
double alpha;  //Recurrence - Parameter for Chebyshev - Polynomials
double beta;   //Recurrence - Parameter for Chebyshev - Polynomials
double* Gamma;  //Chebyshev-Coefficients


double approximant(double x) {
//  double y = mainPoly->evaluatePolynomial(x, 0);
  double nu = mainPoly->getPartPolyNormalization(0, mainPolyDegree);
  Complex z(nu, 0.0);
  Complex c(x, 0.0);
  
  for (int I=0; I<mainPolyDegree; I++) {
    z = z * (c - mainPolyRoots[I]);  
  }
  double y = z.x;
  
  return 1 / sqrt(y);
}


void clearData() {
  coeffCount = 0;
  alpha = NaN;
  beta = NaN; 
  delete[] Gamma;
  Gamma = NULL;
}


double evaluatePolynomial(double x) {
  if (coeffCount <= 0) return 0;
  double res = 0;
    
  res = Gamma[0];
  if (coeffCount == 1) return res;
  res += 0.5*(alpha*x+beta)*Gamma[1];    
  if (coeffCount == 2) return res;
  double tim1 = 1.0;
  double ti = 0.5*(alpha*x+beta);
  for (int I=2; I<coeffCount; I++) {
    double tip1 = alpha*x*ti + beta*ti - tim1;
    tim1 = ti;
    ti = tip1;
    res += Gamma[I]*ti;    
  }

  return res;
}


void calcApproximationStandard(double (*func)(double x), double minX, double maxX, int degree) {
  if (LogLevel>2) printf("Calculating general Chebyshev-Approximation Polynomial by standard method with minX = %f, maxX = %f, and degree = %d\n",minX,maxX,degree);

  clearData();  
  coeffCount = degree+1;  
  Gamma = new double[coeffCount];
  for (int I=0; I<coeffCount; I++) {
    Gamma[I] = 0;
  }

  double A = (maxX-minX)/2.0;
  double B = A + minX;
  alpha = 2.0 / A;
  beta = -2.0*B/A;

  for (int I=0; I<coeffCount; I++) {
    Gamma[I] = 0;
    for (int I2=0; I2<coeffCount; I2++) {
      double x = cos(pi*(I2+0.5) / coeffCount);
      Gamma[I] += (*func)(A*x+B) * cos(pi*I*(I2+0.5) / coeffCount);
    }
    Gamma[I] *= 2.0/coeffCount;
  }
  Gamma[0] /= 2.0;
}




void calcApproximationKrylovBased(double (*func)(double x), double minX, double maxX, int degree) {
  if (LogLevel>2) printf("Calculating general Krylov-based Chebyshev-Approximation Polynomial with minX = %f, maxX = %f, and degree = %d\n",minX,maxX,degree);
  double TOL = 1E-10;
  
  clearData();
  
  coeffCount = degree+1;  
  if ((coeffCount % 2)==1) {
    printf("ERROR in calcApproximationKrylovBased\n");
    exit(0);
  }

  ComplexMatrix KrylovOp(coeffCount);
  KrylovOp.setZero();
    
  for (int I=0; I<coeffCount-1; I++) {
    KrylovOp.matrix[I][I+1].x = 0.5;
    KrylovOp.matrix[I+1][I].x = 0.5;
  }
  KrylovOp.matrix[1][0].x = 1.0;
      
  bool b = KrylovOp.calcEigenvaluesAndEigenvectors(); 
  if (!b) {
    printf("ERROR: GeneralChebyshevApproximation::calcApproximation cannot calculate eigenvalues!!!\n");
    exit(0);    
  }
    
  for (int I=0; I<coeffCount; I++) {
    if (abs(KrylovOp.eigenvalues[I].y) > TOL) {
      b = false;
      printf("%e\n",KrylovOp.eigenvalues[I].y);
    }
    for (int I2=0; I2<coeffCount; I2++) {
      if (abs(KrylovOp.rightEigenVectors[I]->vectorElements[I2].y) > TOL) {
      b = false;
      printf("%e\n",KrylovOp.rightEigenVectors[I]->vectorElements[I2].y);	  
    }
    }
  }
  if (!b) {
    printf("ERROR: GeneralChebyshevApproximation::calcApproximation imaginary part not zero!!!\n");
    exit(0);    
  }
    
  Gamma = new double[coeffCount];
  ComplexMatrix BasisTransformation(coeffCount);
  ComplexVector unitVector(coeffCount);
    
  for (int I=0; I<coeffCount; I++) {
    for (int I2=0; I2<coeffCount; I2++) {
      BasisTransformation.matrix[I2][I] = KrylovOp.rightEigenVectors[I]->vectorElements[I2];
    }
  }
    
  b = BasisTransformation.invert();
  if (!b) {
    printf("ERROR: GeneralChebyshevApproximation::Cannot invert BasisTransformation!!!\n");
    exit(0);    
  }
    
  unitVector.setZero();
  unitVector.vectorElements[0].x = 1.0;
    
  ComplexVector dummyVec(1);
  dummyVec = BasisTransformation * unitVector;

  double A = (maxX-minX)/2.0;
  double B = A + minX;
  for (int I=0; I<coeffCount; I++) {
    double f = (*func)(A*KrylovOp.eigenvalues[I].x+B);
    dummyVec.vectorElements[I] = f * dummyVec.vectorElements[I];
    if (abs(dummyVec.vectorElements[I].y) > TOL) {
      b = false;
      printf("%e %e\n",dummyVec.vectorElements[I].x,dummyVec.vectorElements[I].y);
    }
  }
  if (!b) {
    printf("ERROR: GeneralChebyshevApproximation::Imaginary Part not zero (II)!!!\n");
    exit(0);    
  }
    
  for (int I=0; I<coeffCount; I++) {
    Gamma[I] = 0;
  }
    
  for (int I=0; I<coeffCount; I++) {
    for (int I2=0; I2<coeffCount; I2++) {
      Gamma[I2] += dummyVec.vectorElements[I].x * KrylovOp.rightEigenVectors[I]->vectorElements[I2].x;      
    }
  }

  alpha = 2.0 / A;
  beta = -2.0*B/A;
}


void calcApproximationRemezUntuned(double (*func)(double x), double minX, double maxX, int degree) {
  if (LogLevel>2) printf("Calculating general Chebyshev-Approximation Polynomial by REMEZ method with minX = %f, maxX = %f, and degree = %d\n",minX,maxX,degree);

  clearData();
  
  coeffCount = degree+1;  

  int N = coeffCount-1;
  Gamma = new double[coeffCount];
  for (int I=0; I<coeffCount; I++) {
    Gamma[I] = 0;
  }

  double A = (maxX-minX)/2.0;
  double B = A + minX;
  alpha = 2.0 / A;
  beta = -2.0*B/A;

  double* x = new double[N+2];
  for (int I=0; I<=N+1; I++) {
    x[I] = cos(I*pi/(N+1));
  }

  ComplexMatrix mat(N+2);

  for (int I=0; I<=N+1; I++) {
    mat.matrix[I][0] = Complex(1.0, 0.0);
    if (N>0) mat.matrix[I][1] = Complex(x[I], 0.0);
    for (int I2=2; I2<=N; I2++) {
      mat.matrix[I][I2] = 2*x[I]*mat.matrix[I][I2-1] - mat.matrix[I][I2-2];
    }
    mat.matrix[I][N+1].x = 1.0;
    mat.matrix[I][N+1].y = 0.0;
    if ((I%2)==1) mat.matrix[I][N+1].x = -1.0;
  }

  bool succ = mat.invert();
  if (!succ) {
    printf("ERROR: GeneralChebyshevApproximation::Cannot invert Chebyschev-Matrix!!!\n");
    exit(0);    
  }

  ComplexVector b(N+2);
  for (int I=0; I<=N+1; I++) {
    b.vectorElements[I].x = (*func)(A*x[I]+B);
    b.vectorElements[I].y = 0;
  }
  delete[] x;

  ComplexVector res = mat * b;

  for (int I=0; I<coeffCount; I++) {
    Gamma[I] = res.vectorElements[I].x;
  }
}





double calcRelAccuracy(int scanPoints, double minX, double maxX) {
  double relAccuracy = 0;
  for (int I=0; I<scanPoints; I++) {
    double x = I*(maxX-minX)/scanPoints + minX;
    double f = approximant(x);
    double p = evaluatePolynomial(x);
    double r = (p-f) / f;
    if (abs(r)>relAccuracy) {
      relAccuracy = abs(r);
    }
  }
    
  return relAccuracy;
}






void appendPolynomialAccuracies(int degree, char* fileName) {
  FILE* file = fopen(fileName, "a");
  int scanPoints = 100000;
  double relAcc1 = NaN;
  double relAcc2 = NaN;
  double relAcc3 = NaN;
  
  
  calcApproximationStandard(&approximant, 0, mainPolLambda, degree);
  relAcc1 = calcRelAccuracy(scanPoints, 0, mainPolLambda);
  
//  calcApproximationKrylovBased(&approximant, 0, mainPolLambda, degree);
//  relAcc2 = calcRelAccuracy(scanPoints, 0, mainPolLambda);

  calcApproximationRemezUntuned(&approximant, 0, mainPolLambda, degree);
  relAcc3 = calcRelAccuracy(scanPoints, 0, mainPolLambda);
  
  fprintf(file, "%d %1.15e %1.15e %1.15e \n", degree, relAcc1, relAcc2, relAcc3);
  
  fclose(file);
}




int main(int argc,char **argv) {
  iniTools(5517);
  LogLevel = 3;
  
  
  mainPoly = new PolynomialApproximation(0, 0, &mainPolyDegree, 1000, 0.5, 0, 0, &mainPolyEpsilon, &mainPolLambda);
  mainPoly->plotApproxPolynomials();
  mainPolyRoots = mainPoly->getApproxPolyRoots(0);


  char* fileName = new char[1000];
  snprintf(fileName, 1000, "RemezComparison_Deg%d_Eps%1.2e_Lam%1.2e.dat",mainPolyDegree,mainPolyEpsilon,mainPolLambda);
  for (int I=670; I<700; I+=1) {
    appendPolynomialAccuracies(I, fileName);
  }
  

  GeneralChebyshevApproximation* generalCheb = new GeneralChebyshevApproximation();
  generalCheb->calcApproximation(&approximant, 0, 1, 1E-4, 100000);
  
  


  delete mainPoly;
}



