#include "GeneralChebyshevApproximation.h"
#include "Tools.h"

GeneralChebyshevApproximation::GeneralChebyshevApproximation() { 
  coeffCount = 0;
  alpha = NaN;
  beta = NaN; 
  gamma = NULL;
  relAccuracy = NaN;
  
  clearData();
}



GeneralChebyshevApproximation::~GeneralChebyshevApproximation() { 
  clearData();
}


void GeneralChebyshevApproximation::clearData() {
  coeffCount = 0;
  alpha = NaN;
  beta = NaN; 
  delete[] gamma;
  gamma = NULL;
  relAccuracy = NaN;
}


double GeneralChebyshevApproximation::evaluatePolynomial(double x, bool withRenorm) {
  if (coeffCount <= 0) return 0;
  double res = 0;
    
  if (!withRenorm) {
    res = gamma[0];
    if (coeffCount == 1) return res;
    res += x*gamma[1];    
    if (coeffCount == 2) return res;
    double tim1 = 1.0;
    double ti = x;
    for (int I=2; I<coeffCount; I++) {
      double tip1 = 2*x*ti - tim1;
      tim1 = ti;
      ti = tip1;
      res += gamma[I]*ti;    
    }
  } else {
    res = gamma[0];
    if (coeffCount == 1) return res;
    res += 0.5*(alpha*x+beta)*gamma[1];    
    if (coeffCount == 2) return res;
    double tim1 = 1.0;
    double ti = 0.5*(alpha*x+beta);
    for (int I=2; I<coeffCount; I++) {
      double tip1 = alpha*x*ti + beta*ti - tim1;
      tim1 = ti;
      ti = tip1;
      res += gamma[I]*ti;    
    }
  }

  return res;
}


double GeneralChebyshevApproximation::evaluatePolynomial(double x) {
  return evaluatePolynomial(x);
}


void GeneralChebyshevApproximation::calcApproximation(double (*func)(double x), double minX, double maxX, double relAcc, int scanPoints) {
/*  int maxCoeffCount = 2000;
  int approxType = -1;
  bool b = calcApproximationRemez(func, minX, maxX, relAcc, scanPoints, maxCoeffCount);
  if (b && (maxCoeffCount>0)) approxType = 1;
  clearData();
  int dummy = maxCoeffCount;
  b = calcApproximationStandard(func, minX, maxX, relAcc, scanPoints, dummy);
  if (b && (dummy>0) && (dummy<maxCoeffCount)) {
    approxType = 2;
    maxCoeffCount = dummy;
  }
  clearData();
  b = calcApproximationKrylovBased(func, minX, maxX, relAcc, scanPoints, dummy);
  if (b && (dummy>0) && (dummy<maxCoeffCount)) {
    approxType = 3;
    maxCoeffCount = dummy;
  }
  clearData();

  maxCoeffCount = -1;
  if (approxType == 1) {
    b = calcApproximationRemez(func, minX, maxX, relAcc, scanPoints, maxCoeffCount);
    if (LogLevel > 2) printf("Final polynom constructed via Remez-Algorithm.\n");
  } else if (approxType == 2) {
    b = calcApproximationStandard(func, minX, maxX, relAcc, scanPoints, maxCoeffCount);
    if (LogLevel > 2) printf("Final polynom constructed via standard Algorithm.\n");
  } else if (approxType == 3) {
    b = calcApproximationKrylovBased(func, minX, maxX, relAcc, scanPoints, maxCoeffCount);
    if (LogLevel > 2) printf("Final polynom constructed via Krylov based Algorithm.\n");
  } else {
    printf("ERROR in calcApproximation: Could not find working approximation scheme!\n");
    exit(0);
  }
  if (!b) {
    printf("ERROR in calcApproximation: Could not find working approximation scheme!\n");
    exit(0);
  }*/
  
  int maxCoeffCount = 4000;
  bool b = calcApproximationRemez(func, minX, maxX, relAcc, scanPoints, maxCoeffCount);
  if (b && (maxCoeffCount>0)) {
    if (LogLevel > 2) printf("Final polynom constructed via Remez-Algorithm.\n");  
  } else {
    clearData();
    maxCoeffCount = 2000;
    bool b2 = calcApproximationStandard(func, minX, maxX, relAcc, scanPoints, maxCoeffCount);
    if (b2 && (maxCoeffCount>0)) {
      if (LogLevel > 2) printf("Final polynom constructed via standard Algorithm.\n");      
    } else {
      printf("ERROR in calcApproximation: Could not find working approximation scheme!\n");
      exit(0);  
    }
  } 
}


bool GeneralChebyshevApproximation::calcApproximationKrylovBased(double (*func)(double x), double minX, double maxX, double relAcc, int scanPoints, int &maxCoeffCount) {
  if (LogLevel>2) printf("Calculating general Krylov-based Chebyshev-Approximation Polynomial with minX = %f, maxX = %f, relAcc = %1.2e, scanPoints = %d, and maxCoeffCount = %d\n",minX,maxX,relAcc,scanPoints, maxCoeffCount);
  relAccuracy = 1E10;
  double TOL = 1E-10;
  
  int counter = 0;
  int counterAdd = 40;
  while ((counterAdd > 1) || (relAccuracy > relAcc)) {
    counter += counterAdd;
    clearData();
  
    coeffCount = 2*counter;  
    if ((maxCoeffCount >= 0) && (coeffCount > maxCoeffCount+counterAdd)) {
      if (LogLevel>2) printf("This scheme is inferior to prior attempts! Aborting.\n");
      return false;
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
    
    gamma = new double[coeffCount];
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
      gamma[I] = 0;
    }
    
    for (int I=0; I<coeffCount; I++) {
      for (int I2=0; I2<coeffCount; I2++) {
        gamma[I2] += dummyVec.vectorElements[I].x * KrylovOp.rightEigenVectors[I]->vectorElements[I2].x;      
      }
    }

    alpha = 2.0 / A;
    beta = -2.0*B/A;
    
    relAccuracy = 0;
    for (int I=0; I<scanPoints; I++) {
      double x = I*(maxX-minX)/scanPoints + minX;
      double f = (*func)(x);
      double p = evaluatePolynomial(x, true);
      double r = (p-f) / f;
      if (abs(r)>relAccuracy) {
        relAccuracy = abs(r);
      }
    }
    
    if (LogLevel>2) printf("  ...polynomial with degree = %d ==> relAcc = %1.2e\n",coeffCount-1,relAccuracy);   
    if ((counterAdd > 1) && (relAccuracy <= relAcc)) {
      counter -= counterAdd;
      counterAdd /= 6;
      if (counterAdd<=1) counterAdd = 1;
      relAccuracy = relAcc + 1;
    }
  } 
  maxCoeffCount = coeffCount;
  return true;
}


bool GeneralChebyshevApproximation::calcApproximationStandard(double (*func)(double x), double minX, double maxX, double relAcc, int scanPoints, int &maxCoeffCount) {
  if (LogLevel>2) printf("Calculating general Chebyshev-Approximation Polynomial by standard method with minX = %f, maxX = %f, relAcc = %1.2e, and scanPoints = %d, and maxCoeffCount = %d\n",minX,maxX,relAcc,scanPoints, maxCoeffCount);
  relAccuracy = 1E10;
  
  int counter = 0;
  int counterAdd = 50;
  while ((counterAdd > 1) || (relAccuracy > relAcc)) {
    counter += counterAdd;
    clearData();
  
    coeffCount = counter;  
    if ((maxCoeffCount >= 0) && (coeffCount > maxCoeffCount+counterAdd)) {
      if (LogLevel>2) printf("This scheme is inferior to prior attempts! Aborting.\n");
      return false;
    }
    gamma = new double[coeffCount];
    for (int I=0; I<coeffCount; I++) {
      gamma[I] = 0;
    }

    double A = (maxX-minX)/2.0;
    double B = A + minX;
    alpha = 2.0 / A;
    beta = -2.0*B/A;

    for (int I=0; I<coeffCount; I++) {
      gamma[I] = 0;
      for (int I2=0; I2<coeffCount; I2++) {
        double x = cos(pi*(I2+0.5) / coeffCount);
        gamma[I] += (*func)(A*x+B) * cos(pi*I*(I2+0.5) / coeffCount);
      }
      gamma[I] *= 2.0/coeffCount;
    }
    gamma[0] /= 2.0;

    relAccuracy = 0;
    for (int I=0; I<scanPoints; I++) {
      double x = I*(maxX-minX)/scanPoints + minX;
      double f = (*func)(x);
      double p = evaluatePolynomial(x, true);
      double r = (p-f) / f;
      if (abs(r)>relAccuracy) {
        relAccuracy = abs(r);
      }
    }
    
    if (LogLevel>2) printf("  ...polynomial with degree = %d ==> relAcc = %1.2e\n",coeffCount-1,relAccuracy);   
    if ((counterAdd > 1) && (relAccuracy <= relAcc)) {
      counter -= counterAdd;
      counterAdd /= 7;
      if (counterAdd<=1) counterAdd = 1;
      relAccuracy = relAcc + 1;
    }
  } 
  maxCoeffCount = coeffCount;
  return true;
}


bool GeneralChebyshevApproximation::calcApproximationRemez(double (*func)(double x), double minX, double maxX, double relAcc, int scanPoints, int &maxCoeffCount) {
  if (LogLevel>2) printf("Calculating general Chebyshev-Approximation Polynomial by REMEZ method with minX = %f, maxX = %f, relAcc = %1.2e, scanPoints = %d, and maxCoeffCount = %d\n",minX,maxX,relAcc,scanPoints, maxCoeffCount);
  relAccuracy = 1E10;
  
  int counter = 0;
  int counterAdd = 200;
  while ((counterAdd > 1) || (relAccuracy > relAcc)) {
    counter += counterAdd;
    clearData();
  
    coeffCount = counter;  
    if ((maxCoeffCount >= 0) && (coeffCount > maxCoeffCount+counterAdd)) {
      if (LogLevel>2) printf("This scheme is inferior to prior attempts! Aborting.\n");
      return false;
    }
    int N = coeffCount-1;
    gamma = new double[coeffCount];
    for (int I=0; I<coeffCount; I++) {
      gamma[I] = 0;
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
      gamma[I] = res.vectorElements[I].x;
    }

    relAccuracy = 0;
    for (int I=0; I<scanPoints; I++) {
      double x = I*(maxX-minX)/scanPoints + minX;
      double f = (*func)(x);
      double p = evaluatePolynomial(x, true);
      double r = (p-f) / f;
      if (abs(r)>relAccuracy) {
        relAccuracy = abs(r);
      }
    }
    
    if (LogLevel>2) printf("  ...polynomial with degree = %d ==> relAcc = %1.2e\n",coeffCount-1,relAccuracy);   
    if ((counterAdd > 1) && (relAccuracy <= relAcc)) {
      counter -= counterAdd;
      counterAdd /= 5;
      if (counterAdd<=1) counterAdd = 1;
      relAccuracy = relAcc + 1;
    }
  } 
  maxCoeffCount = coeffCount;  
  return true;
}


void GeneralChebyshevApproximation::plotPolynomial(char*fileName, double (*func)(double x), double minX, double maxX, int scanPoints) {
  FILE* file = fopen(fileName,"w");

  fprintf(file,"# Degree: %d\n", coeffCount-1);
  fprintf(file,"# relative accuracy(worst): %1.2e\n", relAccuracy);
  fprintf(file,"# alpha = %1.15f\n", alpha);
  fprintf(file,"# beta = %1.15f\n", beta);  
  for (int I=0; I<coeffCount; I++) {
    fprintf(file,"# Coefficient Gamma[%d] =  %1.15e\n", I, gamma[I]);
  }
  fprintf(file,"# \n");  
  
  for (int I=0; I<scanPoints; I++) {
    double x = I*(maxX-minX)/scanPoints + minX;
    double f = (*func)(x);
    double p = evaluatePolynomial(x, true);
    double r = (p-f) / f;
    fprintf(file,"%1.15f %1.15e %1.15e %1.15e\n",x,f,p,r);       
  }
  fclose(file);
}
 
  
int GeneralChebyshevApproximation::getCoeffCount() {
  return coeffCount;
}


double GeneralChebyshevApproximation::getAlpha() {
  return alpha;
}


double GeneralChebyshevApproximation::getBeta() {
  return beta;
}


double* GeneralChebyshevApproximation::getGammas() {
  return gamma;
}


double GeneralChebyshevApproximation::getRelAccuracy() {
  return relAccuracy;
}
