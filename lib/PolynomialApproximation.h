#ifndef PolynomialApproximation_included
#define PolynomialApproximation_included

#include <errno.h>
//#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <fstream>
#include <time.h>
#include <sys/stat.h> 
#include <unistd.h>
#include <cmath>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/shm.h>
#include <signal.h>
#include <cln/number.h>
#include <cln/float.h>
#include <cln/complex.h>
#include <cln/real.h>
#include <cln/io.h>
#include <cln/integer_io.h>
#include <cln/float_io.h>
#include <cln/complex_io.h>
#include <cln/univpoly.h>
#include <cln/univpoly_integer.h>
#include <cln/univpoly_complex.h>
#include <cln/univpoly_real.h>


#include "Global.h"
#include "Complex.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ios;

/**
*  This class is based on hep-lat/0302025.
**/
class PolynomialApproximation {
protected:
  int ID;
  int DIGIT; 
  int subPolyCount;
  int* polyDegree;
  bool Mass1Zero;
  bool Mass2Zero;  
  cln::float_format_t clnDIGIT;  
  cln::cl_F  ONE;
  cln::cl_F  TWO;
  cln::cl_F  ZERO;
  cln::cl_F  HALF;
  cln::cl_F  HUND;
  cln::cl_F  ALPHA;
  cln::cl_F  WEIGHTEXP;  
  cln::cl_F  Mass1;
  cln::cl_F  Mass2;
  
  cln::cl_F* EPSILON;
  cln::cl_F* LAMBDA;
  cln::cl_F** Recb;
  cln::cl_F** Recg;
  cln::cl_F** Orth;
  cln::cl_F** Coed;
  
  cln::cl_F** Coef;
  cln::cl_N** Roots;
  cln::cl_N** RootRoots;
  cln::cl_F* RootRepTotalNorm;
  cln::cl_F* RootRepSingleNorm;
  

  void ini(int subPolCnt, int id, int* deg, int digit, double alpha, double M1, double M2, double* eps, double* lam);
  void BaseIntS(int Maxpow, cln::cl_F Epsilon, cln::cl_F Lambda, cln::cl_F* Sint);
  void BaseIntT(int Maxpow, cln::cl_F Epsilon, cln::cl_F Lambda, cln::cl_F* Tint);
  cln::cl_F Recurev(int Maxpow, cln::cl_F Xval, cln::cl_F* Recb,cln::cl_F* Recg,cln::cl_F* Orth,cln::cl_F* Coed);
  void WriteAssign(std::ostream& Ostream, char* Text, cln::cl_F* Wlist, int Leng, char* Arrayname);
  void WriteRoots(char * filename, cln::cl_N *list, const int length);
  void WriteCoefs(cln::cl_N *list, const int length);
  void Quadropt(int Maxpow, cln::cl_F Epsilon, cln::cl_F Lambda, cln::cl_F* Recb, cln::cl_F* Recg, cln::cl_F* Orth, cln::cl_F* Coed, char* Filename);
  void RestorePolyCoef(int Maxpow, cln::cl_F* Recb, cln::cl_F* Recg, cln::cl_F* Orth, cln::cl_F* Coed, cln::cl_F* Coef);
  cln::cl_N EvalPoly(cln::cl_N* Poly, int Maxpow, cln::cl_N Valu);
  cln::cl_N EvalPolyProd(cln::cl_N* roots, const int Maxpow, const cln::cl_N x);
  cln::cl_N Lasolv(cln::cl_N* Poly, int Maxpow, cln::cl_N root, const int itemax=100);
  void Polyrootc(cln::cl_N* Poly, int Maxpow, cln::cl_N* Root);
  void Optimord(int Maxpow, cln::cl_N* Root, cln::cl_F Coef0, cln::cl_F Epsilon, cln::cl_F Lambda, int Ncheck);
  void OptimordPair(int Maxpow, cln::cl_N* Root, cln::cl_F Coef0, cln::cl_F Epsilon, cln::cl_F Lambda, int Ncheck);
  bool CompSortNaiv(const cln::cl_N & x, const cln::cl_N & y);
  int bitReversalRepresentation(const int idx, const int digits);
  void NaiveOrder(cln::cl_N * Roots, const int Maxpow);
  void BitReversalOrder(cln::cl_N *Roots, const int Maxpow);
  void BitReversalOrderOfRootPairs(cln::cl_N *Roots, const int Maxpow);
  void ApproxiRootr(int Maxpow, cln::cl_F Epsilon, cln::cl_F Lambda, cln::cl_F* Coef, cln::cl_N* Roots, cln::cl_N* RootRoots, cln::cl_F& RootRepTotalNorm, cln::cl_F& RootRepSingleNorm, char* Filename1,char* Filename2);
  void quicksort(const int n, cln::cl_N arr[], int idx[]);
  void calcApproxPolys(); 
  void enforceHermiticity(cln::cl_N *Roots, const int Maxpow);

/*  void ChebyCoeff(const int order, const cln::cl_F &epsilon, const cln::cl_F &lambda, cln::cl_F * Coeff, cln::cl_F * c);
  template <class T> cln::cl_F func(T &x);
  cln::cl_N EvalCheby(const int order, cln::cl_N * Coeff, cln::cl_N &x, const cln::cl_N & epsilon, const cln::cl_N lambda);*/
 
public:
  PolynomialApproximation(int subPolCnt, int id, int* deg, int digit, double alpha, double M1, double M2, double* eps, double* lam);
  void plotApproxPolynomials();
  Complex* getApproxPolyRoots(int polynomSlot);  
  int getApproxPolyDegree(int polynomSlot);
  double getPartPolyNormalization(int polynomSlot, int partPolyDeg);
  double evaluatePolynomial(double x, int polynomSlot);
};

#endif
