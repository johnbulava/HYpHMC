#ifndef SecondOrderBosonicEffectivePotential_included
#define SecondOrderBosonicEffectivePotential_included

#include <math.h>
#include "Global.h"
#include FFTWIncludeFile
#include "Complex.h"
#include "Tools.h"
#include "LatticeMomentumBins.h"
#include "Quat.h"
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include "FermionMatrixOperations.h"
#include "SimulationParameterSet.h"
#include "NeubergerMatrix.h"


#define SecondOrderBosonicEffectivePotential_MAXDERIVE 2

class SecondOrderBosonicEffectivePotential {
private:
  int L0,L1,L2,L3;
  int Vol;
  int perturbativeOrder;
  double lam0;
  double lam6;
  double lam8;
  double lam10; 
  double vev;
  double m0Sqr; 
  double YukT;
  double YukB;
  
  bool m0InDet;
  bool lambdasInDet;
  double* logBosDetDn;
  double* FermionBubble;

  double* pHatSqr;
  double** xSpaceHPropDn;
  double** xSpaceGPropDn;
  ComplexMatrix**** xSpaceFPropDn;
  double** pSpaceHPropDn;
  double** pSpaceGPropDn;
  ComplexMatrix**** pSpaceFPropDn;
  
  
  long int** CouplingTermCoeffcientCount;
  long int**** CouplingTermCoeffcients;  // N1 x N2 x termCount x (fac, doubleSum, v, GH0. GHxy, GG0, GGxy)
  int xSpacePropSumDBcount;
  int** xSpacePropSumDBid;
  double* xSpacePropSumDBres;
  

  fftw_plan FFTplan;
  Complex* FFTdummyArray;
  
  void calculatePHatSqr();
  void calculateProp();
  void calculateFermionProp(double rho, double r);  
  void calculateLogBosDet();
  void calculateFermionBubble();
  void calcCoefficientsOfCouplingTermRecursion(int n, int N, int* facLine, int &resultCount, long int* resultStore);
  void calcCoefficientsOfCouplingTermCombination(int N1, int N2);
  void evaluateCouplingTermDn(int N1, int N2, double* resDn, bool ignoreDisconnected);
  double numericalDerivativeOfEffPot_dv(double v, int order, int baseOrder, double h);
  

  
public:
  SecondOrderBosonicEffectivePotential(int l0, int l1, int l2, int l3, double v, double yt, double yb, double m0sqr, double la0, double la6, double la8, double la10, bool m0iD, bool lamInDet, int pertOrd); 
  ~SecondOrderBosonicEffectivePotential();
  
  void printCoefficientsOfCouplingTermCombination(int N1, int N2);
  void plotxSpacePropagators(int deriveNr);
  void plotpSpacePropagators(int deriveNr);

  void setLambdas(double la0, double la6, double la8, double la10); 
  void setVeV(double v);
  void setM0Sqr(double m0sqr);
  void setYukawa(double yt, double yb);
  
  
  double calcEffectivePotDn(int n);

};


#endif
