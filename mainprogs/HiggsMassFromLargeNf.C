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
#include "LatticeMomentumBins.h"
#include "SimulationParameterSet.h"
#include "SecondOrderBosonicEffectivePotential.h"


#define USEIMPROVEMENT_NO
#define USEBOSONICDETERMINANT_NO
#define BOSONICDETERMINANTPERTORDER 2
#define BOSONICDETERMINANTM0INDET 1
#define BOSONICDETERMINANTLAMBDASINDET 1
#define USEXTRTREAMMENTOFM0SQR
#define IMPROVEMENT_V4CORRECTION_NO
#define INCLUDEZEROMOMENTUM_NO
#define CALCULATEFERMIONMASSES_NO

#define INTEGRATORACCURACY 1E-10

//using namespace std;

struct ParameterType {
  double yt0;
  double yb0;
  double lam0;
  double lam6;
  double lam8;
  double lam10;  
  double CutoffInGev;
  int Nf;
  int L0;
  int L1;
  int L2;
  int L3;
  double rho;
  double r;
};
ParameterType  Parameters;
ParameterType  LastFermionContributionParameters;
ParameterType  LastBosonicContributionParameters;
double  LastBosonicContributionParameter_m0Sqr = NaN;
double LastFermionUContributionResult = NaN;
double LastFermionUprimeContributionResult = NaN;
double LastFermionUprimeprimeContributionResult = NaN;
double LastBosonicUContributionResult = NaN;
double LastBosonicUprimeContributionResult = NaN;
double LastBosonicUprimeprimeContributionResult = NaN;
ParameterType  LastPSqrParameters;
double* pSqr = NULL;
ParameterType  LastGoldstonePropagatorSumParameters;
double  LastGoldstonePropagatorSumResult = NaN;
ParameterType  LastHiggsPropagatorSumParameters;
double LastHiggsPropagatorSumParameter_MSqr = NaN;
double  LastHiggsPropagatorSumResult = NaN;
ParameterType  LastSquaredHiggsPropagatorSumParameters;
double LastSquaredHiggsPropagatorSumParameter_MSqr = NaN;
double  LastSquaredHiggsPropagatorSumResult = NaN;
NeubergerMatrix* diracOp = NULL;
vector4D calcFermionMassesExternalMomentumP;
ParameterType  LastDiracOperatorEigenvaluesAndEigenvectorsParameters;
Complex* DiracOpEigenvalues = NULL;
ComplexVector** DiracOpEigenvectors = NULL;
bool QuietMode = false;
int CouplingTermCoeffcientsCount[20];
int CouplingTermCoeffcients[20][20][5][4];  //Order Phi^N, term-count, derivative d/dv, term v,h,g,fac
SecondOrderBosonicEffectivePotential* secondOrderBosonicPotential = NULL;
ParameterType LastSecondOrderBosonicPotInitializationParameters;


double LambdaFunctionLuescherHelperK3(double p, double para) {
  return p*p*p*exp(-p*p-para*(p*p-1)*(p*p-1));
}
double LambdaFunctionLuescherHelperK5(double p, double para) {
  return p*p*p*p*p*exp(-p*p-para*(p*p-1)*(p*p-1));
}
double LambdaFunctionLuescherHelperK7(double p, double para) {
  return p*p*p*p*p*p*p*exp(-p*p-para*(p*p-1)*(p*p-1));
}
double LambdaFunctionLuescher(double x, double para) {
  if (isNaN(x)) return 1;
  double J3 = 0;
  double J5 = 0;
  double J7 = 0;
  int addCount = 1;
  
  int pos=0;
  double oldRes = 0;
  int count = addCount;
  while (count>0) {
    double res = J3 + integrate(&LambdaFunctionLuescherHelperK3, x, pos, pos+1, 1E-12);
    J3 = res;
    if (fabs((res-oldRes)/res)<1E-12) count--;
    oldRes = res;
    pos++;
  }

  pos=0;
  oldRes = 0;
  count = addCount;
  while (count>0) {
    double res = J5 + integrate(&LambdaFunctionLuescherHelperK5, x, pos, pos+1, 1E-12);
    J5 = res;
    if (fabs((res-oldRes)/res)<1E-12) count--;
    oldRes = res;
    pos++;
  }

  pos=0;
  oldRes = 0;
  count = addCount;
  while (count>0) {
    double res = J7 + integrate(&LambdaFunctionLuescherHelperK7, x, pos, pos+1, 1E-12);
    J7 = res;
    if (fabs((res-oldRes)/res)<1E-12) count--;
    oldRes = res;
    pos++;
  }

  return 3-2*J7*J3/(J5*J5);
}


double LambdaFunction1(double x, double para) {
  return x / (para + x);
}


double LambdaFunction2(double x, double para) {
  return  (exp(x/(para+0.001))-1) / exp(x/(para+0.001));
}


bool areParametersEqualWithRespectToSecondOrderBosonicPotInitialization(ParameterType p1, ParameterType p2) {
  if (p1.L0 != p2.L0) return false;
  if (p1.L1 != p2.L1) return false;
  if (p1.L2 != p2.L2) return false;
  if (p1.L3 != p2.L3) return false;

  return true;
}


void makeSecondOrderBosonicPotentialAvail(ParameterType p) {
  if ((secondOrderBosonicPotential == NULL) || (!areParametersEqualWithRespectToSecondOrderBosonicPotInitialization(p, LastSecondOrderBosonicPotInitializationParameters))) {
    delete secondOrderBosonicPotential;
//LogLevel=3;
    secondOrderBosonicPotential = new SecondOrderBosonicEffectivePotential(p.L0, p.L1, p.L2, p.L3, 1.0, p.yt0, p.yb0, 0.5, p.lam0, p.lam6, p.lam8, p.lam10, BOSONICDETERMINANTM0INDET, BOSONICDETERMINANTLAMBDASINDET, BOSONICDETERMINANTPERTORDER);
    LastSecondOrderBosonicPotInitializationParameters = p;
  }
}


bool areParametersEqualWithRespectToPSqr(ParameterType p1, ParameterType p2) {
  if (p1.L0 != p2.L0) return false;
  if (p1.L1 != p2.L1) return false;
  if (p1.L2 != p2.L2) return false;
  if (p1.L3 != p2.L3) return false;

  return true;
}


char* constructFileName(ParameterType p, char* des) {
  char* fileName = new char[1000];
  #ifdef USEBOSONICDETERMINANT
      snprintf(fileName, 1000,"HiggsMassFromLargeNfBOSONICDET_m0InDet%d_LamInDet%d_PertOrder%d_L%dx%dx%dx%dyt%1.3fyb%1.3flam%1.5fl6am%1.5fl8am%1.5fl10am%1.5fCut%1.1fNf%dRho%1.3fr%1.3f_%s.dat", BOSONICDETERMINANTM0INDET, BOSONICDETERMINANTLAMBDASINDET, BOSONICDETERMINANTPERTORDER, p.L0,p.L1,p.L2,p.L3, p.yt0,p.yb0, p.lam0,p.lam6,p.lam8,p.lam10,p.CutoffInGev,p.Nf,p.rho,p.r,des);  
  #else
    #ifdef USEIMPROVEMENT
      snprintf(fileName, 1000, "HiggsMassFromLargeNf_L%dx%dx%dx%dyt%1.3fyb%1.3flam%1.5fCut%1.1fNf%dRho%1.3fr%1.3f_%s.dat",p.L0,p.L1,p.L2,p.L3, p.yt0,p.yb0,p.lam0,p.CutoffInGev,p.Nf,p.rho,p.r,des);
    #else
      snprintf(fileName, 1000, "HiggsMassFromLargeNfUNIMPROVED_L%dx%dx%dx%dyt%1.3fyb%1.3flam%1.5fCut%1.1fNf%dRho%1.3fr%1.3f_%s.dat",p.L0,p.L1,p.L2,p.L3, p.yt0,p.yb0,p.lam0,p.CutoffInGev,p.Nf,p.rho,p.r,des);
    #endif
  #endif
  return fileName;  
}


char* constructFileNameShort(ParameterType p, char* des) {
  char* fileName = new char[1000];
  #ifdef USEBOSONICDETERMINANT
      snprintf(fileName, 1000,"HiggsMassFromLargeNfBOSONICDET_m0InDet%d_LamInDet%d_PertOrder%d_L%dx%dx%dx%dyt%1.3fyb%1.3fCut%1.1fNf%dRho%1.3fr%1.3f_%s.dat", BOSONICDETERMINANTM0INDET, BOSONICDETERMINANTLAMBDASINDET, BOSONICDETERMINANTPERTORDER,p.L0,p.L1,p.L2,p.L3, p.yt0,p.yb0,p.CutoffInGev,p.Nf,p.rho,p.r,des);  
  #else
    #ifdef USEIMPROVEMENT
      snprintf(fileName, 1000, "HiggsMassFromLargeNf_L%dx%dx%dx%dyt%1.3fyb%1.3fCut%1.1fNf%dRho%1.3fr%1.3f_%s.dat",p.L0,p.L1,p.L2,p.L3, p.yt0,p.yb0,p.CutoffInGev,p.Nf,p.rho,p.r,des);
    #else
      snprintf(fileName, 1000, "HiggsMassFromLargeNfUNIMPROVED_L%dx%dx%dx%dyt%1.3fyb%1.3fCut%1.1fNf%dRho%1.3fr%1.3f_%s.dat",p.L0,p.L1,p.L2,p.L3, p.yt0,p.yb0,p.CutoffInGev,p.Nf,p.rho,p.r,des);
    #endif
  #endif
  return fileName;  
}


double calcCouplingTermActionContributions(int N, int der, double v, double hProp, double gProp) {
  double res = 0;
  for (int I=0; I<CouplingTermCoeffcientsCount[N]; I++) {
    double dummy = 1;
    for (int I2=0; I2<CouplingTermCoeffcients[N][I][der][0]; I2++) dummy *= v;
    for (int I2=0; I2<CouplingTermCoeffcients[N][I][der][1]/2; I2++) dummy *= hProp;
    for (int I2=0; I2<CouplingTermCoeffcients[N][I][der][2]/2; I2++) dummy *= gProp;    
    
    res += CouplingTermCoeffcients[N][I][der][3] * dummy;
  }
  return res;
}


double calcCouplingTermActionContributionsImproved(int N, int derNr, double* derivatives, bool ignoreExactTerms) {
  int* tableCount = new int[derNr+1];
  int*** table = new int**[derNr+1];
  int colMax = 1+3*(1+derNr);
  
  //Initialize derivative table
  for (int I=0; I<1+derNr; I++) {
    tableCount[I] = 0;
    int lineMax = CouplingTermCoeffcientsCount[N];
    for (int I2=0; I2<I; I2++) lineMax *= (3+I2);
    table[I] = new int*[lineMax];
    for (int I2=0; I2<lineMax; I2++) {
      table[I][I2] = new int[colMax];
      for (int I3=0; I3<colMax; I3++) {
        table[I][I2][I3] = 0;
      }
    }
  }
  
  //Load base polynom
  for (int I=0; I<CouplingTermCoeffcientsCount[N]; I++) { 
    bool ignoreTerm = false;
    if ((ignoreExactTerms) && (CouplingTermCoeffcients[N][I][0][1]==0) && (CouplingTermCoeffcients[N][I][0][2]==0)) ignoreTerm = true;
    if ((ignoreExactTerms) && (CouplingTermCoeffcients[N][I][0][1]==2) && (CouplingTermCoeffcients[N][I][0][2]==0)) ignoreTerm = true;
    if ((ignoreExactTerms) && (CouplingTermCoeffcients[N][I][0][1]==0) && (CouplingTermCoeffcients[N][I][0][2]==2)) ignoreTerm = true;
    if (!ignoreTerm) {
      table[0][tableCount[0]][0] = CouplingTermCoeffcients[N][I][0][3];
      table[0][tableCount[0]][1] = CouplingTermCoeffcients[N][I][0][0];
      table[0][tableCount[0]][2] = CouplingTermCoeffcients[N][I][0][1]/2;
      table[0][tableCount[0]][3] = CouplingTermCoeffcients[N][I][0][2]/2;
      tableCount[0]++;
    }
  }
  
  //Calculate derivative table
  for (int I=1; I<1+derNr; I++) {
    for (int I2=0; I2<tableCount[I-1]; I2++) {
      for (int I3=1; I3<colMax; I3++) {
        if (table[I-1][I2][I3] > 0) {
          for (int I4=0; I4<colMax; I4++) {
 	    table[I][tableCount[I]][I4] = table[I-1][I2][I4];
	  }
	  table[I][tableCount[I]][0] *= table[I-1][I2][I3];
	  table[I][tableCount[I]][I3] -= 1;
	  table[I][tableCount[I]][I3+3] += 1;
	  tableCount[I]++;
	}
      }
    }
  }
  
  //Evaluate derivative table
  double res = 0;
  for (int I=0; I<tableCount[derNr]; I++) {
    double contrib = table[derNr][I][0];
    for (int I2=1; I2<colMax; I2++) {
      for (int I3=0; I3<table[derNr][I][I2]; I3++) {
        contrib *= derivatives[I2-1];
      }
    }
    res += contrib;
  } 

  //Desinitialize derivative table
  for (int I=0; I<1+derNr; I++) {
    int lineMax = CouplingTermCoeffcientsCount[N];
    for (int I2=0; I2<I; I2++) lineMax *= (3+I2);
    for (int I2=0; I2<lineMax; I2++) {
      delete[] table[I][I2];
    }
    delete[] table[I];
  }
  delete[] table;

  return res;
}


void makePSqrAvailableForParameters(ParameterType p) {
  if ((pSqr==NULL) || (!areParametersEqualWithRespectToPSqr(p, LastPSqrParameters))) {
    printf("Generating PSqr\n");
    delete[] pSqr;
    pSqr = new double[p.L0 * p.L1 * p.L2 * p.L3];
    LatticeMomentumBins* latBin = new LatticeMomentumBins(p.L0, p.L1, p.L2, p.L3);
    for (int I=0; I<p.L0 * p.L1 * p.L2 * p.L3; I++) {
      pSqr[I] = latBin->getLatMomSqrFromIndex(I);
    }
    delete latBin;
    LastPSqrParameters = p;  
  }
}


inline double sinPSqr(double p0, double p1, double p2, double p3) {
  return 4.0 * (sqr(sin(0.5*p0)) + sqr(sin(0.5*p1)) + sqr(sin(0.5*p2)) + sqr(sin(0.5*p3)));
}


double calcPropagatorSums_Integrand(double p0, double p1, double p2, double p3, double para) {
  return 1.0 / (sinPSqr(p0, p1, p2, p3) + para);
}


void calcPropagatorSums(ParameterType p, double mSqr, double &higgsSum, double &goldSum) {
  if ((p.L0>=0) && (p.L1>=0) && (p.L2>=0) && (p.L3>=0)) {
    makePSqrAvailableForParameters(p);
    if ((!isNaN(LastHiggsPropagatorSumResult)) && (LastHiggsPropagatorSumParameter_MSqr==mSqr) && (areParametersEqualWithRespectToPSqr(p, LastHiggsPropagatorSumParameters))) {
      higgsSum = LastHiggsPropagatorSumResult;
    } else {    
      higgsSum = 0;
      int startP = 1;
      #ifdef INCLUDEZEROMOMENTUM
        startP = 0;
      #endif
      
      for (int I=startP; I<p.L0 * p.L1 * p.L2 * p.L3; I++) {
        higgsSum += 1.0 / (pSqr[I] + mSqr);
      }
      higgsSum /= p.L0 * p.L1 * p.L2 * p.L3;
    
      LastHiggsPropagatorSumResult = higgsSum;
      LastHiggsPropagatorSumParameter_MSqr = mSqr;
      LastHiggsPropagatorSumParameters = p;
    }
    
    if ((!isNaN(LastGoldstonePropagatorSumResult)) && (areParametersEqualWithRespectToPSqr(p, LastGoldstonePropagatorSumParameters))) {
      goldSum = LastGoldstonePropagatorSumResult;
    } else {
      goldSum = 0;  
      for (int I=1; I<p.L0 * p.L1 * p.L2 * p.L3; I++) {
        goldSum += 1.0 / (pSqr[I] + 0);    
      }
      goldSum /= p.L0 * p.L1 * p.L2 * p.L3;  
    
      LastGoldstonePropagatorSumResult = goldSum;
      LastGoldstonePropagatorSumParameters = p;
    }
  } else {  
    vector4D start;
    vector4D end;
    for (int I=0; I<4; I++) {
      start[I] = 0;
      end[I] = pi;
    }
  
    if ((!isNaN(LastHiggsPropagatorSumResult)) && (LastHiggsPropagatorSumParameter_MSqr==mSqr) && (areParametersEqualWithRespectToPSqr(p, LastHiggsPropagatorSumParameters))) {
      higgsSum = LastHiggsPropagatorSumResult;
    } else {
      higgsSum = sqr(sqr(1.0/pi)) * integrate(&calcPropagatorSums_Integrand, mSqr, start, end, INTEGRATORACCURACY);
    
      LastHiggsPropagatorSumResult = higgsSum;
      LastHiggsPropagatorSumParameter_MSqr = mSqr;
      LastHiggsPropagatorSumParameters = p;
    }
    
    if ((!isNaN(LastGoldstonePropagatorSumResult)) && (areParametersEqualWithRespectToPSqr(p, LastGoldstonePropagatorSumParameters))) {
      goldSum = LastGoldstonePropagatorSumResult;
    } else {
      printf("Integrating Goldstone - Modes\n");
      goldSum = sqr(sqr(1.0/pi)) * integrate(&calcPropagatorSums_Integrand, 0, start, end, INTEGRATORACCURACY);  
    
      LastGoldstonePropagatorSumResult = goldSum;
      LastGoldstonePropagatorSumParameters = p;
    }
  }  
}


double calcSelfCouplingCorrectionCoefficientForEffectivePotential(ParameterType p, double mSqr) {
  double higgsSum = NaN;
  double goldSum = NaN;

  calcPropagatorSums(p, mSqr, higgsSum, goldSum);
  return 6*higgsSum + 6*goldSum;   // Vertauschungs-Multiplizitaet
}


double calcSquaredHiggsPropagatorSum_Integrand(double p0, double p1, double p2, double p3, double para) {
  double dummy = sinPSqr(p0, p1, p2, p3) + para;
  return 1.0 / (dummy*dummy);
}


double calcSquaredHiggsPropagatorSum(ParameterType p, double mSqr) {
  double higgsSum = NaN;
  if ((p.L0>=0) && (p.L1>=0) && (p.L2>=0) && (p.L3>=0)) {
    makePSqrAvailableForParameters(p);
    if ((!isNaN(LastSquaredHiggsPropagatorSumResult)) && (LastSquaredHiggsPropagatorSumParameter_MSqr==mSqr) && (areParametersEqualWithRespectToPSqr(p, LastSquaredHiggsPropagatorSumParameters))) {
      higgsSum = LastSquaredHiggsPropagatorSumResult;
    } else {    
      higgsSum = 0;
      int startP = 1;
      #ifdef INCLUDEZEROMOMENTUM
        startP = 0;
      #endif
      for (int I=startP; I<p.L0 * p.L1 * p.L2 * p.L3; I++) {
        higgsSum += 1.0 / ((pSqr[I] + mSqr) * (pSqr[I] + mSqr));

      }
      higgsSum /= p.L0 * p.L1 * p.L2 * p.L3;

      LastSquaredHiggsPropagatorSumResult = higgsSum;
      LastSquaredHiggsPropagatorSumParameter_MSqr = mSqr;
      LastSquaredHiggsPropagatorSumParameters = p;
    }    
  } else {  
    vector4D start;
    vector4D end;
    for (int I=0; I<4; I++) {
      start[I] = 0;
      end[I] = pi;
    }
  
    if ((!isNaN(LastSquaredHiggsPropagatorSumResult)) && (LastSquaredHiggsPropagatorSumParameter_MSqr==mSqr) && (areParametersEqualWithRespectToPSqr(p, LastSquaredHiggsPropagatorSumParameters))) {
      higgsSum = LastSquaredHiggsPropagatorSumResult;
    } else {
      higgsSum = sqr(sqr(1.0/pi)) * integrate(&calcSquaredHiggsPropagatorSum_Integrand, mSqr, start, end, INTEGRATORACCURACY);
    
      LastSquaredHiggsPropagatorSumResult = higgsSum;
      LastSquaredHiggsPropagatorSumParameter_MSqr = mSqr;
      LastSquaredHiggsPropagatorSumParameters = p;
    }
  }  
  return higgsSum;
}


double calcSelfCouplingV4CorrectionCoefficientForEffectivePotential(ParameterType p, double mSqr) {
  #ifdef IMPROVEMENT_V4CORRECTION
    return 1*calcSquaredHiggsPropagatorSum(p, mSqr);   // Vertauschungs-Multiplizitaet  
  #else
    return 1E-30;  
  #endif
}


bool areParametersEqualWithRespectToFermionContribution(ParameterType p1, ParameterType p2) {
  if (p1.L0 != p2.L0) return false;
  if (p1.L1 != p2.L1) return false;
  if (p1.L2 != p2.L2) return false;
  if (p1.L3 != p2.L3) return false;
  if (p1.yt0 != p2.yt0) return false;
  if (p1.yb0 != p2.yb0) return false;
  if (p1.CutoffInGev != p2.CutoffInGev) return false;
  if (p1.rho != p2.rho) return false;
  if (p1.r != p2.r) return false;

  return true;
}


double calcFermionUContribution_IntegrandU(double p0, double p1, double p2, double p3, double para) {
  vector4D k;  
  k[0] = p0;
  k[1] = p1;
  k[2] = p2;
  k[3] = p3;
    
  double fac = 0.5 / Parameters.rho;
  double v = Physical_VEV_GeV / Parameters.CutoffInGev;
  double yt = Parameters.yt0;
  double yb = Parameters.yb0;    
  Complex ew = diracOp->analyticalEigenvalue(k);
  
  Complex wt = ComplexUnity - (fac*ew);
  Complex zt = ew + yt*v*wt;
  Complex zb = ew + yb*v*wt;
  return -2.0*log(zt.x*zt.x+zt.y*zt.y) -2.0*log(zb.x*zb.x+zb.y*zb.y);  
}


double calcFermionUContribution_IntegrandUprime(double p0, double p1, double p2, double p3, double para) {
  vector4D k;  
  k[0] = p0;
  k[1] = p1;
  k[2] = p2;
  k[3] = p3;
    
  double fac = 0.5 / Parameters.rho;
  double v = Physical_VEV_GeV / Parameters.CutoffInGev;
  double yt = Parameters.yt0;
  double yb = Parameters.yb0;    
  Complex ew = diracOp->analyticalEigenvalue(k);
  
  Complex w = ComplexUnity - (fac*ew);
  Complex zt = ew + yt*v*w;
  Complex zb = ew + yb*v*w;
  return -4.0*(yt*w/zt).x -4.0*(yb*w/zb).x;
}


double calcFermionUContribution_IntegrandUprimeprime(double p0, double p1, double p2, double p3, double para) {
  vector4D k;  
  k[0] = p0;
  k[1] = p1;
  k[2] = p2;
  k[3] = p3;
    
  double fac = 0.5 / Parameters.rho;
  double v = Physical_VEV_GeV / Parameters.CutoffInGev;
  double yt = Parameters.yt0;
  double yb = Parameters.yb0;    
  double ytSqr = Parameters.yt0*Parameters.yt0;
  double ybSqr = Parameters.yb0*Parameters.yb0;  
  Complex ew = diracOp->analyticalEigenvalue(k);
  
  Complex w = ComplexUnity - (fac*ew);
  Complex zt = ew + yt*v*w;
  Complex zb = ew + yb*v*w;
  return 4.0*(ytSqr*w*w/(zt*zt)).x + 4.0*(ybSqr*w*w/(zb*zb)).x;
}


void calcFermionUContribution(ParameterType p, double &U, double &Uprime, double &Uprimeprime) {
  if ((p.yt0==0) && (p.yb0==0)) {
    U = 0;
    Uprime = 0;
    Uprimeprime = 0;
    return;
  }

  if ((!isNaN(LastFermionUContributionResult)) && (!isNaN(LastFermionUprimeContributionResult)) && (!isNaN(LastFermionUprimeprimeContributionResult))) {
    if (areParametersEqualWithRespectToFermionContribution(p, LastFermionContributionParameters)) {
      U = LastFermionUContributionResult;
      Uprime = LastFermionUprimeContributionResult;
      Uprimeprime = LastFermionUprimeprimeContributionResult;
      return;
    }
  }
  
  if (!QuietMode) printf("Calculating Fermion Contribution\n");
  U = 0;
  Uprime = 0;
  Uprimeprime = 0;
  vector4D k;  
  double fac = 0.5 / p.rho;
  double v = Physical_VEV_GeV / p.CutoffInGev;
  double ytSqr = p.yt0*p.yt0;
  double ybSqr = p.yb0*p.yb0;  
  
  if ((p.L0>=0) && (p.L1>=0) && (p.L2>=0) && (p.L3>=0)) {
    for (int i0=0; i0<p.L0; i0++) {
      k[0] = 2*pi*i0 / p.L0;
      for (int i1=0; i1<p.L1; i1++) {
        k[1] = 2*pi*i1 / p.L1;
        for (int i2=0; i2<p.L2; i2++) {
          k[2] = 2*pi*i2 / p.L2;
          for (int i3=0; i3<p.L3; i3++) {
            k[3] = 2*pi*i3 / p.L3;

            Complex ew = diracOp->analyticalEigenvalue(k);
	    //if ((fabs(fac*ew.x-1)>1E-14) || (fabs(ew.y)>1E-14)) {
  	     /* 
				Complex z = ew / (ComplexUnity - (fac*ew));
	      double nz = z.x*z.x + z.y*z.y;
	   
				//Also disagreement w/ Phillip's thesis here? 
  	      double at = ytSqr*v*v + nz;
	      double ab = ybSqr*v*v + nz;
	      U += -2*log(at) -2*log(ab);
	  
	      Uprime += -4*ytSqr*v/at  -4*ybSqr*v/ab;
	   
	      Uprimeprime += -4*ytSqr*(nz - ytSqr*v*v) / (at*at)  -4*ybSqr*(nz - ybSqr*v*v) / (ab*ab);
	   		*/
				U += calcFermionUContribution_IntegrandU(k[0], k[1], k[2], k[3], 0.0); 
				Uprime += calcFermionUContribution_IntegrandUprime(k[0], k[1], k[2], k[3], 0.0); 
				Uprimeprime += calcFermionUContribution_IntegrandUprimeprime(k[0], k[1], k[2], k[3], 0.0); 
				//Expression from his thesis
				//Complex tz = (z + p.yt0*v*ComplexUnity)*(ComplexUnity - (fac*ew));
				//double at2 = tz.x*tz.x + tz.y*tz.y; 

				//std::cout << "at = " << at << std::endl
				//std::cout << "at2 = " << at2 << std::endl; 
				//std::cout << "at - 2at2 = " << at - 2*at2 << std::endl<<std::endl; 

			
		//	} 
   	  }
        }
      }
    }
    U /= p.L0*p.L1*p.L2*p.L3;  
    Uprime /= p.L0*p.L1*p.L2*p.L3;  
    Uprimeprime /= p.L0*p.L1*p.L2*p.L3;  
  } else {
    printf("Using Integration\n");
    vector4D start;
    vector4D end;
    for (int I=0; I<4; I++) {
      start[I] = 0;
      end[I] = pi;
    }

    U = sqr(sqr(1.0/pi)) * integrate(&calcFermionUContribution_IntegrandU, NaN, start, end, INTEGRATORACCURACY);
    Uprime = sqr(sqr(1.0/pi)) * integrate(&calcFermionUContribution_IntegrandUprime, NaN, start, end, INTEGRATORACCURACY);
    Uprimeprime = sqr(sqr(1.0/pi)) * integrate(&calcFermionUContribution_IntegrandUprimeprime, NaN, start, end, INTEGRATORACCURACY);
  }
  
  U *= p.Nf;
  Uprime *= p.Nf;
  Uprimeprime *= p.Nf;
  
  LastFermionContributionParameters = p;
  LastFermionUContributionResult = U;
  LastFermionUprimeContributionResult = Uprime;
  LastFermionUprimeprimeContributionResult = Uprimeprime;
}


bool areParametersEqualWithRespectToBosonicContribution(ParameterType p1, ParameterType p2) {
  if (p1.L0 != p2.L0) return false;
  if (p1.L1 != p2.L1) return false;
  if (p1.L2 != p2.L2) return false;
  if (p1.L3 != p2.L3) return false;
  if (p1.CutoffInGev != p2.CutoffInGev) return false;
  if (p1.lam0 != p2.lam0) return false;
  if (p1.lam6 != p2.lam6) return false;
  if (p1.lam8 != p2.lam8) return false;
  if (p1.lam10 != p2.lam10) return false;

  return true;
}


double calcBosonicUContribution_IntegrandU(double p0, double p1, double p2, double p3, double para) {
  vector4D k;  
  k[0] = p0;
  k[1] = p1;
  k[2] = p2;
  k[3] = p3;
    
  double v = Physical_VEV_GeV / Parameters.CutoffInGev;
  double lam0 = Parameters.lam0;
  double lam6 = Parameters.lam6;
  double lam8 = Parameters.lam8;
  double lam10 = Parameters.lam10;
  double pS = sinPSqr(p0, p1, p2, p3);

  double zh   = lam0  * (12*v*v)
              + lam6  * (30*v*v*v*v)
  	      + lam8  * (56*v*v*v*v*v*v)
              + lam10 * (90*v*v*v*v*v*v*v*v);
  double zg   = lam0  * ( 4*v*v)
              + lam6  * ( 6*v*v*v*v)
	      + lam8  * ( 8*v*v*v*v*v*v)
	      + lam10 * (10*v*v*v*v*v*v*v*v);

  return 0.5*log(pS + para + zh) + 1.5*log(pS + para + zg);
}

      
double calcBosonicUContribution_IntegrandUprime(double p0, double p1, double p2, double p3, double para) {
  vector4D k;  
  k[0] = p0;
  k[1] = p1;
  k[2] = p2;
  k[3] = p3;
    
  double v = Physical_VEV_GeV / Parameters.CutoffInGev;
  double lam0 = Parameters.lam0;
  double lam6 = Parameters.lam6;
  double lam8 = Parameters.lam8;
  double lam10 = Parameters.lam10;
  double pS = sinPSqr(p0, p1, p2, p3);

  double zh   = lam0  * (12*v*v)
              + lam6  * (30*v*v*v*v)
   	      + lam8  * (56*v*v*v*v*v*v)
	      + lam10 * (90*v*v*v*v*v*v*v*v);
  double zg   = lam0  * ( 4*v*v)
              + lam6  * ( 6*v*v*v*v)
	      + lam8  * ( 8*v*v*v*v*v*v)
	      + lam10 * (10*v*v*v*v*v*v*v*v);
  double dzh  = lam0  * (24*v)
              + lam6  * (120*v*v*v)
	      + lam8  * (336*v*v*v*v*v)
	      + lam10 * (720*v*v*v*v*v*v*v);
  double dzg  = lam0  * ( 8*v)
              + lam6  * (24*v*v*v)
	      + lam8  * (48*v*v*v*v*v)
   	      + lam10 * (80*v*v*v*v*v*v*v);

  return 0.5*dzh / (pS + para + zh) + 1.5*dzg / (pS + para + zg);
}


double calcBosonicUContribution_IntegrandUprimeprime(double p0, double p1, double p2, double p3, double para) {
  vector4D k;  
  k[0] = p0;
  k[1] = p1;
  k[2] = p2;
  k[3] = p3;
    
  double v = Physical_VEV_GeV / Parameters.CutoffInGev;
  double lam0 = Parameters.lam0;
  double lam6 = Parameters.lam6;
  double lam8 = Parameters.lam8;
  double lam10 = Parameters.lam10;
  double pS = sinPSqr(p0, p1, p2, p3);

  double zh   = lam0  * (12*v*v)
              + lam6  * (30*v*v*v*v)
	      + lam8  * (56*v*v*v*v*v*v)
	      + lam10 * (90*v*v*v*v*v*v*v*v);
  double zg   = lam0  * ( 4*v*v)
              + lam6  * ( 6*v*v*v*v)
	      + lam8  * ( 8*v*v*v*v*v*v)
	      + lam10 * (10*v*v*v*v*v*v*v*v);
  double dzh  = lam0  * (24*v)
              + lam6  * (120*v*v*v)
	      + lam8  * (336*v*v*v*v*v)
	      + lam10 * (720*v*v*v*v*v*v*v);
  double dzg  = lam0  * ( 8*v)
              + lam6  * (24*v*v*v)
	      + lam8  * (48*v*v*v*v*v)
	      + lam10 * (80*v*v*v*v*v*v*v);
  double ddzh = lam0  * (24)
              + lam6  * (360*v*v)
	      + lam8  * (1680*v*v*v*v)
	      + lam10 * (5040*v*v*v*v*v*v);
  double ddzg = lam0  * (8)
              + lam6  * (72*v*v)
	      + lam8  * (240*v*v*v*v)
	      + lam10 * (560*v*v*v*v*v*v);  
		  
  return -0.5*(dzh*dzh - ddzh*(pS + para + zh)) / sqr(pS + para + zh) - 1.5*(dzg*dzg - ddzg*(pS + para + zg)) / sqr(pS + para + zg);
}


void calcBosonicUContribution(ParameterType p, double m0Sqr, double &U, double &Uprime, double &Uprimeprime) {
  if ((p.lam0==0) && (p.lam6==0) && (p.lam8==0) && (p.lam10==0)) {
    U = 0;
    Uprime = 0;
    Uprimeprime = 0;
    return;
  }

  if ((!isNaN(LastBosonicUContributionResult)) && (!isNaN(LastBosonicUprimeContributionResult)) && (!isNaN(LastBosonicUprimeprimeContributionResult)) && (!isNaN(LastBosonicContributionParameter_m0Sqr))) {
    if ((areParametersEqualWithRespectToBosonicContribution(p, LastBosonicContributionParameters) && (m0Sqr == LastBosonicContributionParameter_m0Sqr))) {
      U = LastBosonicUContributionResult;
      Uprime = LastBosonicUprimeContributionResult;
      Uprimeprime = LastBosonicUprimeprimeContributionResult;
      return;
    }
  }
  
  if (!QuietMode) printf("Calculating Bosonic-Determinant Contribution\n");
  U = 0;
  Uprime = 0;
  Uprimeprime = 0;
  double v = Physical_VEV_GeV / p.CutoffInGev;
 
  double xtrM0SqrSave = m0Sqr;
  double xtrTrigger = -1;
  #ifdef USEXTRTREAMMENTOFM0SQR
    xtrTrigger = 20.6;
  #endif

  if ((p.L0>=0) && (p.L1>=0) && (p.L2>=0) && (p.L3>=0)) {
    makePSqrAvailableForParameters(p);
    double Tg = 0;
    double Tgprime = 0;
    double Tgprimeprime = 0;
    double Th = 0;
    double Thprime = 0;
    double Thprimeprime = 0;
    
    double xtrTg = 0;
    double xtrTgprime = 0;
    double xtrTgprimeprime = 0;
    double xtrTh = 0;
    double xtrThprime = 0;
    double xtrThprimeprime = 0;
    

    for (int I=1; I<p.L0 * p.L1 * p.L2 * p.L3; I++) {
      if (pSqr[I]>xtrTrigger) {
        m0Sqr = xtrM0SqrSave;
      } else {
        m0Sqr = 0;
      }    
    
      double zh   = p.lam0  * (12*v*v)
                  + p.lam6  * (30*v*v*v*v)
	 	  + p.lam8  * (56*v*v*v*v*v*v)
		  + p.lam10 * (90*v*v*v*v*v*v*v*v);
      double zg   = p.lam0  * ( 4*v*v)
                  + p.lam6  * ( 6*v*v*v*v)
		  + p.lam8  * ( 8*v*v*v*v*v*v)
		  + p.lam10 * (10*v*v*v*v*v*v*v*v);
      double dzh  = p.lam0  * (24*v)
                  + p.lam6  * (120*v*v*v)
		  + p.lam8  * (336*v*v*v*v*v)
		  + p.lam10 * (720*v*v*v*v*v*v*v);
      double dzg  = p.lam0  * ( 8*v)
                  + p.lam6  * (24*v*v*v)
		  + p.lam8  * (48*v*v*v*v*v)
		  + p.lam10 * (80*v*v*v*v*v*v*v);
      double ddzh = p.lam0  * (24)
                  + p.lam6  * (360*v*v)
		  + p.lam8  * (1680*v*v*v*v)
		  + p.lam10 * (5040*v*v*v*v*v*v);
      double ddzg = p.lam0  * (8)
                  + p.lam6  * (72*v*v)
	          + p.lam8  * (240*v*v*v*v)
	          + p.lam10 * (560*v*v*v*v*v*v);

      U += 0.5*log(pSqr[I] + m0Sqr + zh);
      U += 1.5*log(pSqr[I] + m0Sqr + zg);
      Th += 1 / (pSqr[I] + m0Sqr + zh);
      Tg += 1 / (pSqr[I] + m0Sqr + zg);

      Uprime += 0.5*dzh / (pSqr[I] + m0Sqr + zh);
      Uprime += 1.5*dzg / (pSqr[I] + m0Sqr + zg);
      Thprime += -2*12*p.lam0*v / ((pSqr[I] + m0Sqr + zh)*(pSqr[I] + m0Sqr + zh));
      Tgprime += -2*4*p.lam0*v  / ((pSqr[I] + m0Sqr + zg)*(pSqr[I] + m0Sqr + zg));
      
      Uprimeprime += -0.5*(dzh*dzh - ddzh*(pSqr[I] + m0Sqr + zh)) / sqr(pSqr[I] + m0Sqr + zh);
      Uprimeprime += -1.5*(dzg*dzg - ddzg*(pSqr[I] + m0Sqr + zg)) / sqr(pSqr[I] + m0Sqr + zg);
      Thprimeprime += (8*12*p.lam0*v*p.lam0*v - 2*12*p.lam0*(pSqr[I] + m0Sqr + zg)) / ((pSqr[I] + m0Sqr + zh)*(pSqr[I] + m0Sqr + zh)*(pSqr[I] + m0Sqr + zh));      
      Tgprimeprime += (8*4*p.lam0*v*p.lam0*v - 2*4*p.lam0*(pSqr[I] + m0Sqr + zg))  / ((pSqr[I] + m0Sqr + zg)*(pSqr[I] + m0Sqr + zg)*(pSqr[I] + m0Sqr + zg));
      
      if (pSqr[I]<=xtrTrigger) {
        xtrTh += 1 / (pSqr[I] + m0Sqr + zh);
        xtrTg += 1 / (pSqr[I] + m0Sqr + zg);

        xtrThprime += -2*12*p.lam0*v / ((pSqr[I] + m0Sqr + zh)*(pSqr[I] + m0Sqr + zh));
        xtrTgprime += -2*4*p.lam0*v  / ((pSqr[I] + m0Sqr + zg)*(pSqr[I] + m0Sqr + zg));
      
        xtrThprimeprime += (8*12*p.lam0*v*p.lam0*v - 2*12*p.lam0*(pSqr[I] + m0Sqr + zg)) / ((pSqr[I] + m0Sqr + zh)*(pSqr[I] + m0Sqr + zh)*(pSqr[I] + m0Sqr + zh));      
        xtrTgprimeprime += (8*4*p.lam0*v*p.lam0*v - 2*4*p.lam0*(pSqr[I] + m0Sqr + zg))  / ((pSqr[I] + m0Sqr + zg)*(pSqr[I] + m0Sqr + zg)*(pSqr[I] + m0Sqr + zg));
      }
      
      m0Sqr = xtrM0SqrSave;
    }
    U /= p.L0*p.L1*p.L2*p.L3;  
    Uprime /= p.L0*p.L1*p.L2*p.L3;  
    Uprimeprime /= p.L0*p.L1*p.L2*p.L3;  
    Th /= p.L0*p.L1*p.L2*p.L3;
    Tg /= p.L0*p.L1*p.L2*p.L3;
    Thprime /= p.L0*p.L1*p.L2*p.L3;
    Tgprime /= p.L0*p.L1*p.L2*p.L3;
    Thprimeprime /= p.L0*p.L1*p.L2*p.L3;
    Tgprimeprime /= p.L0*p.L1*p.L2*p.L3;
    xtrTh /= p.L0*p.L1*p.L2*p.L3;
    xtrTg /= p.L0*p.L1*p.L2*p.L3;
    xtrThprime /= p.L0*p.L1*p.L2*p.L3;
    xtrTgprime /= p.L0*p.L1*p.L2*p.L3;
    xtrThprimeprime /= p.L0*p.L1*p.L2*p.L3;
    xtrTgprimeprime /= p.L0*p.L1*p.L2*p.L3;

    double derivatives[9];
    derivatives[0] = v;
    derivatives[1] = Th;
    derivatives[2] = Tg;
    derivatives[3] = 1;
    derivatives[4] = Thprime;
    derivatives[5] = Tgprime;
    derivatives[6] = 0;
    derivatives[7] = Thprimeprime;
    derivatives[8] = Tgprimeprime;

    double xtrU = p.lam0*calcCouplingTermActionContributionsImproved(4, 0, &(derivatives[0]), true);
    double xtrUprime = p.lam0*calcCouplingTermActionContributionsImproved(4, 1, &(derivatives[0]), true);
    double xtrUprimeprime = p.lam0*calcCouplingTermActionContributionsImproved(4, 2, &(derivatives[0]), true);
    
    xtrU += p.lam6*calcCouplingTermActionContributionsImproved(6, 0, &(derivatives[0]), true);
    xtrUprime += p.lam6*calcCouplingTermActionContributionsImproved(6, 1, &(derivatives[0]), true);
    xtrUprimeprime += p.lam6*calcCouplingTermActionContributionsImproved(6, 2, &(derivatives[0]), true);
    
    xtrU += p.lam8*calcCouplingTermActionContributionsImproved(8, 0, &(derivatives[0]), true);
    xtrUprime += p.lam8*calcCouplingTermActionContributionsImproved(8, 1, &(derivatives[0]), true);
    xtrUprimeprime += p.lam8*calcCouplingTermActionContributionsImproved(8, 2, &(derivatives[0]), true);
    
    xtrU += p.lam10*calcCouplingTermActionContributionsImproved(10, 0, &(derivatives[0]), true);
    xtrUprime += p.lam10*calcCouplingTermActionContributionsImproved(10, 1, &(derivatives[0]), true);
    xtrUprimeprime += p.lam10*calcCouplingTermActionContributionsImproved(10, 2, &(derivatives[0]), true);
    
    xtrU += 0.5*m0Sqr*(xtrTh+3*xtrTg);
    xtrUprime += 0.5*m0Sqr*(xtrThprime+3*xtrTgprime);
    xtrUprimeprime += 0.5*m0Sqr*(xtrThprimeprime+3*xtrTgprimeprime);

    U += xtrU;
    Uprime += xtrUprime;
    Uprimeprime += xtrUprimeprime;        
    
    LastPSqrParameters = p;        
  } else {
    printf("Using Integration\n");
    vector4D start;
    vector4D end;
    for (int I=0; I<4; I++) {
      start[I] = 0;
      end[I] = pi;
    }

    U = sqr(sqr(1.0/pi)) * integrate(&calcBosonicUContribution_IntegrandU, m0Sqr, start, end, INTEGRATORACCURACY);
    Uprime = sqr(sqr(1.0/pi)) * integrate(&calcBosonicUContribution_IntegrandUprime, m0Sqr, start, end, INTEGRATORACCURACY);
    Uprimeprime = sqr(sqr(1.0/pi)) * integrate(&calcBosonicUContribution_IntegrandUprimeprime, m0Sqr, start, end, INTEGRATORACCURACY);
  }  

  LastBosonicContributionParameters = p;
  LastBosonicContributionParameter_m0Sqr = m0Sqr;
  LastBosonicUContributionResult = U;
  LastBosonicUprimeContributionResult = Uprime;
  LastBosonicUprimeprimeContributionResult = Uprimeprime;
}


double calcRenormalizedQurticCoupling(ParameterType p, double m0Sqr, double mSqr) {
  double h = 1E-3;
  double v = Physical_VEV_GeV / p.CutoffInGev;
  double C0 = Physical_VEV_GeV / (v-0.5*h);
  double C1 = Physical_VEV_GeV / (v);
  double C2 = Physical_VEV_GeV / (v+0.5*h);
  
  double FU0, FUprime0, FUprimeprime0;
  double FU1, FUprime1, FUprimeprime1;
  double FU2, FUprime2, FUprimeprime2;

  p.CutoffInGev = C2;
  calcFermionUContribution(p, FU2, FUprime2, FUprimeprime2);
  p.CutoffInGev = C0;
  calcFermionUContribution(p, FU0, FUprime0, FUprimeprime0);
  p.CutoffInGev = C1;
  calcFermionUContribution(p, FU1, FUprime1, FUprimeprime1);
  
  double lamRen = (FUprimeprime2+FUprimeprime0-2*FUprimeprime1)/sqr(h);
  
  #ifdef USEBOSONICDETERMINANT
    double BU0, BUprime0, BUprimeprime0;
    double BU1, BUprime1, BUprimeprime1;
    double BU2, BUprime2, BUprimeprime2;
    
    p.CutoffInGev = C2;
    calcBosonicUContribution(p, m0Sqr, BU2, BUprime2, BUprimeprime2);
    p.CutoffInGev = C0;
    calcBosonicUContribution(p, m0Sqr, BU0, BUprime0, BUprimeprime0);
    p.CutoffInGev = C1;
    calcBosonicUContribution(p, m0Sqr, BU1, BUprime1, BUprimeprime1);
  
    lamRen += (BUprimeprime2+BUprimeprime0-2*BUprimeprime1)/sqr(h);
    lamRen += 24*p.lam0 + 360*p.lam6*v*v + 1680*p.lam8*v*v*v*v + 5040*p.lam10*v*v*v*v*v*v;
  #endif

  return lamRen/24;
}


void calcHiggsMasses(ParameterType p, double &m0Sqr, double &mSqr) {
  double FU, FUprime, FUprimeprime;
  double v = Physical_VEV_GeV / p.CutoffInGev;

  calcFermionUContribution(p, FU, FUprime, FUprimeprime);
  m0Sqr = 0.0;
  mSqr = 0;
 
  #ifdef USEBOSONICDETERMINANT
//    double BU, BUprime, BUprimeprime;
    double m0SqrOld = 0;
    int count = 0;
    int countMAX = 500;

    makeSecondOrderBosonicPotentialAvail(p);
    secondOrderBosonicPotential->setVeV(v);
    secondOrderBosonicPotential->setYukawa(p.yt0, p.yb0);    
    secondOrderBosonicPotential->setLambdas(p.lam0, p.lam6, p.lam8, p.lam10); 

    while (true) {
/*      calcBosonicUContribution(p, m0Sqr, BU, BUprime, BUprimeprime);
      m0Sqr = -4*p.lam0*v*v -6*p.lam6*v*v*v*v -8*p.lam8*v*v*v*v*v*v -10*p.lam10*v*v*v*v*v*v*v*v - FUprime/v - BUprime/v;
      mSqr = m0Sqr + 12*p.lam0*v*v +30*p.lam6*v*v*v*v + 56*p.lam8*v*v*v*v*v*v +90*p.lam10*v*v*v*v*v*v*v*v + FUprimeprime + BUprimeprime;*/
//    printf("%f %f\n",0.5*m0Sqr*v*v+p.lam0*v*v*v*v+p.lam6*v*v*v*v*v*v+p.lam8*v*v*v*v*v*v*v*v+p.lam10*v*v*v*v*v*v*v*v*v*v+BU, secondOrderBosonicPotential->calcEffectivePotDn(0));
      
      secondOrderBosonicPotential->setM0Sqr(m0Sqr);
      mSqr = FUprimeprime + secondOrderBosonicPotential->calcEffectivePotDn(2);
      m0Sqr = -(FUprime + secondOrderBosonicPotential->calcEffectivePotDn(1) - m0Sqr*v)/v;
      
      if ((fabs((m0Sqr-m0SqrOld)/m0Sqr))<1E-10) break;
      m0SqrOld = m0Sqr;
      count++;
      if (count>countMAX) {
        m0Sqr = NaN;
        mSqr = NaN;
	return;	
      }
    }
  #else
/*    double S1 = NaN;
    double S2 = NaN;
    double dummy = NaN;
  
    calcPropagatorSums(p, 0, S1, dummy);
    S2 = calcSquaredHiggsPropagatorSum(p, 0);
    
    double fac = 1 + 48*S1*S1*p.lam0 + 12*p.lam0*S2;
  
  
    m0Sqr = (-1.0/fac)*(4*p.lam0*v*v + 24*p.lam0*S1 - 288*p.lam0*p.lam0*v*v*S1*S1 - 96*p.lam0*p.lam0*v*v*S2 + FUprime/v);
    mSqr  = m0Sqr + 12*p.lam0*v*v + 24*p.lam0*S1 + S1*S1*(48*p.lam0*m0Sqr - 864*p.lam0*p.lam0*v*v)
          + S2*(12*p.lam0*m0Sqr - 288*p.lam0*p.lam0*v*v) + FUprimeprime;*/
  
    double SelfCouplingCorrection = 0;
    double SelfCouplingV4Correction = 0;    
    double hPropSum = 0;
    double gPropSum = 0;
    double oldCorr = NaN;   //Wenn dieser Term konvergiert ist, so auch alle einzelnen Propagator-Summen
    double oldV4Corr = NaN;
    int count = 0;
    int countMAX = 500;
        
    while (true) {
      double lam6prime = calcCouplingTermActionContributions(6, 1, v, hPropSum, gPropSum);
      double lam6primeprime = calcCouplingTermActionContributions(6, 2, v, hPropSum, gPropSum);
      double lam8prime = calcCouplingTermActionContributions(8, 1, v, hPropSum, gPropSum);
      double lam8primeprime = calcCouplingTermActionContributions(8, 2, v, hPropSum, gPropSum);
      double lam10prime = calcCouplingTermActionContributions(10, 1, v, hPropSum, gPropSum);
      double lam10primeprime = calcCouplingTermActionContributions(10, 2, v, hPropSum, gPropSum);
    
      m0Sqr = -4*(p.lam0-p.lam0*p.lam0*SelfCouplingV4Correction)*v*v - FUprime/v 
            - 2*p.lam0*SelfCouplingCorrection
	    - p.lam6*(lam6prime/v) - p.lam8*(lam8prime/v) - p.lam10*(lam10prime/v);

			//Discrepancy with the factor here. Phillips thesis says 8 rather than 12. 
      mSqr = m0Sqr + 12*(p.lam0-p.lam0*p.lam0*SelfCouplingV4Correction)*v*v 
           + FUprimeprime + 2*p.lam0*SelfCouplingCorrection
	   + p.lam6*(lam6primeprime) + p.lam8*(lam8primeprime) + p.lam10*(lam10primeprime);
   
      #ifdef USEIMPROVEMENT
        SelfCouplingCorrection = calcSelfCouplingCorrectionCoefficientForEffectivePotential(p, mSqr);      
        SelfCouplingV4Correction = calcSelfCouplingV4CorrectionCoefficientForEffectivePotential(p, mSqr);
	calcPropagatorSums(p, mSqr, hPropSum, gPropSum);	
        if ((!isNaN(oldCorr)) && ((fabs(oldCorr-SelfCouplingCorrection)/SelfCouplingCorrection)<1E-10) && (!isNaN(oldV4Corr)) && ((fabs(oldV4Corr-SelfCouplingV4Correction)/SelfCouplingV4Correction)<1E-10)) break;
      #else
        break;
      #endif
      oldCorr = SelfCouplingCorrection;
      oldV4Corr = SelfCouplingV4Correction;
      count++;
      if (count>countMAX) {
        m0Sqr = NaN;
        mSqr = NaN;
	return;	
      }
    }  
  #endif
}


double FindLocalMinimum_Helper_m0Sqr = NaN;
double FindLocalMinimum_Helper_mSqr = NaN;
double FindLocalMinimum_Helper(double* x) {
  double FU, FUprime, FUprimeprime;
  double CutoffSaver = Parameters.CutoffInGev;
  double v = *x;
  Parameters.CutoffInGev = Physical_VEV_GeV / v;
  double S = NaN;
  calcFermionUContribution(Parameters, FU, FUprime, FUprimeprime);

  #ifdef USEBOSONICDETERMINANT
    double BU, BUprime, BUprimeprime;
    calcBosonicUContribution(Parameters, FindLocalMinimum_Helper_m0Sqr, BU, BUprime, BUprimeprime);
    S = FU + BU + 0.5*FindLocalMinimum_Helper_m0Sqr*v*v + Parameters.lam0*v*v*v*v + Parameters.lam6*v*v*v*v*v*v + Parameters.lam8*v*v*v*v*v*v*v*v + Parameters.lam10*v*v*v*v*v*v*v*v*v*v;
  #else
    #ifdef USEIMPROVEMENT
      double hPropSum = NaN;
      double gPropSum = NaN;
      calcPropagatorSums(Parameters, FindLocalMinimum_Helper_mSqr, hPropSum, gPropSum);
      double SelfCouplingCorrection = calcSelfCouplingCorrectionCoefficientForEffectivePotential(Parameters, FindLocalMinimum_Helper_mSqr);      
      double SelfCouplingV4Correction = calcSelfCouplingV4CorrectionCoefficientForEffectivePotential(Parameters, FindLocalMinimum_Helper_mSqr);

      double lam6f = calcCouplingTermActionContributions(6, 0, v, hPropSum, gPropSum);
      double lam8f = calcCouplingTermActionContributions(8, 0, v, hPropSum, gPropSum);
      double lam10f = calcCouplingTermActionContributions(10, 0, v, hPropSum, gPropSum);
    
      S = FU + 0.5*FindLocalMinimum_Helper_m0Sqr*v*v + Parameters.lam0*(v*v*v*v + SelfCouplingCorrection*v*v)
        - Parameters.lam0*Parameters.lam0*SelfCouplingV4Correction*v*v*v*v
        + Parameters.lam6*(lam6f)
        + Parameters.lam8*(lam8f) + Parameters.lam10*(lam10f);
    #else
      S = FU + 0.5*FindLocalMinimum_Helper_m0Sqr*v*v + Parameters.lam0*v*v*v*v + Parameters.lam6*v*v*v*v*v*v + Parameters.lam8*v*v*v*v*v*v*v*v + Parameters.lam10*v*v*v*v*v*v*v*v*v*v;
    #endif
  #endif      
  
  Parameters.CutoffInGev = CutoffSaver;
  return S;
}


double FindAbsoluteMinimum(double m0Sqr, double mSqr, double maxV, int tries) {  
  double bestS = NaN;
  double bestV = NaN;
  for (int I=0; I<=tries; I++) {
    double pos = I*(maxV/tries);
    double bounds1 = 0;
    double bounds2 = 5;
    FindLocalMinimum_Helper_m0Sqr = m0Sqr;
    FindLocalMinimum_Helper_mSqr = mSqr;
    GradientMinimization(&FindLocalMinimum_Helper, 1, 1E-2, 1E-4, 1E-6, &pos, &bounds1, &bounds2, NULL, 1, 1000);
    double S = FindLocalMinimum_Helper(&pos);
    if ((!isNaN(pos)) && (!isNaN(S))) {
      if ((isNaN(bestS)) || (isNaN(bestV))) {
        bestS = S;
	bestV = pos;
      }
      if (S<bestS) {
        bestS = S;
	bestV = pos;
      }
    }
  }
  return bestV;
}


void plotActionDependenceOnV() {
  printf("Plotting action to file...\n");
  double FU, FUprime, FUprimeprime;
  FILE* file = fopen("AnalyticActionPlot.dat","w");

  double m0Sqr, mSqr;
  calcHiggsMasses(Parameters, m0Sqr, mSqr);
  
  printf("... with m0Sqr:%f, mSqr:%f\n",m0Sqr, mSqr);
  double CutoffSaver = Parameters.CutoffInGev;
    
  for (double fac=0.3; fac<13.3; fac*=1.01) {
    Parameters.CutoffInGev = CutoffSaver*fac;
    double v = Physical_VEV_GeV / Parameters.CutoffInGev;
    double S = NaN;
    double dS = NaN;
    double ddS = NaN;
      
    calcFermionUContribution(Parameters, FU, FUprime, FUprimeprime);
    #ifdef USEBOSONICDETERMINANT
      double BU, BUprime, BUprimeprime;
      calcBosonicUContribution(Parameters, m0Sqr, BU, BUprime, BUprimeprime);
      
      S = FU + BU + 0.5*m0Sqr*v*v + Parameters.lam0*v*v*v*v + Parameters.lam6*v*v*v*v*v*v + Parameters.lam8*v*v*v*v*v*v*v*v + Parameters.lam10*v*v*v*v*v*v*v*v*v*v;
      dS = FUprime + BUprime + m0Sqr*v + 4*Parameters.lam0*v*v*v + 6*Parameters.lam6*v*v*v*v*v + 8*Parameters.lam8*v*v*v*v*v*v*v + 10*Parameters.lam10*v*v*v*v*v*v*v*v*v;
      ddS = FUprimeprime + BUprimeprime + m0Sqr + 12*Parameters.lam0*v*v + 24*Parameters.lam6*v*v*v*v + 56*Parameters.lam8*v*v*v*v*v*v + 90*Parameters.lam10*v*v*v*v*v*v*v*v;
    #else
      #ifdef USEIMPROVEMENT
        double hPropSum = NaN;
        double gPropSum = NaN;
        calcPropagatorSums(Parameters, mSqr, hPropSum, gPropSum);
        double SelfCouplingCorrection = calcSelfCouplingCorrectionCoefficientForEffectivePotential(Parameters, mSqr);      
        double SelfCouplingV4Correction = calcSelfCouplingV4CorrectionCoefficientForEffectivePotential(Parameters, mSqr);

        double lam6f = calcCouplingTermActionContributions(6, 0, v, hPropSum, gPropSum);
        double lam6prime = calcCouplingTermActionContributions(6, 1, v, hPropSum, gPropSum);
        double lam6primeprime = calcCouplingTermActionContributions(6, 2, v, hPropSum, gPropSum);
        double lam8f = calcCouplingTermActionContributions(8, 0, v, hPropSum, gPropSum);
        double lam8prime = calcCouplingTermActionContributions(8, 1, v, hPropSum, gPropSum);
        double lam8primeprime = calcCouplingTermActionContributions(8, 2, v, hPropSum, gPropSum);
        double lam10f = calcCouplingTermActionContributions(10, 0, v, hPropSum, gPropSum);
        double lam10prime = calcCouplingTermActionContributions(10, 1, v, hPropSum, gPropSum);
        double lam10primeprime = calcCouplingTermActionContributions(10, 2, v, hPropSum, gPropSum);
    
        S = FU + 0.5*m0Sqr*v*v + Parameters.lam0*(v*v*v*v + SelfCouplingCorrection*v*v)
          - Parameters.lam0*Parameters.lam0*SelfCouplingV4Correction*v*v*v*v
          + Parameters.lam6*(lam6f)
          + Parameters.lam8*(lam8f) + Parameters.lam10*(lam10f);
   
        dS = FUprime + m0Sqr*v + Parameters.lam0*(4*v*v*v + 2*SelfCouplingCorrection*v)
            - 4*Parameters.lam0*Parameters.lam0*SelfCouplingV4Correction*v*v*v	
	    + Parameters.lam6*(lam6prime)
   	    + Parameters.lam8*(lam8prime) + Parameters.lam10*(lam10prime);
	ddS = FUprimeprime + m0Sqr + Parameters.lam0*(12*v*v + 2*SelfCouplingCorrection)	
            - 12*Parameters.lam0*Parameters.lam0*SelfCouplingV4Correction*v*v	
	    + Parameters.lam6*(lam6primeprime)
            + Parameters.lam8*(lam8primeprime) + Parameters.lam10*(lam10primeprime);

      #else
        S = FU + 0.5*m0Sqr*v*v + Parameters.lam0*v*v*v*v + Parameters.lam6*v*v*v*v*v*v + Parameters.lam8*v*v*v*v*v*v*v*v + Parameters.lam10*v*v*v*v*v*v*v*v*v*v;
        dS = FUprime + m0Sqr*v + 4*Parameters.lam0*v*v*v + 6*Parameters.lam6*v*v*v*v*v + 8*Parameters.lam8*v*v*v*v*v*v*v + 10*Parameters.lam10*v*v*v*v*v*v*v*v*v;
        ddS = FUprimeprime + m0Sqr + 12*Parameters.lam0*v*v + 24*Parameters.lam6*v*v*v*v + 56*Parameters.lam8*v*v*v*v*v*v + 90*Parameters.lam10*v*v*v*v*v*v*v*v;
      #endif
    #endif      
    fprintf(file, "%1.15f %1.15f %1.15f %1.15f\n",v,S,dS,ddS);        
  }
  Parameters.CutoffInGev = CutoffSaver;
  fclose(file);
}


void randomScanParameterSpace(double targetMassGeV,double targetL6, double targetL8, int fixLevel) {
  double FU, FUprime, FUprimeprime;
  double v = Physical_VEV_GeV / Parameters.CutoffInGev;
  
  #ifndef USEBOSONICDETERMINANT
    #ifdef USEIMPROVEMENT
      double FU2, FUprime2, FUprimeprime2;
      double CutoffSaver = Parameters.CutoffInGev;
      double h = 1E-5;
      Parameters.CutoffInGev = Physical_VEV_GeV/(v+0.5*h);
      calcFermionUContribution(Parameters, FU2, FUprime2, FUprimeprime2);
      Parameters.CutoffInGev = Physical_VEV_GeV/(v-0.5*h);
      calcFermionUContribution(Parameters, FU, FUprime, FUprimeprime);
      double FUprimeprimeprime = (FUprimeprime2-FUprimeprime) / h;  
      Parameters.CutoffInGev = CutoffSaver;
    #endif
  #endif
    
  double smallestMass = NaN;
  calcFermionUContribution(Parameters, FU, FUprime, FUprimeprime);  

//  while (true) {
    Parameters.lam0 = (AdvancedZufall(AdvancedSeed)-0.5) * 0.001;
    Parameters.lam6 = (AdvancedZufall(AdvancedSeed)-0.0) * 0.001;
    Parameters.lam8 = (AdvancedZufall(AdvancedSeed)-0.0) * 0.00;
    Parameters.lam10 = (AdvancedZufall(AdvancedSeed)) * 0.000;
    
Parameters.lam6 = targetL6;
Parameters.lam8 = targetL8;
    
    
    double m0Sqr = 0.0;
    double mSqr = 0;    
    calcHiggsMasses(Parameters, m0Sqr, mSqr);    
    
    #ifndef USEBOSONICDETERMINANT
      #ifdef USEIMPROVEMENT
        int FixationLevel = fixLevel;
	if (FixationLevel>0) {
	  double targetMSqr = sqr(targetMassGeV/Parameters.CutoffInGev);
          double hPropSum = 0;
          double gPropSum = 0;
          double SelfCouplingCorrection = calcSelfCouplingCorrectionCoefficientForEffectivePotential(Parameters, targetMSqr);      
   	  calcPropagatorSums(Parameters, targetMSqr, hPropSum, gPropSum);

          double lam6prime = calcCouplingTermActionContributions(6, 1, v, hPropSum, gPropSum);
          double lam6primeprime = calcCouplingTermActionContributions(6, 2, v, hPropSum, gPropSum);
          double lam6primeprimeprime = calcCouplingTermActionContributions(6, 3, v, hPropSum, gPropSum);	
          double lam8prime = calcCouplingTermActionContributions(8, 1, v, hPropSum, gPropSum);
          double lam8primeprime = calcCouplingTermActionContributions(8, 2, v, hPropSum, gPropSum);
          double lam8primeprimeprime = calcCouplingTermActionContributions(8, 3, v, hPropSum, gPropSum);
          double lam10prime = calcCouplingTermActionContributions(10, 1, v, hPropSum, gPropSum);
          double lam10primeprime = calcCouplingTermActionContributions(10, 2, v, hPropSum, gPropSum);
          double lam10primeprimeprime = calcCouplingTermActionContributions(10, 3, v, hPropSum, gPropSum);
    
          ComplexVector vec(1+FixationLevel);
  	  vec.setZero();
	  ComplexMatrix mat(1+FixationLevel);
  	  mat.setZero();
	  if (FixationLevel==1) {
            vec.vectorElements[0].x = -FUprime/v - Parameters.lam6*(lam6prime/v) - Parameters.lam8*(lam8prime/v) - Parameters.lam10*(lam10prime/v);
            vec.vectorElements[1].x = +targetMSqr -FUprimeprime - Parameters.lam6*(lam6primeprime) - Parameters.lam8*(lam8primeprime) - Parameters.lam10*(lam10primeprime);

            mat.matrix[0][0].x = 1;
            mat.matrix[0][1].x = 4*v*v + 2*SelfCouplingCorrection;
            mat.matrix[1][0].x = 1;
            mat.matrix[1][1].x = 12*v*v + 2*SelfCouplingCorrection; 
	  }
  	  if (FixationLevel==2) {
            vec.vectorElements[0].x = -FUprime/v - Parameters.lam8*(lam8prime/v) - Parameters.lam10*(lam10prime/v);
            vec.vectorElements[1].x = +targetMSqr -FUprimeprime - Parameters.lam8*(lam8primeprime) - Parameters.lam10*(lam10primeprime);
            vec.vectorElements[2].x = -FUprimeprimeprime  - Parameters.lam8*(lam8primeprimeprime) - Parameters.lam10*(lam10primeprimeprime);

            mat.matrix[0][0].x = 1;
            mat.matrix[0][1].x = 4*v*v + 2*SelfCouplingCorrection;
            mat.matrix[0][2].x = lam6prime/v;
	 
            mat.matrix[1][0].x = 1;
            mat.matrix[1][1].x = 12*v*v + 2*SelfCouplingCorrection;
            mat.matrix[1][2].x = lam6primeprime;
	 
            mat.matrix[2][0].x = 0;
            mat.matrix[2][1].x = 24*v;
            mat.matrix[2][2].x = lam6primeprimeprime;
  	  }
	
  	  mat.invert();
	
 	  ComplexVector res = mat*vec;
	  m0Sqr = res.vectorElements[0].x;
	  Parameters.lam0 = res.vectorElements[1].x;
	  if (FixationLevel>=2) Parameters.lam6 = res.vectorElements[2].x;
	

          mSqr = m0Sqr + 12*(Parameters.lam0)*v*v 
               + FUprimeprime + 2*Parameters.lam0*SelfCouplingCorrection
	       + Parameters.lam6*(lam6primeprime) + Parameters.lam8*(lam8primeprime) + Parameters.lam10*(lam10primeprime);
        }   
      #endif
    #endif
    printf("\nSelected point is v=%f, mSqr=%f\n",v,mSqr);
    printf("  lam0:  %1.5e\n",Parameters.lam0);
    printf("  lam6:  %1.5e\n",Parameters.lam6);
    printf("  lam8:  %1.5e\n",Parameters.lam8);
    printf("  lam10: %1.5e\n",Parameters.lam10);
    
    if (((isNaN(smallestMass)) || (mSqr<smallestMass)) && (mSqr>=0) && (!isNaN(mSqr))) {
      double realMin = FindAbsoluteMinimum(m0Sqr, mSqr, 2, 5);
      double CutoffSaver = Parameters.CutoffInGev;
      Parameters.CutoffInGev = Physical_VEV_GeV / realMin;
      double FU2, FUprime2, FUprimeprime2;      
      calcFermionUContribution(Parameters, FU2, FUprime2, FUprimeprime2);

      #ifdef USEBOSONICDETERMINANT
        double BU, BUprime, BUprimeprime;
        calcBosonicUContribution(Parameters, FindLocalMinimum_Helper_m0Sqr, BU, BUprime, BUprimeprime);
        mSqr = m0Sqr + 12*Parameters.lam0*realMin*realMin +30*Parameters.lam6*realMin*realMin*realMin*realMin + 56*Parameters.lam8*realMin*realMin*realMin*realMin*realMin*realMin +90*Parameters.lam10*realMin*realMin*realMin*realMin*realMin*realMin*realMin*realMin + FUprimeprime2 + BUprimeprime;
      #else
        #ifdef USEIMPROVEMENT
          double SelfCouplingCorrection = 0;
          double SelfCouplingV4Correction = 0;
          double hPropSum = 0;
          double gPropSum = 0;
          double oldCorr = NaN;   //Wenn dieser Term konvergiert ist, so auch alle einzelnen Propagator-Summen
          double oldV4Corr = NaN;
	  int count = 0;
	  int countMAX=20;

          while (true) {
            double lam6primeprime = calcCouplingTermActionContributions(6, 2, realMin, hPropSum, gPropSum);
            double lam8primeprime = calcCouplingTermActionContributions(8, 2, realMin, hPropSum, gPropSum);
            double lam10primeprime = calcCouplingTermActionContributions(10, 2, realMin, hPropSum, gPropSum);

            mSqr = m0Sqr + 12*(Parameters.lam0-Parameters.lam0*Parameters.lam0*SelfCouplingV4Correction)*realMin*realMin 
                 + FUprimeprime + 2*Parameters.lam0*SelfCouplingCorrection
	         + Parameters.lam6*(lam6primeprime) + Parameters.lam8*(lam8primeprime) + Parameters.lam10*(lam10primeprime);

            SelfCouplingCorrection = calcSelfCouplingCorrectionCoefficientForEffectivePotential(Parameters, mSqr);      
            SelfCouplingV4Correction = calcSelfCouplingV4CorrectionCoefficientForEffectivePotential(Parameters, mSqr);
   	    calcPropagatorSums(Parameters, mSqr, hPropSum, gPropSum);
            if ((!isNaN(oldCorr)) && ((fabs(oldCorr-SelfCouplingCorrection)/SelfCouplingCorrection)<1E-10) && (!isNaN(oldV4Corr)) && ((fabs(oldV4Corr-SelfCouplingV4Correction)/SelfCouplingV4Correction)<1E-10)) break;
            oldCorr = SelfCouplingCorrection;
            oldV4Corr = SelfCouplingV4Correction;
	    count++;
	    
	    if (count>countMAX) {
	      mSqr = NaN;
	      break;
	    }
          }  
        #else
          mSqr = m0Sqr + 12*Parameters.lam0*realMin*realMin +30*Parameters.lam6*realMin*realMin*realMin*realMin + 56*Parameters.lam8*realMin*realMin*realMin*realMin*realMin*realMin +90*Parameters.lam10*realMin*realMin*realMin*realMin*realMin*realMin*realMin*realMin + FUprimeprime2;
        #endif
      #endif
      printf("\nReal Minimum is v=%f, mSqr=%f\n",realMin,mSqr);

      if ((!isNaN(mSqr)) && (mSqr>=0)) {
        if ((!isNaN(realMin)) && ((fabs(realMin-v)/v)<0.05)) {
          if (isNaN(smallestMass)) smallestMass = mSqr;
          if (mSqr<smallestMass) smallestMass = mSqr;  
        }
      }
      double cutoff = Parameters.CutoffInGev;
      Parameters.CutoffInGev = CutoffSaver;  

      if ((!isNaN(mSqr)) && (mSqr>=0)) {
        char* fileName = constructFileNameShort(Parameters, "RandomParameterSpaceScan");      
        FILE* file = fopen(fileName,"a"); 
        if ((!isNaN(cutoff)) && (!isNaN(sqrt(mSqr)*Parameters.CutoffInGev))) {
          fprintf(file, "%1.15f %1.15f %1.15f %1.15f %1.15f %1.15f %1.15f %1.15f %1.15f \n", cutoff, sqrt(mSqr)*Parameters.CutoffInGev, realMin, mSqr, m0Sqr, Parameters.lam0, Parameters.lam6, Parameters.lam8, Parameters.lam10);
        }
        
        fclose(file);
        delete[] fileName;
      }
    }
//  }
}


void calcEffectiveFermionMassesFromMomentumPropagator(int Lt, double* tresProp, double* &effMasses, double* &corr) {
  corr = new double[Lt+1];
  effMasses = new double[Lt];
  for (int I=0; I<Lt+1; I++) corr[I] = 0;
  
  Complex alpha(0,2.0*pi/Lt);
  
  for (int dt=0; dt<=Lt; dt++) {
    double DeltaT = dt;
    Complex res(0,0);
    for (int k=0; k<Lt; k++) {
      res = res + tresProp[k] * exp(k*DeltaT * alpha);
    }    
    
    corr[Lt-dt] += 0.5 * res.x/Lt;
    corr[dt] += 0.5 * res.x/Lt;
  }

  for (int dt=0; dt<Lt; dt++) {
    if (dt<Lt/2) {
      effMasses[dt] = EffectiveMassSolver(dt, dt+1, corr[dt], corr[dt+1], Lt);
    } else {
      effMasses[dt] = EffectiveMassSolver(dt+1, dt, corr[dt+1], corr[dt], Lt);
    }
  }
}


bool areParametersEqualWithRespectToDiracOpEValsandEVecs(ParameterType p1, ParameterType p2) {
  if (p1.L0 != p2.L0) return false;
  if (p1.L1 != p2.L1) return false;
  if (p1.L2 != p2.L2) return false;
  if (p1.L3 != p2.L3) return false;
  if (p1.rho != p2.rho) return false;
  if (p1.r != p2.r) return false;

  return true;
}


void makeDiracOpEigenvaluesAndEigenvectorsAvailable(ParameterType p) {
  if ((DiracOpEigenvalues == NULL) || (DiracOpEigenvectors == NULL) || (!areParametersEqualWithRespectToDiracOpEValsandEVecs(p, LastDiracOperatorEigenvaluesAndEigenvectorsParameters))) {
    printf("Generating eigenvalues and eigenvectors of Dirac-Operator...\n");
    delete[] DiracOpEigenvalues;
    if (DiracOpEigenvectors!=NULL) {
      for (int I=0; I<p.L0*p.L1*p.L2*p.L3; I++) {
        delete[] DiracOpEigenvectors[I];
      }
    }
    delete[] DiracOpEigenvectors;
    
    DiracOpEigenvalues = new Complex[p.L0*p.L1*p.L2*p.L3];
    DiracOpEigenvectors = new ComplexVector*[p.L0*p.L1*p.L2*p.L3];

    vector4D k;    
    int count = 0;
    for (int i0=0; i0<p.L0; i0++) {
      k[0] = 2*pi*i0 / p.L0;
      for (int i1=0; i1<p.L1; i1++) {
        k[1] = 2*pi*i1 / p.L1;
        for (int i2=0; i2<p.L2; i2++) {
          k[2] = 2*pi*i2 / p.L2;
          for (int i3=0; i3<p.L3; i3++) {
            k[3] = 2*pi*i3 / p.L3;

            DiracOpEigenvalues[count] = diracOp->analyticalEigenvalue(k);
	    DiracOpEigenvectors[count] = new ComplexVector[4];
	    DiracOpEigenvectors[count][0].resize(4);
	    DiracOpEigenvectors[count][1].resize(4);
	    DiracOpEigenvectors[count][2].resize(4);
	    DiracOpEigenvectors[count][3].resize(4);
	    
            diracOp->analyticalEigenvectors(k, DiracOpEigenvectors[count]);

            count++;	      
	  }
	}
      }
    }
    LastDiracOperatorEigenvaluesAndEigenvectorsParameters = p;
  }
}


ComplexMatrix calcPropagatorMatrixSandwhichProduct(ComplexMatrix prop, ComplexMatrix T) {
  int mSize = prop.matrixSize;
  ComplexMatrix res(mSize);
  res.setZero();
  
  for (int alpha=0; alpha<mSize; alpha++) {
    for (int beta=0; beta<mSize; beta++) {
  
      for (int gamma=0; gamma<mSize; gamma++) {
        for (int delta=0; delta<mSize; delta++) {
          if ((alpha != delta) && (beta != gamma)) {
            res.matrix[alpha][beta] = res.matrix[alpha][beta] + (prop.matrix[alpha][gamma] * T.matrix[gamma][delta] * prop.matrix[delta][beta]);
	  }
	}
      }
    }
  }
  return res;
}


void calcFermionMasses(ParameterType p, double m0Sqr, double mHiggsSqr, double &topMass, double &bottomMass) {
  topMass = 0;
  bottomMass = 0;  
  if ((p.yt0==0) && (p.yb0==0)) return;
  printf("Calculating Fermion Masses\n");

  if ((p.L0>=0) && (p.L1>=0) && (p.L2>=0) && (p.L3>=0)) {
    makePSqrAvailableForParameters(p);
    makeDiracOpEigenvaluesAndEigenvectorsAvailable(p);
    
    double yt = Parameters.yt0;
    double yb = Parameters.yb0;  
    double v = Physical_VEV_GeV / Parameters.CutoffInGev;
    double fac = 0.5 / Parameters.rho;

    int PSqrIndex = 0;
    vector4D k;
    ComplexVector vK[4];
    vK[0].resize(4);
    vK[1].resize(4);
    vK[2].resize(4);
    vK[3].resize(4);
    ComplexVector vP[4];
    vP[0].resize(4);
    vP[1].resize(4);
    vP[2].resize(4);
    vP[3].resize(4);
    
    ComplexMatrix ProjPlusSmall = getProjectorMatrix(+1);
    ComplexMatrix ProjMinusSmall = getProjectorMatrix(-1);    
    ComplexMatrix ThetaSmall[4];
    ComplexMatrix BMat[4];
    ComplexMatrix diagY(2);
    diagY.setZero();
    diagY.matrix[0][0].x = yt;
    diagY.matrix[1][1].x = yb;
    for (int I=0; I<4; I++) {
      ThetaSmall[I] = getThetaMatrix(I);
      BMat[I] = ComplexMatrix(diagY*ThetaSmall[I], ProjMinusSmall);
      ThetaSmall[I].dagger();
      BMat[I] = BMat[I] + ComplexMatrix(ThetaSmall[I]*diagY, ProjPlusSmall);      
      ThetaSmall[I].dagger();
    }

    double* topPropTres = new double[p.L3];
    double* bottomPropTres = new double[p.L3];
    
    for (int outerMomentumP3=0; outerMomentumP3<=p.L3/2; outerMomentumP3++) {
      calcFermionMassesExternalMomentumP[0] = 0;
      calcFermionMassesExternalMomentumP[1] = 0;
      calcFermionMassesExternalMomentumP[2] = 0;
      calcFermionMassesExternalMomentumP[3] = 2*pi*outerMomentumP3 / p.L3;

      Complex ewP = diracOp->analyticalEigenvalue(calcFermionMassesExternalMomentumP);
      diracOp->analyticalEigenvectors(calcFermionMassesExternalMomentumP, vP);

      ComplexMatrix mP1(vP[0]);  
      ComplexMatrix mP2(vP[1]);  
      ComplexMatrix mP3(vP[2]);  
      ComplexMatrix mP4(vP[3]);  
      Complex ewGamP = (ComplexUnity-fac*ewP);
      ComplexMatrix GamPSmall = ewGamP*mP1 + ewGamP*mP2 + adj(ewGamP)*mP3 + adj(ewGamP)*mP4;
      ComplexMatrix GamP(8);
      GamP.setZero();
      GamP.insertMatrix(GamPSmall, 0, 0);
      GamP.insertMatrix(GamPSmall, 4, 4);
	      
      Complex ewDinvPtop    = ComplexUnity / (ewP + yt*v*(ComplexUnity-fac*ewP));
      Complex ewDinvPbottom = ComplexUnity / (ewP + yb*v*(ComplexUnity-fac*ewP));
      ComplexMatrix DinvPtop    = ewDinvPtop*mP1 + ewDinvPtop*mP2 + adj(ewDinvPtop)*mP3 + adj(ewDinvPtop)*mP4;
      ComplexMatrix DinvPbottom = ewDinvPbottom*mP1 + ewDinvPbottom*mP2 + adj(ewDinvPbottom)*mP3 + adj(ewDinvPbottom)*mP4;
      ComplexMatrix DinvP(8);
      DinvP.setZero();
      DinvP.insertMatrix(DinvPtop, 0, 0);
      DinvP.insertMatrix(DinvPbottom, 4, 4);
      
      ComplexMatrix DP = DinvP;
      DP.invert();
 
      ComplexMatrix HiggsContribution(8);
      HiggsContribution.setZero();
      ComplexMatrix Goldstone1Contribution(8);
      Goldstone1Contribution.setZero();
      ComplexMatrix Goldstone2Contribution(8);
      Goldstone2Contribution.setZero();
      ComplexMatrix Goldstone3Contribution(8);
      Goldstone3Contribution.setZero();
      ComplexMatrix TadPoleContribution(8);
      TadPoleContribution.setZero();
      ComplexMatrix CounterTadPoleContribution(8);
      CounterTadPoleContribution.setZero();

      int count = 0;
      for (int i0=0; i0<p.L0; i0++) {
        k[0] = 2*pi*i0 / p.L0;
        for (int i1=0; i1<p.L1; i1++) {
          k[1] = 2*pi*i1 / p.L1;
          for (int i2=0; i2<p.L2; i2++) {
            k[2] = 2*pi*i2 / p.L2;
            for (int i3=0; i3<p.L3; i3++) {
              k[3] = 2*pi*i3 / p.L3;
	      PSqrIndex = i3-outerMomentumP3;
	      if (PSqrIndex<0) PSqrIndex = -PSqrIndex;
  	      PSqrIndex = PSqrIndex % p.L3;
	      PSqrIndex = PSqrIndex + i2*p.L3 + i1*p.L2*p.L3 + i0*p.L1*p.L2*p.L3;
              double HiggsProp = 1.0 / (pSqr[PSqrIndex] + mHiggsSqr);
              double GoldstoneProp = 1.0 / (pSqr[PSqrIndex] + 0);    

//              Complex ewK = diracOp->analyticalEigenvalue(k);
//              diracOp->analyticalEigenvectors(k, vK);
              Complex ewK = DiracOpEigenvalues[count];
	      for (int I=0; I<4; I++) {
  	        vK[I] = DiracOpEigenvectors[count][I];
              }
	      
              ComplexMatrix mK1(vK[0]);  
              ComplexMatrix mK2(vK[1]);  
              ComplexMatrix mK3(vK[2]);  
              ComplexMatrix mK4(vK[3]);  
              Complex ewGamDinvKtop    = (ComplexUnity-fac*ewK) / (ewK + yt*v*(ComplexUnity-fac*ewK));
              Complex ewGamDinvKbottom = (ComplexUnity-fac*ewK) / (ewK + yb*v*(ComplexUnity-fac*ewK));
              ComplexMatrix GamDinvKtop    = ewGamDinvKtop*mK1 + ewGamDinvKtop*mK2 + adj(ewGamDinvKtop)*mK3 + adj(ewGamDinvKtop)*mK4;
              ComplexMatrix GamDinvKbottom = ewGamDinvKbottom*mK1 + ewGamDinvKbottom*mK2 + adj(ewGamDinvKbottom)*mK3 + adj(ewGamDinvKbottom)*mK4;
              ComplexMatrix GamDinvK(8);
              GamDinvK.setZero();
              GamDinvK.insertMatrix(GamDinvKtop, 0, 0);
              GamDinvK.insertMatrix(GamDinvKbottom, 4, 4);
	      
              HiggsContribution = HiggsContribution + HiggsProp * BMat[0] * GamDinvK * BMat[0] * GamP;
	      if (PSqrIndex!=0) {
 	        Goldstone1Contribution = Goldstone1Contribution + GoldstoneProp * BMat[1] * GamDinvK * BMat[1] * GamP;
 	        Goldstone2Contribution = Goldstone2Contribution + GoldstoneProp * BMat[2] * GamDinvK * BMat[2] * GamP;
 	        Goldstone3Contribution = Goldstone3Contribution + GoldstoneProp * BMat[3] * GamDinvK * BMat[3] * GamP;
	      }
	      TadPoleContribution = TadPoleContribution + (-1.0/mHiggsSqr) * (GamDinvK * (BMat[0])).tres() * (BMat[0]) * GamP;
	      CounterTadPoleContribution = CounterTadPoleContribution + (1.0/mHiggsSqr) * m0Sqr*v * (BMat[0]) * GamP;
	      
	      count++;
            }
	  }
        }      
      }
      
      ComplexMatrix TotalContribution = (1.0/(p.L0*p.L1*p.L2*p.L3))* (HiggsContribution + Goldstone1Contribution + Goldstone2Contribution + Goldstone3Contribution /* + TadPoleContribution + CounterTadPoleContribution*/);
 
      ComplexMatrix res = DP + (-1.0) * TotalContribution;
      ComplexMatrix topProp = res.getSubMatrix(0,0, 4);
      ComplexMatrix bottomProp = res.getSubMatrix(4,4, 4);
      topProp.invert();
      bottomProp.invert();
      
      topPropTres[outerMomentumP3] = topProp.tres().x;
      bottomPropTres[outerMomentumP3] = bottomProp.tres().x;
    }
    for (int outerMomentumP3=1+p.L3/2; outerMomentumP3<p.L3; outerMomentumP3++) {
      topPropTres[outerMomentumP3] = topPropTres[p.L3-outerMomentumP3];
      bottomPropTres[outerMomentumP3] = bottomPropTres[p.L3-outerMomentumP3];
    }

    double* topEffMasses = NULL;
    double* topCorr = NULL;
    double* bottomEffMasses = NULL;
    double* bottomCorr = NULL;
    calcEffectiveFermionMassesFromMomentumPropagator(p.L3, topPropTres, topEffMasses, topCorr);
    calcEffectiveFermionMassesFromMomentumPropagator(p.L3, bottomPropTres, bottomEffMasses, bottomCorr);
    topMass = topEffMasses[p.L3/2-1];
    bottomMass = bottomEffMasses[p.L3/2-1];

/*printf("TopMass: %f\n", topMass);    
printf("BottomMass: %f\n", bottomMass);    
    FILE* file = fopen("FermionCorrelator.dat", "w");
    for (int I=0; I<p.L3+1; I++) {
      fprintf(file, "%e %e %e\n", (double)I, topCorr[I], bottomCorr[I]);
    }
    fclose(file);
    file = fopen("FermionEffectiveMasses.dat", "w");
    for (int I=0; I<p.L3; I++) {
      fprintf(file, "%e %e %e %e %e\n", I+0.5, topEffMasses[I], yt*v, bottomEffMasses[I], yb*v);
    }
    fclose(file);
exit(0); */
 
    delete[] topEffMasses;
    delete[] topCorr;
    delete[] bottomEffMasses;
    delete[] bottomCorr;
    delete[] topPropTres;
    delete[] bottomPropTres;
  }
}


void startFile(char* fileName) {
  FILE* file = fopen(fileName, "w");
  fprintf(file,"# Higgs masses from large Nf analysis\n");
  fclose(file);
}


void calcDerivativeOhPhysicalMhWithRespectToRenCoup(ParameterType p, int n, double h, double &derive, double &atLamRen) {
  double L0 = p.lam0-0.5*h;
  double L1 = p.lam0;
  double L2 = p.lam0+0.5*h;

  if (n==0) {
    double m0Sqr = NaN;
    double mSqr = NaN;
    double lamRen = NaN;
    calcHiggsMasses(p, m0Sqr, mSqr);
    lamRen = calcRenormalizedQurticCoupling(p, m0Sqr, mSqr);
    atLamRen = lamRen;
    derive = p.CutoffInGev*sqrt(mSqr);
  } else if (n==1) {
    double m0Sqr0 = NaN;
    double mSqr0 = NaN;
    double lamRen0 = NaN;
    double m0Sqr2 = NaN;
    double mSqr2 = NaN;
    double lamRen2 = NaN;

    p.lam0 = L2;
    calcHiggsMasses(p, m0Sqr2, mSqr2);
    lamRen2 = calcRenormalizedQurticCoupling(p, m0Sqr2, mSqr2);

    p.lam0 = L0;
    calcHiggsMasses(p, m0Sqr0, mSqr0);
    lamRen0 = calcRenormalizedQurticCoupling(p, m0Sqr0, mSqr0);
    p.lam0 = L1;

    atLamRen = 0.5*(lamRen2+lamRen0);
    derive = (p.CutoffInGev*sqrt(mSqr2) - p.CutoffInGev*sqrt(mSqr0)) / (lamRen2-lamRen0);
  } else {
    double derive2, derive0, atLamRen2, atLamRen0;
    p.lam0 = L2;
    calcDerivativeOhPhysicalMhWithRespectToRenCoup(p, n-1, h, derive2, atLamRen2);

    p.lam0 = L0;
    calcDerivativeOhPhysicalMhWithRespectToRenCoup(p, n-1, h, derive0, atLamRen0);
    p.lam0 = L1;
    
    atLamRen = 0.5*(atLamRen2+atLamRen0);
    derive = (derive2-derive0) / (atLamRen2-atLamRen0);
  }
}


void calcDerivativeOfLatticeMHWithRespectToLambdaFunction(double (*func)(double x, double para), double LamFuncPara, ParameterType p, int n, double h, double &derive, double &atLamFunc) {
  double L0 = p.lam0-0.5*h;
  double L1 = p.lam0;
  double L2 = p.lam0+0.5*h;

  if (n==0) {
    double m0Sqr = NaN;
    double mSqr = NaN;
    double lamFunc = NaN;
    calcHiggsMasses(p, m0Sqr, mSqr);
    
    lamFunc = (*func) (p.lam0, LamFuncPara);
    atLamFunc = lamFunc;
    derive = sqrt(mSqr);
  } else if (n==1) {
    double m0Sqr0 = NaN;
    double mSqr0 = NaN;
    double lamFunc0 = NaN;
    double m0Sqr2 = NaN;
    double mSqr2 = NaN;
    double lamFunc2 = NaN;

    p.lam0 = L2;
    calcHiggsMasses(p, m0Sqr2, mSqr2);
    lamFunc2 = (*func) (p.lam0, LamFuncPara);

    p.lam0 = L0;
    calcHiggsMasses(p, m0Sqr0, mSqr0);
    lamFunc0 = (*func) (p.lam0, LamFuncPara);
    p.lam0 = L1;

    atLamFunc = 0.5*(lamFunc2+lamFunc0);
    derive = (sqrt(mSqr2)-sqrt(mSqr0)) / (lamFunc2-lamFunc0);    
  } else {
    double derive2, derive0, atLamFunc2, atLamFunc0;
    p.lam0 = L2;
    calcDerivativeOfLatticeMHWithRespectToLambdaFunction(func, LamFuncPara, p, n-1, h, derive2, atLamFunc2);

    p.lam0 = L0;
    calcDerivativeOfLatticeMHWithRespectToLambdaFunction(func, LamFuncPara, p, n-1, h, derive0, atLamFunc0);
    p.lam0 = L1;
    
    atLamFunc = 0.5*(atLamFunc2+atLamFunc0);
    derive = (derive2-derive0) / (atLamFunc2-atLamFunc0);
  }
}


double calcLatticeMhAtLargeQuarticCoupling(ParameterType p, double (*func)(double x, double para), double LargeLambda, double &A0, double &A1, double &A2) {
  printf("Calculating Higgs mass at large lambda=%f\n",LargeLambda);
  double para1 = 0;
  double para2 = 1.0;
  double atLamFunc = NaN;    
  double LamFuncPara = NaN;
  double h = 1E-5;
  double lamMerker = p.lam0;
  p.lam0 = 3*h;
  if (func!=&LambdaFunctionLuescher) {
    while (para2-para1>1E-4) {
      LamFuncPara = 0.5*(para1+para2);
      double derive2 = NaN;
      double derive = NaN;
      calcDerivativeOfLatticeMHWithRespectToLambdaFunction(func, LamFuncPara, p, 2,  h, derive, atLamFunc);
      if (derive == 0) break;
      calcDerivativeOfLatticeMHWithRespectToLambdaFunction(func, para2, p, 2,  h, derive2, atLamFunc);

      if (derive == 0) {
        LamFuncPara = para2;    
        break;
      }
    
      if (derive2*derive < 0) {
        para1 = LamFuncPara;
      } else {
        para2 = LamFuncPara;
      }
    }    
    printf("Optimal setting for Lambda-Function parameter: %f\n", LamFuncPara);
  }

  double lambdaBar = 1;
  if (!isNaN(LargeLambda)) lambdaBar = (*func) (LargeLambda, LamFuncPara);

  printf("Value for lambda bar is %f\n", lambdaBar);
  
  if (func!=&LambdaFunctionLuescher) {
    A0 = NaN;
    calcDerivativeOfLatticeMHWithRespectToLambdaFunction(func, LamFuncPara, p, 0,  h, A0, atLamFunc);
    A1 = NaN;
    calcDerivativeOfLatticeMHWithRespectToLambdaFunction(func, LamFuncPara, p, 1,  h, A1, atLamFunc);
    p.lam0 = lamMerker;
    A2 = NaN;
    return A0 + A1*lambdaBar;
  } else {
    A0 = NaN;
    calcDerivativeOfLatticeMHWithRespectToLambdaFunction(func, LamFuncPara, p, 0,  h, A0, atLamFunc);
    A1 = NaN;
    calcDerivativeOfLatticeMHWithRespectToLambdaFunction(func, LamFuncPara, p, 1,  h, A1, atLamFunc);
    A2 = NaN;
    calcDerivativeOfLatticeMHWithRespectToLambdaFunction(func, LamFuncPara, p, 2,  h, A2, atLamFunc);
    p.lam0 = lamMerker;
    
printf("%f %f  %f \n", A0, A1, A2);    
    
    return A0 + A1*lambdaBar;
//    return A0 + A1*lambdaBar + 0.5*A2*lambdaBar*lambdaBar;
  }
}


void appendMassesToFile(ParameterType p, char* fileName) {
  FILE* file = fopen(fileName, "a");

printf("append...");
  double m0Sqr = NaN;
  double mSqr = NaN;
  calcHiggsMasses(p, m0Sqr, mSqr);
  double lamRen = calcRenormalizedQurticCoupling(p, m0Sqr, mSqr);
  double A0, A1, A2;
  double mHAtLamLarge = calcLatticeMhAtLargeQuarticCoupling(p, &LambdaFunctionLuescher, NaN, A0, A1, A2);
  double lamFuncLuescher = LambdaFunctionLuescher(p.lam0, 0);
  
printf("Masses: %f %f\n ", m0Sqr, mSqr);  
printf("LamRen: %f \n ", lamRen);  
printf("MH at lambda=%f is m=%f\n", NaN, mHAtLamLarge);
printf("Lambda-Function at lambda=%f is m=%f\n", p.lam0, lamFuncLuescher);


  
  
  double topMass = NaN;
  double bottomMass = NaN;
  #ifdef CALCULATEFERMIONMASSES
    calcFermionMasses(p, m0Sqr, mSqr, topMass, bottomMass);
  #endif
  
  SimulationParameterSet* simPara = new SimulationParameterSet(m0Sqr, p.lam0, p.yt0, p.Nf, SimulationParameterSet_ContinuumNotation);
  fprintf(file,"%1.5e %1.5e %d %d %d %d %e %e %e %e %d %e %e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e\n", 
                m0Sqr, mSqr, p.L0,p.L1,p.L2,p.L3, p.yt0,p.yb0,p.lam0,p.CutoffInGev,p.Nf,p.rho,p.r,topMass, bottomMass, 
		simPara->getKappaN(), simPara->getYN(), simPara->getLambdaN(), p.lam6, p.lam8, p.lam10, lamRen, 
		mHAtLamLarge, lamFuncLuescher, A0, A1, A2);
  fclose(file);
  delete simPara;
}


void printParametersToScreen() {
  double m0Sqr = 0;
  double mSqr = 0;
  calcHiggsMasses(Parameters, m0Sqr, mSqr);
  
  SimulationParameterSet* simPara = new SimulationParameterSet(m0Sqr, Parameters.lam0, Parameters.yt0, Parameters.Nf, SimulationParameterSet_ContinuumNotation);  
  double fac = simPara->reparametrize_HiggsField(SimulationParameterSet_NfNotation);
  
  printf("\n   *** Parameters ***\n");
  printf("LatVol   = %dx%dx%dx%d\n",Parameters.L0, Parameters.L1, Parameters.L2, Parameters.L3);    
  printf("Cutoff   = %1.1f\n",Parameters.CutoffInGev);    
  printf("mSqr     = %1.5f,   --> m        = %1.2f GeV\n",mSqr,sqrt(mSqr)*Parameters.CutoffInGev);  
  printf("m0Sqr    = %1.5f,   --> kappa    = %1.8f\n",m0Sqr,simPara->getKappaN());
  printf("y        = %1.5f,   --> y        = %1.8f\n",Parameters.yt0, simPara->getYN());
  printf("lambda   = %1.5f,   --> lambda   = %1.8f\n",Parameters.lam0, simPara->getLambdaN());
  printf("lambda6  = %1.5f,   --> lambda6  = %1.8f\n",Parameters.lam6, Parameters.lam6/(fac*fac*fac*fac*fac*fac));
  printf("lambda8  = %1.5f,   --> lambda8  = %1.10f\n",Parameters.lam8, Parameters.lam8/(fac*fac*fac*fac*fac*fac*fac*fac));
  printf("lambda10 = %1.5f,   --> lambda10 = %1.10f\n\n",Parameters.lam10, Parameters.lam10/(fac*fac*fac*fac*fac*fac*fac*fac*fac*fac));
  
  delete simPara;
}


void SetCutOffToGivenKappa(double kappa, double cutoffBoundLow, double cutoffBoundHigh) {
  while (cutoffBoundHigh-cutoffBoundLow > 1E-10) {
    Parameters.CutoffInGev = 0.5 * (cutoffBoundLow + cutoffBoundHigh);
    double m0Sqr = 0;
    double mSqr = 0;
    calcHiggsMasses(Parameters, m0Sqr, mSqr);
  
    SimulationParameterSet* simPara = new SimulationParameterSet(m0Sqr, Parameters.lam0, Parameters.yt0, Parameters.Nf, SimulationParameterSet_ContinuumNotation);  
    double k = simPara->getKappaN();
    delete simPara;
    if (k>kappa) cutoffBoundLow = Parameters.CutoffInGev;
    if (k<kappa) cutoffBoundHigh = Parameters.CutoffInGev;
    if (k==kappa) break;
  }
}


void SetLambdaToGivenHiggsMass(double mSqrTargeted, double lambdaBoundLow, double lambdaBoundHigh) {
  while (lambdaBoundHigh-lambdaBoundLow > 1E-10) {
    Parameters.lam0 = 0.5 * (lambdaBoundLow + lambdaBoundHigh);
    double m0Sqr = 0;
    double mSqr = 0;
    calcHiggsMasses(Parameters, m0Sqr, mSqr);
  
    if (mSqr>mSqrTargeted) lambdaBoundHigh = Parameters.lam0;
    if (mSqr<mSqrTargeted) lambdaBoundLow = Parameters.lam0;
    if (mSqr==mSqrTargeted) break;
  }
}


void calcCoefficientsOfCouplingTermRecursion(int n, int N, int* facLine, int &resultCount, int* resultStore) {
  if (n==N) {
    int typeCounter[5];
    for (int I=0; I<5; I++) typeCounter[I]=0;
    for (int I=0; I<N; I++) typeCounter[facLine[I]]++;

    bool allEven = true;    
    int ID = typeCounter[0] + 100*typeCounter[1] + 10000*(typeCounter[2]+typeCounter[3]+typeCounter[4]);
    for (int I=0; I<5; I++) {
      if ((typeCounter[I] % 2) == 1) allEven = false;
    }
        
    if (allEven) {
      int mul = 1;
      for (int I=1; I<5; I++) {
        if (typeCounter[I]>0) mul *= NumberOfContractionsForRealScalarField(typeCounter[I]);
      }
      
      bool entryFound = false;
      for (int I=0; I<resultCount; I++) if (resultStore[2*I]==ID) {
        entryFound = true;
        resultStore[2*I+1] += mul;
	break;
      }
      if (!entryFound) {
        resultStore[2*resultCount+0] = ID;
        resultStore[2*resultCount+1] = mul;
        resultCount++;
      }
    }
  
    return;
  }
  facLine[n+0] = 0;
  facLine[n+1] = 0;
  calcCoefficientsOfCouplingTermRecursion(n+2, N, facLine, resultCount, resultStore);
  
  facLine[n+0] = 0;
  facLine[n+1] = 1;
  calcCoefficientsOfCouplingTermRecursion(n+2, N, facLine, resultCount, resultStore);
  
  facLine[n+0] = 1;
  facLine[n+1] = 0;
  calcCoefficientsOfCouplingTermRecursion(n+2, N, facLine, resultCount, resultStore);

  facLine[n+0] = 1;
  facLine[n+1] = 1;
  calcCoefficientsOfCouplingTermRecursion(n+2, N, facLine, resultCount, resultStore);

  facLine[n+0] = 2;
  facLine[n+1] = 2;
  calcCoefficientsOfCouplingTermRecursion(n+2, N, facLine, resultCount, resultStore);

  facLine[n+0] = 3;
  facLine[n+1] = 3;
  calcCoefficientsOfCouplingTermRecursion(n+2, N, facLine, resultCount, resultStore);

  facLine[n+0] = 4;
  facLine[n+1] = 4;
  calcCoefficientsOfCouplingTermRecursion(n+2, N, facLine, resultCount, resultStore);
}


void calcCoefficientsOfCouplingTerm(int N) {
  //printf("Calculating coefficients for coupling term Phi^%d...\n",N);
  int* facLine = new int[N];
  int* resultStore = new int[2000];
  int resultCount = 0;
  calcCoefficientsOfCouplingTermRecursion(0, N, facLine, resultCount, resultStore);
  CouplingTermCoeffcientsCount[N] = resultCount;
  
  for (int I=0; I<resultCount; I++) {
    int expCode = resultStore[2*I+0];
    int expV = expCode % 100;
    expCode /= 100;
    int expH = expCode % 100;
    expCode /= 100;
    int expG = expCode % 100;
    
   /* printf("  ");
    if (expV>0) printf(" v^%d ", expV); else printf("     ");
    if (expH>0) printf(" h^%d ", expH); else printf("     ");
    if (expG>0) printf(" g^%d ", expG); else printf("     ");
    printf(" * %d\n",resultStore[2*I+1]);  
   */ 
    CouplingTermCoeffcients[N][I][0][0] = expV;      
    CouplingTermCoeffcients[N][I][0][1] = expH;      
    CouplingTermCoeffcients[N][I][0][2] = expG;      
    CouplingTermCoeffcients[N][I][0][3] = resultStore[2*I+1];      
    for (int I2=1; I2<5; I2++) {
      if (CouplingTermCoeffcients[N][I][I2-1][0]<=0) {
        CouplingTermCoeffcients[N][I][I2][0] = 0;
        CouplingTermCoeffcients[N][I][I2][1] = CouplingTermCoeffcients[N][I][I2-1][1];    
        CouplingTermCoeffcients[N][I][I2][2] = CouplingTermCoeffcients[N][I][I2-1][2];     
        CouplingTermCoeffcients[N][I][I2][3] = 0;     
      } else {
        CouplingTermCoeffcients[N][I][I2][0] = CouplingTermCoeffcients[N][I][I2-1][0]-1;
        CouplingTermCoeffcients[N][I][I2][1] = CouplingTermCoeffcients[N][I][I2-1][1];    
        CouplingTermCoeffcients[N][I][I2][2] = CouplingTermCoeffcients[N][I][I2-1][2];     
        CouplingTermCoeffcients[N][I][I2][3] = CouplingTermCoeffcients[N][I][I2-1][3] * (CouplingTermCoeffcients[N][I][I2-1][0]);     
      }    
    }
  }  
    
  delete[] facLine;
  delete[] resultStore;
}



double calcFermionTresWithZeroExternalMomentaAndNLegs(ParameterType p, int n) {
  if ((p.yt0==0) && (p.yb0==0)) {
    return 0;
  }

  if (!QuietMode) printf("Calculating Fermion Tres with zero external momenta and %d legs\n", n);
  double res = 0;
  double fac = 0.5 / p.rho;
  double v = Physical_VEV_GeV / p.CutoffInGev;
  vector4D k;
  
  if ((p.L0>=0) && (p.L1>=0) && (p.L2>=0) && (p.L3>=0)) {
    for (int i0=0; i0<p.L0; i0++) {
      k[0] = 2*pi*i0 / p.L0;
      for (int i1=0; i1<p.L1; i1++) {
        k[1] = 2*pi*i1 / p.L1;
        for (int i2=0; i2<p.L2; i2++) {
          k[2] = 2*pi*i2 / p.L2;
          for (int i3=0; i3<p.L3; i3++) {
            k[3] = 2*pi*i3 / p.L3;

            Complex ew = diracOp->analyticalEigenvalue(k);
	    
	    Complex gamma = ComplexUnity - fac*ew;
            Complex dummyt = p.yt0 * gamma / (ew + p.yt0*v*gamma);
            Complex dummyb = p.yb0 * gamma / (ew + p.yb0*v*gamma);
	    Complex mult(1,0);
	    Complex mulb(1,0);
	    for (int I=0; I<n; I++) {
	      mult = mult * dummyt;
	      mulb = mulb * dummyb;	      
	    }
	    res += 4*mult.x + 4*mulb.x;
   	  }
        }
      }
    }
    res /= p.L0*p.L1*p.L2*p.L3;  
  } else {
    printf("Using Integration not implemented\n");
    exit(0);
  }
  
  return res;
}


void calcEffectivePotential(ParameterType p) {

double xxx = calcFermionTresWithZeroExternalMomentaAndNLegs(p, 2);
printf("%f\n",xxx);
xxx = calcFermionTresWithZeroExternalMomentaAndNLegs(p, 3);
printf("%f\n",xxx);
xxx = calcFermionTresWithZeroExternalMomentaAndNLegs(p, 4);
printf("%f\n",xxx);
}


void integrateMomentsOfV(ParameterType p, double start, double end, double Tol, double &v0, double &v1, double &v2) {
  v0 = 0;
  v1 = 0;
  v2 = 0;
  double Vol =p.L0*p.L1*p.L2*p.L3;
  double v0Old = NaN;
  double v1Old = NaN;
  double v2Old = NaN;
  double deriveData[100000][2+SecondOrderBosonicEffectivePotential_MAXDERIVE];
  int deriveDataCount = 0;
  int steps = 10;
  double h = (end-start)/steps;
  double FU, FUprime, FUprimeprime;
  double CutOffMerker = p.CutoffInGev;
  double v = Physical_VEV_GeV / p.CutoffInGev;

  calcFermionUContribution(p, FU, FUprime, FUprimeprime);
  makeSecondOrderBosonicPotentialAvail(p);
  secondOrderBosonicPotential->setVeV(v);

  double S0 = secondOrderBosonicPotential->calcEffectivePotDn(0) + FU;
  
  bool neglectable = true;
  for (int I=0; I<steps; I++) {
    double x = start + (0.5+I)*h;
    p.CutoffInGev = Physical_VEV_GeV / x;
    deriveData[I][0] = x;  
    secondOrderBosonicPotential->setVeV(x);
    calcFermionUContribution(p, FU, FUprime, FUprimeprime);
    
    deriveData[I][1] = exp(-Vol*(FU+secondOrderBosonicPotential->calcEffectivePotDn(0)-S0));
    if (deriveData[I][1] > 1E-7) neglectable = false;
    deriveDataCount++;
  } 

  if (neglectable) {
    p.CutoffInGev = CutOffMerker;    
    return;
  }

  while (true) {
    v0 = 0;
    v1 = 0;
    v2 = 0;  
    for (int I=0; I<deriveDataCount; I++) {
      v0 += h*deriveData[I][1];
      v1 += h*deriveData[I][0]*deriveData[I][1];
      v2 += h*sqr(deriveData[I][0])*deriveData[I][1];
    }

    if ((!isNaN(v0Old)) && (!isNaN(v1Old)) && (!isNaN(v2Old)) && (fabs((v0Old-v0)/v0)<Tol) && (fabs((v1Old-v1)/v1)<Tol) && (fabs((v2Old-v2)/v2)<Tol)) {
      break;
    }

    v0Old = v0;
    v1Old = v1;
    v2Old = v2;

    steps*= 3;
    h /= 3;
    for (int I=0; I<steps; I++) {
      if ((I%3) != 1) {
        double x = start + (0.5+I)*h;
        deriveData[deriveDataCount][0] = x;  
        secondOrderBosonicPotential->setVeV(x);
        p.CutoffInGev = Physical_VEV_GeV / x;	
        calcFermionUContribution(p, FU, FUprime, FUprimeprime);
        deriveData[deriveDataCount][1] = exp(-Vol*(FU+secondOrderBosonicPotential->calcEffectivePotDn(0)-S0));

        deriveDataCount++;
      }
    } 
  }
  p.CutoffInGev = CutOffMerker;
}


void calcActualVeVandMassSqr(ParameterType p, double& actVeV, double &actMSqr) {
  double sum0 = 0;
  double sum1 = 0;
  double sum2 = 0;
  
  for (int I=0; I<20; I++) {
    double v0 = 0;
    double v1 = 0;
    double v2 = 0;
    
    integrateMomentsOfV(p, 0.2*I, 0.2*(I+1), 1E-2, v0, v1, v2);
    sum0 += v0;
    sum1 += v1;
    sum2 += v2;
    
    printf("%f %f %f %f %f \n", I*0.2, (I+1)*0.2, v0, v1, v2); 
  }
  
  actVeV = sum1 / sum0;
  actMSqr = 1.0 / ((p.L0*p.L1*p.L2*p.L3)*(sum2 / sum0 - sqr(actVeV)));
}


void calcLambdaRrunningCoeff(double kappa1, double kappa2, double mag1, double mag2, double mH1, double mH2, double& A, double& B) {
  double v2 = sqrt(2*kappa2) * mag2;
  Parameters.CutoffInGev = 246 / v2;  
  SetLambdaToGivenHiggsMass(sqr(mH2), 0, 2.0);            
  double lam2 = Parameters.lam0;
  
  double v1 = sqrt(2*kappa1) * mag1;
  Parameters.CutoffInGev = 246 / v1;  
  SetLambdaToGivenHiggsMass(sqr(mH1), 0, 2.0);            
  double lam1 = Parameters.lam0;

  double C = (lam1/lam2);
  B = (C*log(1/sqr(v1)) - log(1/sqr(v2))) / (1-C);
  A = lam1*(log(1/sqr(v1))+B);
}



double calcLambdaRenMassDeviation(double* para) {
  FILE* file = fopen("HiggsMassFromLargeNf_LambdaRenMassDeviation.dat", "r");
  if (file==NULL) return 0;
  double res = 0;
  
  int L0, L1, L2, L3;
  double kap, mag, mH, mHerr;
  int count = 0;
  while (fscanf(file, "%d %d %d %d %lf %lf %lf %lf\n", &L0, &L1, &L2, &L3, &kap, &mag, &mH, &mHerr) == 8) {
    Parameters.L0 = L0;
    Parameters.L1 = L1;
    Parameters.L2 = L2;
    Parameters.L3 = L3;
    
    double v = sqrt(2*kap) * mag;
    Parameters.CutoffInGev = Physical_VEV_GeV / v;
    double A = para[0];
    double B = para[1];
    Parameters.lam0 = A/(log(1/sqr(v))+B);
    double m0Sqr, mSqr;
    calcHiggsMasses(Parameters, m0Sqr, mSqr);

    double chiSqr = sqr((sqrt(mSqr)-mH)/mHerr);
    res += chiSqr;
    count++;
  }
  res /= count;
  
  fclose(file);
  printf("ChiScr: %f, A: %f, B:%f\n", res, para[0], para[1]);    
  return res;
}


void calcLambdaRrunningCoeffByMinimization(double& A, double& B) {
  double pos[2];
  pos[0]=A;
  pos[1]=B;
  GradientMinimization(&calcLambdaRenMassDeviation, 2, 2E-2, 1E-2, 1E-4, pos, NULL, NULL, NULL, 3, 500);
  A=pos[0];
  B=pos[1];
}



int main(int argc,char **argv) {
  fftw_init_threads();
  fftw_plan_with_nthreads(1); 

	int L; 
	double mt;
	double cut; 
	
	try {
	if (argc != 4)
		throw std::string("usage: ./<progname> <L> <y> <Lambda>"); 

	if (sscanf(argv[1],"%d",&L)!=1)
		throw std::string("couldn't read L");
	if (sscanf(argv[2],"%lf",&mt)!=1)
		throw std::string("couldn't read y");
	if (sscanf(argv[3],"%lf",&cut)!=1)
		throw std::string("couldn't read cut");
	} 
	catch( const std::string& e) { 
		std::cout << "couldn't read input: " << e << std::endl;
		exit(1);
	}

  iniTools(5517);
  Parameters.L0 = L;
  Parameters.L1 = L;
  Parameters.L2 = L;
  Parameters.L3 = 2*L;

  Parameters.yt0 = mt / Physical_VEV_GeV;
  Parameters.yb0 = mt / Physical_VEV_GeV;
  Parameters.lam0 = 0.00;//1.84492e-02;//5.90109e-03  ;//1E-3;//-0.012943;//-0.010793   ;//1.000000e-03;//-0.004951642465726;//0.00466
  Parameters.lam6 = 0.00;//-8.29066e-03;//-2.90575e-03  ;//0E-3;//-0.001397;// -0.002700;//0.000316605203278;
  Parameters.lam8 = 0.00;//1.00000e-03;//3.00000e-04  ;//0E-3;//0.000743;// 0.000948;
  Parameters.lam10 = 0.00;
  Parameters.CutoffInGev = cut; //400.0;
  Parameters.Nf = 1;
  Parameters.rho = 1.0;
  Parameters.r = 0.5; 
  diracOp = new NeubergerMatrix(Parameters.rho, Parameters.r, 1, 1, 1, 1, 2);
  QuietMode = true;
//0.180173 -0.010793 -0.002700 0.000948 nahezu keine masse bei 8x8x8x16
//0.181642 -0.012943 -0.001397 0.000743 für 10GeV   bei 8x8x8x16


  calcCoefficientsOfCouplingTerm(2);
  calcCoefficientsOfCouplingTerm(4);
  calcCoefficientsOfCouplingTerm(6);
  calcCoefficientsOfCouplingTerm(8);
  calcCoefficientsOfCouplingTerm(10);

//randomScanParameterSpace(5, 1E-3, 0E-4, 1);

//plotActionDependenceOnV();
//SetCutOffToGivenKappa(0.12313, 100, 2000);
//L12 30400:   mH=4.341354147472885e-01    mag = 2.053220197441763e-01
//L12 30274:   mH=3.930393024638761e-01    mag = 1.877641010871061e-01
//L16 30400:   mH=4.251290301035624e-01    mag = 2.025470649697020e-01
//L16 30274:   mH=3.818779456487391e-01    mag = 1.833878260561906e-01



/*Parameters.CutoffInGev = 246 / (sqrt(2*0.30400) * 2.053220197441763e-01);   //L12
SetLambdaToGivenHiggsMass(sqr(4.341354147472885e-01), 0, 2.0);              //L12
double lam1 = Parameters.lam0;
Parameters.CutoffInGev = 246 / (sqrt(2*0.30274) * 1.877641010871061e-01);   //L12
SetLambdaToGivenHiggsMass(sqr(3.930393024638761e-01), 0, 2.0);              //L12
double lam2 = Parameters.lam0;*/


//double A,B;
//calcLambdaRrunningCoeff(0.30400, 0.30274, 2.053220197441763e-01, 1.877641010871061e-01, 4.341354147472885e-01, 3.930393024638761e-01, A, B); //L12 FULL HY
//calcLambdaRrunningCoeff(0.30400, 0.30274, 2.025470649697020e-01, 1.833878260561906e-01, 4.251290301035624e-01, 3.818779456487391e-01, A, B); //L16 FULL HY
//calcLambdaRrunningCoeff(0.30400, 0.30274, 2.004101585811714e-01, 1.810278304854023e-01, 4.177383481583092e-01, 3.714911141141783e-01, A, B); //L20 FULL HY
//calcLambdaRrunningCoeff(0.30400, 0.30274, 1.995396327700017e-01, 1.796833115512425e-01, 4.146274467639159e-01, 3.672473303449806e-01, A, B); //L24 FULL HY
//calcLambdaRrunningCoeff(0.30400, 0.30274, 1.983983360605696e-01, 1.780318696210267e-01, 4.085199355322402e-01, 3.583531597105321e-01, A, B); //L32 FULL HY



//calcLambdaRrunningCoeff(0.31040, 0.30890, 2.055475705117546e-01, 1.829693642889963e-01, 4.174238780105952e-01, 3.626606355246622e-01, A, B); //L12 PurePhi4
//calcLambdaRrunningCoeff(0.31040, 0.30890, 2.069541218219623e-01, 1.838036054983557e-01, 4.276213072880637e-01, 3.691221188346390e-01, A, B); //L16 PurePhi4
//calcLambdaRrunningCoeff(0.31040, 0.30890, 2.065097050440377e-01, 1.836166118295823e-01, 4.230652796395589e-01, 3.687846075918060e-01, A, B); //L20 PurePhi4
//calcLambdaRrunningCoeff(0.31040, 0.30890, 2.058764011627206e-01, 1.829946377983989e-01, 4.213365459290799e-01, 3.653785173955839e-01, A, B); //L24 PurePhi4
//calcLambdaRrunningCoeff(0.31040, 0.30890, 2.051086820373522e-01, 1.818968411433086e-01, 4.174746388162021e-01, 3.615047687490005e-01, A, B); //L32 PurePhi4







//calcLambdaRrunningCoeffByMinimization(A, B);

//printf("A=%f B=%f\n", A,B);


//ChiScr: 1.729763, A: 4.366795, B:1.282895  //L12
//A=4.348978 B=1.291982
//ChiScr: 1.932347, A: 3.553243, B:0.222814  //L16
//A=3.559706 B=0.215183



//Parameters.CutoffInGev = 246 / (sqrt(2*0.30400) * 2.025470649697020e-01);   //L16
//SetLambdaToGivenHiggsMass(sqr(4.251290301035624e-01), 0, 2.0);              //L16


//Parameters.CutoffInGev = 246 / (sqrt(2*0.30148) * 1.617281100482811e-01);   //L16 30148
//SetLambdaToGivenHiggsMass(sqr(3.294603235262912e-01), 0, 2.0);              //L16 30148
//3.559047084127087e-01  mag = 1.694273100596201e-01  //L12 30148
//3.294603235262912e-01  mag = 1.617281100482811e-01  //L16 3-148
printParametersToScreen();
exit(0);

  char* fileName = constructFileName(Parameters, "LambdaScan");
/*  startFile(fileName);
  for (Parameters.lam0=0.0; Parameters.lam0<=0.05; Parameters.lam0+=0.001) {
    appendMassesToFile(Parameters, fileName);
  }*/
  delete[] fileName;


  fileName = constructFileName(Parameters, "Lambda6Scan");
/*  startFile(fileName);
  for (Parameters.lam6=0.0; Parameters.lam6<=0.02; Parameters.lam6+=0.0002) {
    appendMassesToFile(Parameters, fileName);
  }*/
  delete[] fileName;

  fileName = constructFileName(Parameters, "Lambda8Scan");
/*  startFile(fileName);
  for (Parameters.lam8=0.0; Parameters.lam8<=0.02; Parameters.lam8+=0.0002) {
    appendMassesToFile(Parameters, fileName);
  }*/
  delete[] fileName;


  fileName = constructFileName(Parameters, "CutoffScan");
/*  startFile(fileName);
  for (Parameters.CutoffInGev=1000; Parameters.CutoffInGev<=4000; Parameters.CutoffInGev+=50) {
//double v = 246 / Parameters.CutoffInGev;
//Parameters.lam0 = A/(log(1/sqr(v))+B);
    appendMassesToFile(Parameters, fileName);
  }*/
  delete[] fileName;


  
  
  fileName = constructFileName(Parameters, "YbScan");
/*  startFile(fileName);
  for (Parameters.yb0=0.01; Parameters.yb0<=200.0/Physical_VEV_GeV; Parameters.yb0+=0.01) {
    appendMassesToFile(Parameters, fileName);
  }*/
  delete[] fileName;


  fileName = constructFileName(Parameters, "YtScan");
/*  startFile(fileName);
  for (Parameters.yt0=0.01; Parameters.yt0<=200.0/Physical_VEV_GeV; Parameters.yt0+=0.01) {
    appendMassesToFile(Parameters, fileName);
  }*/
  delete[] fileName;


  fileName = constructFileName(Parameters, "YScan");
/*  startFile(fileName);
  for (Parameters.yt0=0.01; Parameters.yt0<=200.0/Physical_VEV_GeV; Parameters.yt0+=0.01) {
    Parameters.yb0 = Parameters.yt0;
    appendMassesToFile(Parameters, fileName);
  }*/
  delete[] fileName;


  fileName = constructFileName(Parameters, "VolumeScan");
/*  startFile(fileName);
  for (Parameters.L0=4; Parameters.L0<=30; Parameters.L0+=2) {
    Parameters.L1 = Parameters.L0;
    Parameters.L2 = Parameters.L0;    
    appendMassesToFile(Parameters, fileName);
  }*/
  delete[] fileName;
  
  
  delete[] pSqr;
  delete diracOp;
  delete[] DiracOpEigenvalues;
  if (DiracOpEigenvectors!=NULL) {
    for (int I=0; I<Parameters.L0*Parameters.L1*Parameters.L2*Parameters.L3; I++) {
      delete[] DiracOpEigenvectors[I];
    }
  }
  delete[] DiracOpEigenvectors;
}
