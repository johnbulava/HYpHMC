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



#define INTEGRATORACCURACY 1E-5


struct ParameterType {
  double m;
  double s;
  double yt0;
  double yb0;
  double lam0;
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
double LastFermionUContributionResult = NaN;
double LastFermionUprimeMContributionResult = NaN;
double LastFermionUprimeSContributionResult = NaN;
ParameterType  LastCosKParameters;
double* CosK = NULL;
NeubergerMatrix* diracOp = NULL;
bool QuietMode = false;



char* constructFileName(ParameterType p, char* des) {
  char* fileName = new char[1000];
  snprintf(fileName, 1000, "PhaseDiagramAtLargeLambda_L%dx%dx%dx%dyt%1.3fyb%1.3flam%1.5ffNf%dRho%1.3fr%1.3f_%s.dat", p.L0,p.L1,p.L2,p.L3, p.yt0,p.yb0, p.lam0,p.Nf,p.rho,p.r,des);  
  return fileName;  
}


char* constructFileNameNoY(ParameterType p, char* des) {
  char* fileName = new char[1000];
  snprintf(fileName, 1000, "PhaseDiagramAtLargeLambda_L%dx%dx%dx%dybytRatio%1.3flam%1.5ffNf%dRho%1.3fr%1.3f_%s.dat", p.L0,p.L1,p.L2,p.L3, p.yb0/p.yt0, p.lam0,p.Nf,p.rho,p.r,des);  
  return fileName;  
}


bool areParametersEqualWithRespectToCosK(ParameterType p1, ParameterType p2) {
  if (p1.L0 != p2.L0) return false;
  if (p1.L1 != p2.L1) return false;
  if (p1.L2 != p2.L2) return false;
  if (p1.L3 != p2.L3) return false;

  return true;
}


void makeCosKAvailableForParameters(ParameterType p) {
  if ((CosK==NULL) || (!areParametersEqualWithRespectToCosK(p, LastCosKParameters))) {
    printf("Generating CosK\n");
    delete[] CosK;
    CosK = new double[p.L0 * p.L1 * p.L2 * p.L3];
    double* c0 = new double[p.L0];
    double* c1 = new double[p.L1];
    double* c2 = new double[p.L2];
    double* c3 = new double[p.L3];
    for (int I0=0; I0<p.L0; I0++) c0[I0] = cos((2*pi*I0)/p.L0);
    for (int I1=0; I1<p.L1; I1++) c1[I1] = cos((2*pi*I1)/p.L1);
    for (int I2=0; I2<p.L2; I2++) c2[I2] = cos((2*pi*I2)/p.L2);
    for (int I3=0; I3<p.L3; I3++) c3[I3] = cos((2*pi*I3)/p.L3);
    
    int count = 0;
    for (int I0=0; I0<p.L0; I0++) {
      for (int I1=0; I1<p.L1; I1++) {
        for (int I2=0; I2<p.L2; I2++) {
          for (int I3=0; I3<p.L3; I3++) {
            CosK[count] = c0[I0] + c1[I1] + c2[I2] + c3[I3];    
	    count++;
          }
        }
      }
    }
    delete[] c0;
    delete[] c1;
    delete[] c2;
    delete[] c3;
    LastCosKParameters = p;  
    printf("Ready\n");
  }
}


bool areParametersEqualWithRespectToFermionContribution(ParameterType p1, ParameterType p2) {
  if (p1.L0 != p2.L0) return false;
  if (p1.L1 != p2.L1) return false;
  if (p1.L2 != p2.L2) return false;
  if (p1.L3 != p2.L3) return false;
  if (p1.yt0 != p2.yt0) return false;
  if (p1.yb0 != p2.yb0) return false;
  if (p1.m != p2.m) return false;
  if (p1.s != p2.s) return false;  
  if (p1.rho != p2.rho) return false;
  if (p1.r != p2.r) return false;

  return true;
}


/*double calcFermionUContribution_IntegrandU(double p0, double p1, double p2, double p3, double para) {
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
}*/


void calcFermionUContribution(ParameterType p, double &U, double &UprimeM, double &UprimeS) {
  if ((p.yt0==0) && (p.yb0==0)) {
    U = 0;
    UprimeM = 0;
    UprimeS = 0;
    return;
  }

  if ((!isNaN(LastFermionUContributionResult)) && (!isNaN(LastFermionUprimeMContributionResult)) && (!isNaN(LastFermionUprimeSContributionResult))) {
    if (areParametersEqualWithRespectToFermionContribution(p, LastFermionContributionParameters)) {
      U = LastFermionUContributionResult;
      UprimeM = LastFermionUprimeMContributionResult;
      UprimeS = LastFermionUprimeSContributionResult;
      return;
    }
  }
  
  if (!QuietMode) printf("Calculating Fermion Contribution\n");
  U = 0;
  UprimeM = 0;
  UprimeS = 0;
  vector4D k;  
  vector4D kpi;  
  Complex tworho = Complex(2*p.rho,0);
  double m =  p.m;
  double s =  p.s;  
  double ytSqrd4rhosqr = p.yt0*p.yt0 / (4*p.rho*p.rho);
  double ybSqrd4rhosqr = p.yb0*p.yb0 / (4*p.rho*p.rho);
  
  if ((p.L0>=0) && (p.L1>=0) && (p.L2>=0) && (p.L3>=0)) {
    for (int i0=0; i0<p.L0; i0++) {
      k[0] = 2*pi*i0 / p.L0;
      kpi[0] = pi+k[0];
      for (int i1=0; i1<p.L1; i1++) {
        k[1] = 2*pi*i1 / p.L1;
        kpi[1] = pi+k[1];
        for (int i2=0; i2<p.L2; i2++) {
          k[2] = 2*pi*i2 / p.L2;
          kpi[2] = pi+k[2];
          for (int i3=0; i3<p.L3; i3++) {
            k[3] = 2*pi*i3 / p.L3;
            kpi[3] = pi+k[3];

            Complex ew = diracOp->analyticalEigenvalue(k);
            Complex ewpi = diracOp->analyticalEigenvalue(kpi);
	    
	    double nk = norm(ew);
	    double nkpi = norm(ewpi);
	    double nkm2rho = norm(ew-tworho);
	    double nkpim2rho = norm(ewpi-tworho);
      
            double a1 = nk*nkpi;
	    double a2 = nkm2rho*nkpim2rho;
	    double a3 = sqr(nkm2rho*nkpi - nkpim2rho*nk);
	    
	    double nt1 = a1 + ytSqrd4rhosqr*(m*m-s*s)*a2;
	    double nb1 = a1 + ybSqrd4rhosqr*(m*m-s*s)*a2;	    
	    double nt2 = sqr(nt1) + m*m*ytSqrd4rhosqr*a3;
	    double nb2 = sqr(nb1) + m*m*ybSqrd4rhosqr*a3;

            U += -log(nt2);
            U += -log(nb2);
	    
	    UprimeM += -(2*nt1*2*m*ytSqrd4rhosqr*a2 + 2*m*ytSqrd4rhosqr*a3) / nt2;
	    UprimeM += -(2*nb1*2*m*ybSqrd4rhosqr*a2 + 2*m*ybSqrd4rhosqr*a3) / nb2;
	    
	    UprimeS += -(-2*nt1*2*s*ytSqrd4rhosqr*a2) / nt2;
	    UprimeS += -(-2*nb1*2*s*ybSqrd4rhosqr*a2) / nb2;
   	  }
        }
      }
    }
    U /= p.L0*p.L1*p.L2*p.L3;  
    UprimeM /= p.L0*p.L1*p.L2*p.L3;  
    UprimeS /= p.L0*p.L1*p.L2*p.L3;  
  } else {
    printf("Using Integration\n");
    vector4D start;
    vector4D end;
    for (int I=0; I<4; I++) {
      start[I] = 0;
      end[I] = pi;
    }

/*    U = sqr(sqr(1.0/pi)) * integrate(&calcFermionUContribution_IntegrandU, NaN, start, end, INTEGRATORACCURACY);
    UprimeM = sqr(sqr(1.0/pi)) * integrate(&calcFermionUContribution_IntegrandUprime, NaN, start, end, INTEGRATORACCURACY);
    UprimeS = sqr(sqr(1.0/pi)) * integrate(&calcFermionUContribution_IntegrandUprimeprime, NaN, start, end, INTEGRATORACCURACY);*/
  }
   
  LastFermionContributionParameters = p;
  LastFermionUContributionResult = U;
  LastFermionUprimeMContributionResult = UprimeM;
  LastFermionUprimeSContributionResult = UprimeS;
}


void startFile(char* fileName) {
  FILE* file = fopen(fileName, "w");
  fprintf(file,"# Phase diagram at large lambda from large Nf analysis\n");
  fclose(file);
}


double calcBosonicGapTerm(ParameterType p, double kappa, double lam) {
  makeCosKAvailableForParameters(p);
  double res = 0;
  for (int I=1; I<p.L0*p.L1*p.L2*p.L3; I++) {
    double c = CosK[I];
    if (c>-4+1E-5) {
      double cont = 1/(-4*kappa*CosK[I] + 2*lam);
      if (cont<0) {
        printf("Negative contribtion in calcBosonicGapTerm (I: %d, kappa:%f, lam: %f)!\n", I, kappa, lam);
	exit(0);
      }
      res += cont;
    } else {
//      printf("neglect\n");
    }
  }
  return 4*res / (p.Nf*p.L0*p.L1*p.L2*p.L3);
}


double findKappa(ParameterType p, double minK, double maxK, double TOL) {  
  if ((fabs(p.m)>1E-10) && (fabs(p.s)>1E-10)) return NaN;
//printf("Min %f, max %f\n",minK,maxK);
  double U, UprimeM, UprimeS;
  calcFermionUContribution(p, U, UprimeM, UprimeS);
  double k0 = 0.5*(minK+maxK);
  if (fabs(minK-maxK)<TOL) return k0;
  
  double minKlam = (-UprimeM + 8*minK*2*p.m) / (2*p.m);
  double maxKlam = (-UprimeM + 8*maxK*2*p.m) / (2*p.m);
  double k0lam = (-UprimeM + 8*k0*2*p.m) / (2*p.m);
  if (fabs(p.s)>1E-10) {
    minKlam = (-UprimeS - 8*minK*2*p.s) / (2*p.s);
    maxKlam = (-UprimeS - 8*maxK*2*p.s) / (2*p.s);
    k0lam = (-UprimeS - 8*k0*2*p.s) / (2*p.s);
  }
  
  double minKres = 1 - p.m*p.m - p.s*p.s - calcBosonicGapTerm(p, minK, minKlam);
  double maxKres = 1 - p.m*p.m - p.s*p.s - calcBosonicGapTerm(p, maxK, maxKlam);
  double k0res = 1 - p.m*p.m - p.s*p.s - calcBosonicGapTerm(p, k0, k0lam);
//printf("res: %f %f %f\n", minKres, maxKres, k0res);  
  
  if (minKres*maxKres > 0) return NaN;
  if (fabs(minKres)<TOL) return minK;
  if (fabs(maxKres)<TOL) return maxK;
  if (fabs(k0res)<TOL) return k0;

  if (k0res*minKres > 0) return findKappa(p, k0, maxK, TOL);
  if (k0res*maxKres > 0) return findKappa(p, minK, k0, TOL);
  return NaN;
}


double calcKappaMRelation(ParameterType p, int N, double** &rel, double startM, double endM, double minK, double maxK) {
  rel = new double*[N];
  for (int I=0; I<N; I++) {
    rel[I] = new double[3];
    p.m = (endM-startM)*I/(N-1) + startM;
    p.s = 0;
    rel[I][0] = p.m;
    rel[I][1] = findKappa(p, minK, maxK, 1E-10);
  }
  double maxDer = 0;
  double kcrit = NaN;
  for (int I=1; I<N-1; I++) {
    rel[I][2] = (rel[I+1][0]-rel[I-1][0]) / (rel[I+1][1]-rel[I-1][1]);
    if (rel[I][2]>maxDer) {
      maxDer = rel[I][2];
      kcrit = rel[I][1];
    }
  }
  rel[0][2] = NaN;
  rel[N-1][2] = NaN;
  return kcrit;
}


void plotPhaseDiagram(ParameterType p, double minY, double maxY, double yRatio) {
  int N=100;
  int Nscan = 1000;
  char* fileName = constructFileNameNoY(p, "PhaseDiagram");
  FILE* file = fopen(fileName, "w");  
  
  for (int I=0; I<N; I++) {
    p.yt0 = (maxY-minY)*I/(N-1)+minY;
    p.yb0 = p.yt0 * yRatio;
    
    double** relation = NULL;
    double kappaCrit = calcKappaMRelation(p, Nscan, relation, 0.10, 0.60, 0.0, 1.0);
    delete[] relation;
    
    fprintf(file, "%1.15f %1.15f\n", p.yt0, kappaCrit);
  }
  fclose(file);
}


void plotMagnetizations(ParameterType p) {
  int Nscan = 1000;
  char* fileName = constructFileName(p, "Magnetizations");
  FILE* file = fopen(fileName, "w");  

  double** relation = NULL;
  calcKappaMRelation(p, Nscan, relation, 0.10, 0.90, 0.0, 1.0);
  for (int I=0; I<Nscan; I++) {
    if (!isNaN(relation[I][1])) fprintf(file, "%1.15f %1.15f\n", relation[I][1], relation[I][0]);
  }

  delete[] relation;
  fclose(file);
}



int main(int argc,char **argv) {
  iniTools(5517);
  Parameters.L0 = 8;
  Parameters.L1 = 8;
  Parameters.L2 = 8;
  Parameters.L3 = 8;

  Parameters.m = 0.2;
  Parameters.s = 0.0;  
  Parameters.yt0 = 1.0;
  Parameters.yb0 = 1.0;
  Parameters.lam0 = NaN;
  Parameters.Nf = 20;
  Parameters.rho = 1.0;
  Parameters.r = 0.5; 
  diracOp = new NeubergerMatrix(Parameters.rho, Parameters.r, 1, 1, 1, 1, 2);

  QuietMode = true;

/*  printf("Result for Kappa: %f\n", findKappa(Parameters, 0.0, 1.00, 1E-5));
  double** rel = NULL;
  printf("critical kappa: %f\n", calcKappaMRelation(Parameters, 10, rel, 0.10, 0.60, 0.0, 1.0));
  
  for (int I=0; I<10; I++) {
    printf("%f %f %f\n",rel[I][0], rel[I][1], rel[I][2]);
  }*/
  
  plotMagnetizations(Parameters);
  plotPhaseDiagram(Parameters, 0.001, 3.0, 1.0);

/*  char* fileName = constructFileName(Parameters, "LambdaScan");
  startFile(fileName);
  for (Parameters.lam0=0.0; Parameters.lam0<=0.03; Parameters.lam0+=0.001) {
    appendMassesToFile(Parameters, fileName);
  }
  delete[] fileName;*/


  
  
  delete[] CosK;
  delete diracOp;
}
