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
  double yt0;
  double yb0;
  double lam0;
  double lam6;
  double lam8;
  double lam10;  
  double m0Sqr;
  int Nf;
  int L0;
  int L1;
  int L2;
  int L3;
  double rho;
  double r;
};
ParameterType  Parameters;
ParameterType  LastPSqrParameters;
double* pSqr = NULL;
NeubergerMatrix* diracOp = NULL;
ParameterType  LastDiracOperatorEigenvaluesAndEigenvectorsParameters;
Complex* DiracOpEigenvalues = NULL;
ComplexVector** DiracOpEigenvectors = NULL;
bool QuietMode = false;



char* constructFileName(ParameterType p, const char* des) {
  char* fileName = new char[1000];
  snprintf(fileName, 1000, "PerturbativePredictions_L%dx%dx%dx%dyt%1.3fyb%1.3flam%1.5fm0Sqr%1.6fNf%dRho%1.3fr%1.3f_%s.dat",p.L0,p.L1,p.L2,p.L3, p.yt0,p.yb0,p.lam0,p.m0Sqr,p.Nf,p.rho,p.r,des);
  return fileName;  
}


void startFile(char* fileName) {
  FILE* file = fopen(fileName, "w");
  fprintf(file,"# Perturbative predictions from lattice perturbation theory\n");
  fclose(file);
}


bool areParametersEqualWithRespectToPSqr(ParameterType p1, ParameterType p2) {
  if (p1.L0 != p2.L0) return false;
  if (p1.L1 != p2.L1) return false;
  if (p1.L2 != p2.L2) return false;
  if (p1.L3 != p2.L3) return false;

  return true;
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


double calcSingleBosonicLoop(ParameterType p, double A) {
  double res = 0;
  if ((p.L0>=0) && (p.L1>=0) && (p.L2>=0) && (p.L3>=0)) {
    makePSqrAvailableForParameters(p);
      
    for (int I=0; I<p.L0 * p.L1 * p.L2 * p.L3; I++) {
      res += 1.0 / (pSqr[I] + A);
    }
    res /= p.L0 * p.L1 * p.L2 * p.L3;
    
  } else {  
    printf("Continous space-time not implemented!\n");
    exit(0);
  }  
  return res;
}


double calcSingleHiggsLoop(ParameterType p, double v) {
  return calcSingleBosonicLoop(p, p.m0Sqr+12*p.lam0*v*v);
}


double calcSingleGoldstoneLoop(ParameterType p, double v) {
  return calcSingleBosonicLoop(p, p.m0Sqr+4*p.lam0*v*v);
}


double calcDoubleBosonicLoop(ParameterType p, double A, int exP0, int exP1, int exP2, int exP3) {
  double res = 0;
  if ((p.L0>=0) && (p.L1>=0) && (p.L2>=0) && (p.L3>=0)) {
    makePSqrAvailableForParameters(p);

    int index1 = 0;
    int index2 = 0;    
    int j0,j1,j2,j3;
    for (int I0=0; I0<p.L0; I0++) {
      j0 = ((I0 - exP0 + 2*p.L0) % p.L0) * p.L1*p.L2*p.L3;
      for (int I1=0; I1<p.L1; I1++) {
        j1 = ((I1 - exP1 + 2*p.L1) % p.L1) * p.L2*p.L3;
        for (int I2=0; I2<p.L2; I2++) {
          j2 = ((I2 - exP2 + 2*p.L2) % p.L2) * p.L3;
          for (int I3=0; I3<p.L3; I3++) {
            j3 = (I3 - exP3 + 2*p.L3) % p.L3;
	    index2 = j0+j1+j2+j3;
	  
            res += 1.0 / ((pSqr[index1] + A) * (pSqr[index2] + A));
	    index1++;
          }
        }
      }
    }
    res /= p.L0 * p.L1 * p.L2 * p.L3;
  } else {  
    printf("Continous space-time not implemented!\n");
    exit(0);
  }  
  return res;
}


double calcDoubleHiggsLoop(ParameterType p, double v, int exP0, int exP1, int exP2, int exP3) {
  return calcDoubleBosonicLoop(p, p.m0Sqr+12*p.lam0*v*v, exP0, exP1, exP2, exP3);
}


double calcDoubleGoldstoneLoop(ParameterType p, double v, int exP0, int exP1, int exP2, int exP3) {
  return calcDoubleBosonicLoop(p, p.m0Sqr+4*p.lam0*v*v, exP0, exP1, exP2, exP3);
}


double calcSingleFermionLoop(ParameterType p, double v) {
  double res = 0;
  if ((p.yt0==0) && (p.yb0==0)) return res;

  if ((p.L0>=0) && (p.L1>=0) && (p.L2>=0) && (p.L3>=0)) {
    makeDiracOpEigenvaluesAndEigenvectorsAvailable(p);
    
    double yt = p.yt0;
    double yb = p.yb0;  
    double fac = 0.5 / p.rho;

    for (int I=0; I<p.L0 * p.L1 * p.L2 * p.L3; I++) {
      Complex ewK = DiracOpEigenvalues[I];
      Complex gamma = ComplexUnity - fac*ewK;
      Complex at = 4*yt*gamma / (ewK + yt*v*gamma);
      Complex ab = 4*yb*gamma / (ewK + yb*v*gamma);
      res += at.x + ab.x;
    }
    res /= p.L0 * p.L1 * p.L2 * p.L3;    
  } else {      
    printf("Continous space-time not implemented!\n");
    exit(0);
  }
  return res;
}


double calcDoubleFermionLoop(ParameterType p, double v, int Bindex, int exP0, int exP1, int exP2, int exP3) {
  double res = 0;
  if ((p.yt0==0) && (p.yb0==0)) return res;

  if ((p.L0>=0) && (p.L1>=0) && (p.L2>=0) && (p.L3>=0)) {
    makeDiracOpEigenvaluesAndEigenvectorsAvailable(p);
    
    double yt = Parameters.yt0;
    double yb = Parameters.yb0;  
    double fac = 0.5 / Parameters.rho;
    
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

    ComplexVector vK[4];
    vK[0].resize(4);
    vK[1].resize(4);
    vK[2].resize(4);
    vK[3].resize(4);
    ComplexVector vKP[4];
    vKP[0].resize(4);
    vKP[1].resize(4);
    vKP[2].resize(4);
    vKP[3].resize(4);


    int index1 = 0;
    int index2 = 0;    
    int j0,j1,j2,j3;
    Complex sum(0,0);
    for (int I0=0; I0<p.L0; I0++) {
      j0 = ((I0 - exP0 + 2*p.L0) % p.L0) * p.L1*p.L2*p.L3;
      for (int I1=0; I1<p.L1; I1++) {
        j1 = ((I1 - exP1 + 2*p.L1) % p.L1) * p.L2*p.L3;
        for (int I2=0; I2<p.L2; I2++) {
          j2 = ((I2 - exP2 + 2*p.L2) % p.L2) * p.L3;
          for (int I3=0; I3<p.L3; I3++) {
            j3 = (I3 - exP3 + 2*p.L3) % p.L3;
	    index2 = j0+j1+j2+j3;
	  
            Complex ewK = DiracOpEigenvalues[index1];
            Complex gammaK = ComplexUnity - fac*ewK;
            for (int I=0; I<4; I++) {
  	      vK[I] = DiracOpEigenvectors[index1][I];
            }
	  
            Complex ewKP = DiracOpEigenvalues[index2];
            Complex gammaKP = ComplexUnity - fac*ewKP;
	    for (int I=0; I<4; I++) {
  	      vKP[I] = DiracOpEigenvectors[index2][I];
            }
	  
            ComplexMatrix mK1(vK[0]);  
            ComplexMatrix mK2(vK[1]);  
            ComplexMatrix mK3(vK[2]);  
            ComplexMatrix mK4(vK[3]);  

            Complex ewGamDinvKtop    = gammaK / (ewK + yt*v*gammaK);
            Complex ewGamDinvKbottom = gammaK / (ewK + yb*v*gammaK);
            ComplexMatrix GamDinvKtop    = ewGamDinvKtop*mK1 + ewGamDinvKtop*mK2 + adj(ewGamDinvKtop)*mK3 + adj(ewGamDinvKtop)*mK4;
            ComplexMatrix GamDinvKbottom = ewGamDinvKbottom*mK1 + ewGamDinvKbottom*mK2 + adj(ewGamDinvKbottom)*mK3 + adj(ewGamDinvKbottom)*mK4;
            ComplexMatrix GamDinvK(8);
            GamDinvK.setZero();
            GamDinvK.insertMatrix(GamDinvKtop, 0, 0);
            GamDinvK.insertMatrix(GamDinvKbottom, 4, 4);
	  
            ComplexMatrix mKP1(vKP[0]);  
            ComplexMatrix mKP2(vKP[1]);  
            ComplexMatrix mKP3(vKP[2]);  
            ComplexMatrix mKP4(vKP[3]);  
	  
            Complex ewGamDinvKPtop    = gammaKP / (ewKP + yt*v*gammaKP);
            Complex ewGamDinvKPbottom = gammaKP / (ewKP + yb*v*gammaKP);
            ComplexMatrix GamDinvKPtop    = ewGamDinvKPtop*mKP1 + ewGamDinvKPtop*mKP2 + adj(ewGamDinvKPtop)*mKP3 + adj(ewGamDinvKPtop)*mKP4;
            ComplexMatrix GamDinvKPbottom = ewGamDinvKPbottom*mKP1 + ewGamDinvKPbottom*mKP2 + adj(ewGamDinvKPbottom)*mKP3 + adj(ewGamDinvKPbottom)*mKP4;
            ComplexMatrix GamDinvKP(8);
            GamDinvKP.setZero();
            GamDinvKP.insertMatrix(GamDinvKPtop, 0, 0);
            GamDinvKP.insertMatrix(GamDinvKPbottom, 4, 4);

            sum = sum + (BMat[Bindex]*GamDinvK*BMat[Bindex]*GamDinvKP).tres();
	    index1++;
          }
        }
      }
    }
    res = sum.x;
    res /= p.L0 * p.L1 * p.L2 * p.L3;
  } else {      
    printf("Continous space-time not implemented!\n");
    exit(0);
  }
  return res;
}


double calcExpectationOfHiggsField(ParameterType p, double v) {
  double fermionLoop = calcSingleFermionLoop(p, v);
  double HiggsLoop = calcSingleHiggsLoop(p, v);
  double GoldstoneLoop = calcSingleGoldstoneLoop(p, v);
  
  return (p.m0Sqr*v + 4*p.lam0*v*v*v - fermionLoop + 4*p.lam0*v*(3*HiggsLoop + 3*GoldstoneLoop) ) / (p.m0Sqr+12*p.lam0*v*v);
}


double findVEV(ParameterType p) {
  double v = 1.0;
  double vOLD = 0;
  while (fabs((v-vOLD)/v)>1E-10) {
    vOLD = v;
    double h = calcExpectationOfHiggsField(p, v);
    v -= h;    
  }
  return v;
}


double calcSelfEnergy(ParameterType p, double v, int exP0, int exP1, int exP2, int exP3) {
  double HiggsLoop = calcSingleHiggsLoop(p, v);
  double GoldstoneLoop = calcSingleGoldstoneLoop(p, v);
  double fermionDoubleLoop = calcDoubleFermionLoop(p, v, 0, exP0, exP1, exP2, exP3);
  double HiggsDoubleLoop = calcDoubleHiggsLoop(p, v, exP0, exP1, exP2, exP3);
  double GoldstoneDoubleLoop = calcDoubleHiggsLoop(p, v, exP0, exP1, exP2, exP3);
 
  return -p.lam0*(12*HiggsLoop + 12*GoldstoneLoop) - fermionDoubleLoop 
         + 16*v*v*p.lam0*p.lam0*(18*HiggsDoubleLoop + 6*GoldstoneDoubleLoop);
}


void plotSelfEnergy(ParameterType p, char* fileName, int scanMode) {
  FILE* file = fopen(fileName, "w");
  
  SimulationParameterSet simPara(p.m0Sqr, p.lam0, p.yt0, p.Nf, SimulationParameterSet_ContinuumNotation);
  double vev = findVEV(p);

  if (scanMode==1) {
    for (int I=0; I<p.L3/2; I++) {
      int exP0 = (I*p.L0)/p.L3;
      int exP1 = (I*p.L1)/p.L3;
      int exP2 = (I*p.L2)/p.L3;
      int exP3 = I;

      double pHatSqr = 0;  
      pHatSqr += 4*sqr(sin(exP0*pi/p.L0));
      pHatSqr += 4*sqr(sin(exP1*pi/p.L1));
      pHatSqr += 4*sqr(sin(exP2*pi/p.L2));
      pHatSqr += 4*sqr(sin(exP3*pi/p.L3));

      double selfEnergy = calcSelfEnergy(p, vev, exP0, exP1, exP2, exP3);

      fprintf(file," %1.15e  %1.15e %1.15e %1.15e %d %d %d %d %e %e %e %d %e %e %1.5e %1.5e %1.5e %1.5e %1.5e\n", pHatSqr, selfEnergy, vev, p.m0Sqr, simPara.getKappaN(), p.L0,p.L1,p.L2,p.L3, p.yt0,p.yb0,p.lam0,p.Nf,p.rho,p.r, simPara.getYN(), simPara.getLambdaN(), p.lam6, p.lam8, p.lam10);
    }
  } else if (scanMode==0) {
    for (int I=0; I<p.L3/2; I++) {
      int exP0 = 0;
      int exP1 = 0;
      int exP2 = 0;
      int exP3 = I;

      double pHatSqr = 0;  
      pHatSqr += 4*sqr(sin(exP0*pi/p.L0));
      pHatSqr += 4*sqr(sin(exP1*pi/p.L1));
      pHatSqr += 4*sqr(sin(exP2*pi/p.L2));
      pHatSqr += 4*sqr(sin(exP3*pi/p.L3));

      double selfEnergy = calcSelfEnergy(p, vev, exP0, exP1, exP2, exP3);

      fprintf(file," %1.15e  %1.15e %1.15e %1.15e %d %d %d %d %e %e %e %d %e %e %1.5e %1.5e %1.5e %1.5e %1.5e\n", pHatSqr, selfEnergy, vev, p.m0Sqr, simPara.getKappaN(), p.L0,p.L1,p.L2,p.L3, p.yt0,p.yb0,p.lam0,p.Nf,p.rho,p.r, simPara.getYN(), simPara.getLambdaN(), p.lam6, p.lam8, p.lam10);
    }
  } else if (scanMode==2) {
    double results[100000];
    int resCount = 0;
  
    for (int I=0; I<p.L3/2; I++) {
      int exP0 = 0;
      int exP1 = 0;
      int exP2 = 0;
      int exP3 = I;

      double pHatSqr = 0;  
      pHatSqr += 4*sqr(sin(exP0*pi/p.L0));
      pHatSqr += 4*sqr(sin(exP1*pi/p.L1));
      pHatSqr += 4*sqr(sin(exP2*pi/p.L2));
      pHatSqr += 4*sqr(sin(exP3*pi/p.L3));

      double selfEnergy = calcSelfEnergy(p, vev, exP0, exP1, exP2, exP3);
      results[2*resCount+0] = pHatSqr;
      results[2*resCount+1] = selfEnergy;
      resCount++;
    }
    for (int I=0; I<p.L3/2; I++) {
      int exP0 = (I*p.L0)/p.L3;
      int exP1 = (I*p.L1)/p.L3;
      int exP2 = (I*p.L2)/p.L3;
      int exP3 = I;

      double pHatSqr = 0;  
      pHatSqr += 4*sqr(sin(exP0*pi/p.L0));
      pHatSqr += 4*sqr(sin(exP1*pi/p.L1));
      pHatSqr += 4*sqr(sin(exP2*pi/p.L2));
      pHatSqr += 4*sqr(sin(exP3*pi/p.L3));


      bool found = false;
      for (int I=0; I<resCount; I++) {
        if (fabs(pHatSqr-results[2*I])<1E-10) found = true;
      }
      if (!found) {
        double selfEnergy = calcSelfEnergy(p, vev, exP0, exP1, exP2, exP3);
        results[2*resCount+0] = pHatSqr;
        results[2*resCount+1] = selfEnergy;
        resCount++;
      }
    }
  
    for (int exP0=0; exP0<p.L0; exP0++) {
      for (int exP1=0; exP1<p.L1; exP1++) {
        for (int exP2=0; exP2<p.L2; exP2++) {
          for (int exP3=0; exP3<p.L3; exP3++) {
	    if ((exP0<=p.L0/2) && (exP1<=p.L1/2) && (exP2<=p.L2/2) && (exP3<=p.L3/2)) {
              double pHatSqr = 0;  
              pHatSqr += 4*sqr(sin(exP0*pi/p.L0));
              pHatSqr += 4*sqr(sin(exP1*pi/p.L1));
              pHatSqr += 4*sqr(sin(exP2*pi/p.L2));
              pHatSqr += 4*sqr(sin(exP3*pi/p.L3));
	    
  	      if (pHatSqr<=1.1) {
	        bool found = false;
                for (int I=0; I<resCount; I++) {
		  if (fabs(pHatSqr-results[2*I])<1E-10) found = true;
		}
		
		if (!found) {
                  double selfEnergy = calcSelfEnergy(p, vev, exP0, exP1, exP2, exP3);
                  results[2*resCount+0] = pHatSqr;
                  results[2*resCount+1] = selfEnergy;
                  resCount++;
		}
	      }
	    }
	  }
	}
      }
    }
 
    //Sort Data
    bool change = true;
    while (change) {
      change = false;
      for (int I=0; I<resCount-1; I++) {
        if (results[2*I+0]>results[2*I+2]) {
	  double dummy1 = results[2*I+0];
	  double dummy2 = results[2*I+1];
	  results[2*I+0] = results[2*I+2];
	  results[2*I+1] = results[2*I+3];
	  results[2*I+2] = dummy1;
	  results[2*I+3] = dummy2;	  
	  change = true;
	}      
      }    
    }

    for (int I=0; I<resCount; I++) {
      fprintf(file," %1.15e  %1.15e %1.15e %1.15e %d %d %d %d %e %e %e %d %e %e %1.5e %1.5e %1.5e %1.5e %1.5e\n", results[2*I+0], results[2*I+1], vev, p.m0Sqr, simPara.getKappaN(), p.L0,p.L1,p.L2,p.L3, p.yt0,p.yb0,p.lam0,p.Nf,p.rho,p.r, simPara.getYN(), simPara.getLambdaN(), p.lam6, p.lam8, p.lam10);
    }
  }

  fclose(file);
}



void appendVEVToFile(ParameterType p, char* fileName) {
  FILE* file = fopen(fileName, "a");
  
  SimulationParameterSet simPara(p.m0Sqr, p.lam0, p.yt0, p.Nf, SimulationParameterSet_ContinuumNotation);
  double vev = findVEV(p);

  fprintf(file,"%1.15e %1.15e %1.15e %d %d %d %d %e %e %e %d %e %e %1.5e %1.5e %1.5e %1.5e %1.5e\n", vev, p.m0Sqr, simPara.getKappaN(), p.L0,p.L1,p.L2,p.L3, p.yt0,p.yb0,p.lam0,p.Nf,p.rho,p.r, simPara.getYN(), simPara.getLambdaN(), p.lam6, p.lam8, p.lam10);
  fclose(file);
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


void calcFermionMasses(ParameterType p, double vev, double mHiggsSqr, double mGoldSqr, double &topMass, double &bottomMass, const char* fileNameCorr, const char* fileNameEffMa) {
  topMass = 0;
  bottomMass = 0;  
  if ((p.yt0==0) && (p.yb0==0)) return;
  printf("Calculating Fermion Masses with Higgs mass %f and Goldstone mass %f and vev=%f\n", mHiggsSqr, mGoldSqr, vev);

  if ((p.L0>=0) && (p.L1>=0) && (p.L2>=0) && (p.L3>=0)) {
    makePSqrAvailableForParameters(p);
    makeDiracOpEigenvaluesAndEigenvectorsAvailable(p);
    
    double m0Sqr = p.m0Sqr;
    double yt = p.yt0;
    double yb = p.yb0;  
    double v = vev;
    double fac = 0.5 / p.rho;
    double facProj = 1.0 / p.rho;
    
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
      vector4D calcFermionMassesExternalMomentumP;
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

      Complex ewIDm1oRhoD = (ComplexUnity-facProj*ewP);
      ComplexMatrix IDm1oRhoD = ewIDm1oRhoD*mP1 + ewIDm1oRhoD*mP2 + adj(ewIDm1oRhoD)*mP3 + adj(ewIDm1oRhoD)*mP4;
      ComplexMatrix GammaHatP = Gamma5 * IDm1oRhoD;
      ComplexMatrix ProjMinusSmall = getProjectorMatrix(-1);
      ComplexMatrix ProjHatMinusSmall = 0.5*(ID4x4 - GammaHatP);

      ComplexMatrix ProjMinus(8);
      ProjMinus.setZero();
      ProjMinus.insertMatrix(ProjMinusSmall, 0, 0);
      ProjMinus.insertMatrix(ProjMinusSmall, 4, 4);
      
      ComplexMatrix ProjHatMinus(8);
      ProjHatMinus.setZero();
      ProjHatMinus.insertMatrix(ProjHatMinusSmall, 0, 0);
      ProjHatMinus.insertMatrix(ProjHatMinusSmall, 4, 4);

	      
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
              double GoldstoneProp = 1.0 / (pSqr[PSqrIndex] + mGoldSqr);    

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
      
      //Multiplications with Projectors to get Left-Righthanded structure
      topProp = ProjHatMinusSmall * topProp * ProjMinusSmall;
      bottomProp = ProjHatMinusSmall * bottomProp * ProjMinusSmall;
      
      
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

    printf("TopMass: %f\n", topMass);    
    printf("BottomMass: %f\n", bottomMass);    
    FILE* file = fopen(fileNameCorr, "w");
    for (int I=0; I<p.L3+1; I++) {
      fprintf(file, "%e %e %e %e %e %e\n", (double)I, topCorr[I], bottomCorr[I], v, m0Sqr, p.lam0);
    }
    fclose(file);
    file = fopen(fileNameEffMa, "w");
    for (int I=0; I<p.L3; I++) {
      fprintf(file, "%e %e %e %e %e %e %e %e\n", I+0.5, topEffMasses[I], yt*v, bottomEffMasses[I], yb*v, v, m0Sqr, p.lam0);
    }
    fclose(file);
 
    delete[] topEffMasses;
    delete[] topCorr;
    delete[] bottomEffMasses;
    delete[] bottomCorr;
    delete[] topPropTres;
    delete[] bottomPropTres;
  }
}


void appendTopMassToFile(ParameterType p, char* fileName, bool renMode) {
  FILE* file = fopen(fileName, "a");
  
  SimulationParameterSet simPara(p.m0Sqr, p.lam0, p.yt0, p.Nf, SimulationParameterSet_ContinuumNotation);
  double vev = findVEV(p);
  
  double topMass = NaN;
  double &bottomMass = NaN;
  double mHiggsSqr = p.m0Sqr + 12*p.lam0*vev*vev;
  double mGoldSqr = p.m0Sqr + 4*p.lam0*vev*vev; 
  
  if (renMode) {
    double selfEnergyAtP0 = calcSelfEnergy(p, vev, 0, 0, 0, 0);
    mHiggsSqr = p.m0Sqr + 12*p.lam0*vev*vev - selfEnergyAtP0;
    mGoldSqr = 0; 
  }
  
  calcFermionMasses(p, vev, mHiggsSqr, mGoldSqr, topMass, bottomMass, "corr.dat", "effMa.dat");

  fprintf(file,"%1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %d %d %d %d %e %e %e %d %e %e %1.5e %1.5e %1.5e %1.5e %1.5e\n", 
                vev, p.m0Sqr, mHiggsSqr, mGoldSqr, bottomMass, topMass,
                simPara.getKappaN(), p.L0,p.L1,p.L2,p.L3, p.yt0,p.yb0,p.lam0,p.Nf,p.rho,p.r, simPara.getYN(), simPara.getLambdaN(), p.lam6, p.lam8, p.lam10);
  fclose(file);
}


void plotFermionCorrelatorAndEffMasses(ParameterType p, char* fileNameCorr, char* fileNameEffMa, bool renMode) {
//  double vev = findVEV(p);

double vev = sqrt(2*0.12313) * 1.246077235374959e+00;
//double vev = sqrt(2*0.12301) * 4.196570372841069e-01;

  double topMass = NaN;
  double &bottomMass = NaN;
  double mHiggsSqr = p.m0Sqr + 12*p.lam0*vev*vev;
  double mGoldSqr = p.m0Sqr + 4*p.lam0*vev*vev; 
  
  if (renMode) {
    double selfEnergyAtP0 = calcSelfEnergy(p, vev, 0, 0, 0, 0);
    mHiggsSqr = p.m0Sqr + 12*p.lam0*vev*vev - selfEnergyAtP0;
    mGoldSqr = 0; 
  }
  
  calcFermionMasses(p, vev, mHiggsSqr, mGoldSqr, topMass, bottomMass, fileNameCorr, fileNameEffMa);
  printf("Top Mass: %f\n", topMass);  
  printf("Bottom Mass: %f\n", bottomMass);  
}




int main(int argc,char **argv) {
  iniTools(5517);
  Parameters.L0 = 12;
  Parameters.L1 = 12;
  Parameters.L2 = 12;
  Parameters.L3 = 32;

  Parameters.yt0 = 175.0 / Physical_VEV_GeV;
  Parameters.yb0 = 175.0 / Physical_VEV_GeV;
  Parameters.lam0 = 0.00;//0.0016462;
  Parameters.lam6 = 0.00;
  Parameters.lam8 = 0.00;
  Parameters.lam10 = 0.0;
//kappa=0.12301 lam=0  
//  Parameters.m0Sqr = 1.294203723274534e-01;
//kappa=0.12313 lam=0
//  Parameters.m0Sqr = 1.214976041582066e-01;//   0.03750; //0.11328;//0.12150;//0.11302;

//kappa=0.12313 lam=0.001
  Parameters.m0Sqr = 1.204999995949116e-01;

  Parameters.Nf = 1;
  Parameters.rho = 1.0;
  Parameters.r = 0.5; 
  diracOp = new NeubergerMatrix(Parameters.rho, Parameters.r, 1, 1, 1, 1, 2);
  QuietMode = false;



double xxx = calcSingleHiggsLoop(Parameters, 0.6);
double xxx2 = calcSingleGoldstoneLoop(Parameters, 0.6);
double xxx3 = calcSingleFermionLoop(Parameters, 0.6);
double xxx4 = calcDoubleFermionLoop(Parameters, 0.6, 0, 0, 0, 0, 0);
double xxx5 = calcDoubleHiggsLoop(Parameters, 0.6, 0, 0, 0, 0);
double xxx6 = calcDoubleGoldstoneLoop(Parameters, 0.6, 0, 0, 0, 0);
double xxx8 = findVEV(Parameters);
double xxx7 = calcExpectationOfHiggsField(Parameters, xxx8);
double xxx9 = 12*Parameters.lam0*sqr(xxx8)-calcSelfEnergy(Parameters, xxx8, 0,0,0,0);

printf("%f %f %f %f %f %f %f %f %f\n", xxx, xxx2, xxx3, xxx4, xxx5, xxx6, xxx7, xxx8, xxx9);



  char* fileName = constructFileName(Parameters, "M0Scan");
  startFile(fileName);
  for (Parameters.m0Sqr=0.115; Parameters.m0Sqr<=0.1154; Parameters.m0Sqr+=0.0002) {
    appendVEVToFile(Parameters, fileName);
  }
  delete[] fileName;

  fileName = constructFileName(Parameters, "TopMassScan");
  startFile(fileName);
  for (Parameters.m0Sqr=0.1246; Parameters.m0Sqr<=0.1250; Parameters.m0Sqr+=0.0002) {
    appendTopMassToFile(Parameters, fileName, true);
  }
  delete[] fileName;

  char* fileNameCorr = constructFileName(Parameters, "FermionCorrelator");
  char* fileNameEffMa = constructFileName(Parameters, "EffectiveMasses");
 plotFermionCorrelatorAndEffMasses(Parameters, fileNameCorr, fileNameEffMa, false);
  delete[] fileNameCorr;
  delete[] fileNameEffMa;

  fileNameCorr = constructFileName(Parameters, "FermionCorrelatorRenMode");
  fileNameEffMa = constructFileName(Parameters, "EffectiveMassesRenMode");
  plotFermionCorrelatorAndEffMasses(Parameters, fileNameCorr, fileNameEffMa, true);
  delete[] fileNameCorr;
  delete[] fileNameEffMa;

  fileName = constructFileName(Parameters, "SelfEnergyL3Scan");
  plotSelfEnergy(Parameters, fileName, 0);
  delete[] fileName;

  fileName = constructFileName(Parameters, "SelfEnergyDiagScan");
  plotSelfEnergy(Parameters, fileName, 1);
  delete[] fileName;

  fileName = constructFileName(Parameters, "SelfEnergyCombScan");
  plotSelfEnergy(Parameters, fileName, 2);
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
