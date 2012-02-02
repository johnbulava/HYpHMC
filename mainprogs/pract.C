#include "EvaluateObservablePropagatorBase.h"
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
#include <cln/number.h>
#include <cln/float.h>
#include <cln/complex.h>

//#include "StateDescriptorReader.h"

//#include "PTAnalysis.h"


NeubergerMatrix* diracOp = NULL;
ComplexVector v[4];
double Yukawa0;
double HiggsVEV;
FermionMatrixOperations* fermiOps = NULL;


double func(double p, double para) {
  return 1.0/p;

}



HighPrecisionComplexPolynom getAnalyticalContinuationOfPropagatoFit(double A1, double A2, double A3, double A4, int N, int DIGIT) {
  char* xxxStr = new char[1000];
  snprintf(xxxStr, 1000,"0.0e+0_%d",DIGIT);
  cln::cl_F ZERO = xxxStr;
  snprintf(xxxStr, 1000,"1.0e+0_%d",DIGIT);
  cln::cl_F ONE = xxxStr;
  snprintf(xxxStr, 1000,"2.0e+0_%d",DIGIT);
  cln::cl_F TWO = xxxStr;
  snprintf(xxxStr, 1000,"1.0e+0_%d",DIGIT);  
  cln::cl_F fac = xxxStr;
  
  HighPrecisionComplexPolynom StartNenner(3, DIGIT);
  StartNenner.setCoeff(0, Complex(A4,0));
  StartNenner.setCoeff(2, Complex(1,0));
  HighPrecisionComplexPolynom Nenner = StartNenner;
  
  HighPrecisionComplexPolynom TwoP(2, DIGIT);
  TwoP.setCoeff(1, complex(TWO, ZERO));
  HighPrecisionComplexPolynom FacPol(1, DIGIT);
  FacPol.setCoeff(0, complex(ONE, ZERO));
  
  HighPrecisionComplexPolynom Zaehler(3, DIGIT);
  Zaehler.setCoeff(2, Complex(A3/A1,0));
  
  HighPrecisionComplexPolynom anaCon(N+1,DIGIT);
  for (int I=0; I<=N; I++) {
    cln::cl_N n = Nenner.evaluateAt(ZERO);
    cln::cl_N z = Zaehler.evaluateAt(ZERO);
    
//Nenner.evaluateAt(Complex(0,0)).print();
//Zaehler.evaluateAt(Complex(0,0)).print();
    
    
    anaCon.setCoeff(I, (fac * z)/n);
//printf("sa\n");    
    
    snprintf(xxxStr, 1000,"%d.0e+0_%d",I+1,DIGIT);
    cln::cl_F Ip1 = xxxStr;
    
    fac = fac / Ip1;
    
//Zaehler.print();  
//Nenner.print();  
 
//((Zaehler.calcDerivative())*Nenner).print();    
    
    FacPol.setCoeff(0, complex(Ip1, ZERO));
    Zaehler = (Zaehler.calcDerivative())*StartNenner - (Zaehler*TwoP)*FacPol;
    Nenner = Nenner * StartNenner;
  }
  
  anaCon.addCoeff(0, Complex((A2*A2)/A1,0));
  anaCon.addCoeff(2, Complex(1.0/A1,0));
  
  delete[] xxxStr;
  return anaCon;
}


bool loadPhiconfiguration(char* fileName, double* phiField, double &weight, bool &weightAvail) {
  if (LogLevel>2) printf("Loading Configuration: %s...",fileName);  

double f = fabs(-1);
printf("%f\n",f);

  std::fstream confFile;
  confFile.open(fileName, std::ios::in);
  if (!confFile.good()) {
    if (LogLevel>2) printf("ERROR !!!\n");
    return false;
  }

  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();

  confFile.read((char*)phiField, 32*L0*L1*L2*L3);
  if (confFile.eof()) {
    if (LogLevel>2) printf("ERROR !!!\n");
    return false;
  }

  confFile.read((char*)(&weight),8);
  if (!confFile.eof()) {
    if (LogLevel>2) printf(" with weight: %f ",weight);
    weightAvail = true;
  } else {
    if (LogLevel>2) printf(" without weight ");
    weight = NaN;
    weightAvail = false;
  }

  confFile.close();

  if (LogLevel>2) printf("successfully.\n");
  return true;
} 



int main(int argc,char **argv) {
  iniTools(5517);



  double m0 = -0.005;
  double Z = 0.981490;
  int N=2;
  double coeff[4];
  coeff[0] = 0.007277;
  coeff[1] = -0.130292;
  coeff[2] = 0.002426;
  coeff[3] = 0.044063;



/*Complex c1(0,-2);
Complex c2(3,4);

HighPrecisionComplex ccc(1000, c1);
HighPrecisionComplex ccc2(1000, c2);

ccc.print();
ccc2.print();

(ccc+ccc2).print();

(c1/c2).print();
(ccc/ccc2).print();

(exp(c1)).print();
(exp(ccc)).print();

(log(c1)).print();
(log(ccc)).print();

exit(0);*/





Complex p(0.0000,0.0046);

Complex c1 = (calcBosonic1LoopInvPropagatorFit(p, m0, Z, N, coeff));
Complex c2 = (calcBosonic1LoopInvPropagatorFitOnSecondSheet(p, m0, Z, N, coeff));

printf("%1.15f %1.15f\n",c1.x, c1.y);
printf("%1.15f %1.15f\n",c2.x, c2.y);


//findZeroOfBosonic1LoopInvPropagatorOnSecondSheetFit(pole, poleVal, m0, Z, N,  &(coeff[0]), 0.102344); 

//  m0 = 0.06;
/*calcGoldstone1LoopInvPropagator(p, m0, m0, 0.8, Z, 0.456).print();
Complex c = calcGoldstone1LoopInvPropagatorWithP0PartSubtracted(p, m0, m0, 0.8, Z, 0.456);
printf("%1.15f\n", c.x);
calcGoldstone1LoopInvPropagatorHighPrecision(p, m0, m0, 0.8, Z, 0.456).print();
Complex c2 = calcGoldstone1LoopInvPropagatorWithP0PartSubtractedHighPrecision(p, m0, m0, 0.8, Z, 0.456);
printf("%1.15f\n", c2.x);
*/



//calcBosonic1LoopContribution(p, m0).print();
//calcBosonic1LoopContributionOnSecondSheet(p, m0).print();

//calcBosonic1LoopInvPropagatorFit(p, m0, Z, N, coeff).print();
//calcBosonic1LoopInvPropagatorFitOnSecondSheet(p, m0, Z, N, coeff).print();



/*findZeroOfBosonic1LoopInvPropagatorFit(pole, poleVal, m0, Z, N, &(coeff[0]), 1.2*m0);
findZeroOfBosonic1LoopInvPropagatorOnSecondSheetFit(pole, poleVal, m0 - (288+96)*0.000139, Z, N, &(coeff[0]), 1.2*m0);
*/

//pole.print();
//poleVal.print();


exit(0);

 






/*  std::fstream fileX;
  fileX.open(argv[1], std::ios::in);
  for (int I=0; I<120; I++) {
    Complex c;
    fileX.read((char*)(&c.x),8);
    fileX.read((char*)(&c.y),8);
    c.print();
  }
  exit(0);*/



  DebugMode = true;
  LogLevel=5;



  fermiOps = new FermionMatrixOperations(16, 16, 16, 32, 1.0, 0.5, 0.35302);
  double* phiField = (double*)createSuperAlignedComplex(16*16*16*32*4);
  double weight = NaN;
  bool weightAvail = false;
  
  char* fileName = "dataBase/data/results/pHMC/configurations/subFolderL16x16x16x32Nf1Kap0.12313Lam0.00000Y0.35302Rho1.000R0.500PolDeg18PolAl0.500_level12/PhiConfL16x16x16x32Nf1Kap0.12313Lam0.00000Y0.35302Rho1.000R0.500PolDeg18PolAl0.500_level12_32.dat";
  bool b = loadPhiconfiguration(fileName, phiField, weight, weightAvail);
  printf("Load success: %d\n",b);

  vector4D m;
  m[0] = 0;
  m[1] = 0;
  m[2] = 0;
  m[3] = 0;
  double max = 0;
  for (int I=0; I<16*16*16*32; I++) {
    m[0] += phiField[4*I+0];
    m[1] += phiField[4*I+1];
    m[2] += phiField[4*I+2];
    m[3] += phiField[4*I+3];    
    if (fabs(phiField[4*I+0]) > max) max = fabs(phiField[4*I+0]);
    if (fabs(phiField[4*I+1]) > max) max = fabs(phiField[4*I+1]);
    if (fabs(phiField[4*I+2]) > max) max = fabs(phiField[4*I+2]);
    if (fabs(phiField[4*I+3]) > max) max = fabs(phiField[4*I+3]);
    
  }
  double mag = sqrt(sqr(m[0])+sqr(m[1])+sqr(m[2])+sqr(m[3])) / (16*16*16*32);
  printf("mag: %f\n", mag);
  m[0] = 0;
  m[1] = 0;
  m[2] = 0;
  m[3] = 0;
  int count = 0;
  for (int I1=0; I1<16; I1++) {
    for (int I2=0; I2<16; I2++) {
      for (int I3=0; I3<16; I3++) {
        for (int I4=0; I4<32; I4++) {
	  double sign = 1;
	  if (((I1+I2+I3+I4) % 2) == 1) sign = -1;
          m[0] += sign*phiField[4*count+0];
          m[1] += sign*phiField[4*count+1];
          m[2] += sign*phiField[4*count+2];
          m[3] += sign*phiField[4*count+3];   
	  count++;
	}
      }
    }
  }
  double smag = sqrt(sqr(m[0])+sqr(m[1])+sqr(m[2])+sqr(m[3])) / (16*16*16*32);
  printf("5 0.01 smag: %f\n", smag);
  printf("max: %f\n", max);
  fermiOps->setPreconditioner(true, mag,  smag);

//  fermiOps->exactFermionMatrixConditionNumber(phiField, eigMin, eigMax, invCond, false, 1, false);
  
  
Complex* eigenV = fermiOps->calcFermionMatrixARPACKEigenValues(5, 2, phiField,
0.01, false, NULL, false, true);
for (int I=0; I<2; I++) {
  eigenV[I].print();
}
  
  


exit(0);


/*Complex p(+0.0001,0.1);
double mH = 0.2;
double mG=0.01;
double lamRen=0.0005;
double vren=1.00;
int n=1;
double E=0.197;
//(-1*calcBosonic1LoopContribution(p, m0)).print();
//(-1*calcBosonic1LoopContributionOnSecondSheet(p, m0)).print();

Complex HiggsPole = findPoleOfBosonic1LoopPropagatorFromRenPT(mH, mG, lamRen, vren, n);
HiggsPole.print();

Complex invProp = calcBosonic1LoopInvPropagatorFromRenPT(p, HiggsPole, mG, lamRen, vren, n);

invProp.print();

double RHO = calcSpectralFunctionOfBosonic1LoopPropagatorFromRenPT(E, HiggsPole, mG, lamRen, vren, n);
printf("rho: %f\n", RHO);


FILE* fileS = fopen("SpectralFunction.dat", "w");
int N=10000;
for (int I=0; I<N; I++) {
  E = 0.75*I/N;
  double RHO = calcSpectralFunctionOfBosonic1LoopPropagatorFromRenPT(E, HiggsPole, mG, lamRen, vren, n);
  fprintf(fileS, "%1.15f %1.15f\n", E, RHO);
}
fclose(fileS);


exit(0);*/

/*calcBosonic1LoopContribution(Complex(0.000,0), 0.4).print();


FILE* file2=fopen("Correlator.dat", "w");
int L = 32;
double coeff[4];
coeff[0] = -0.20;
coeff[1] = 0.83;
coeff[2] = 0.006;
coeff[3] = 0.032;
for (int DeltaT=0; DeltaT<L; DeltaT++) {
  Complex corr = calcBosonic1LoopCorrelatorFit(DeltaT, L, 0.0, 0.998, 2, &(coeff[0]));
  fprintf(file2, "%d %1.15f %1.15f\n", DeltaT, corr.x, corr.y);
}
fclose(file2);

Complex pole;
Complex poleVal;
bool b = findZeroOfBosonic1LoopInvPropagatorFit(pole, poleVal, 0.017, 0.979, 2, &(coeff[0]), 0.4);
printf("%d\n", b);
pole.print();
poleVal.print();




exit(0);*/
  






/*Complex xxx(10.7,-0.00001);
Complex resxxx = arctanhAsPowerSeries(xxx);
printf("%1.15f %1.15f\n",resxxx.x,resxxx.y);

resxxx = arctanh(xxx);
printf("%1.15f %1.15f\n",resxxx.x,resxxx.y);

Complex p(4.01, 0);
resxxx = calcBosonic1LoopContribution(p, 0.032);
printf("%1.15f %1.15f\n",resxxx.x,resxxx.y);




exit(0);*/





/*ComplexPolynom poly(3);
ComplexPolynom poly2(2);
HighPrecisionComplexPolynom hpoly(3,40);
HighPrecisionComplexPolynom hpoly2(2,40);

poly.setCoeff(0,Complex(1,2));
poly.setCoeff(1,Complex(3,4));
poly.setCoeff(2,Complex(5,6));

poly2.setCoeff(0,Complex(11,12));
poly2.setCoeff(1,Complex(13,14));

hpoly.setCoeff(0,Complex(+1,0));
hpoly.setCoeff(1,Complex(0,0));
hpoly.setCoeff(2,Complex(2,0));

hpoly2.setCoeff(0,Complex(11,12));
hpoly2.setCoeff(1,Complex(13,14));



poly.print();
poly2.print();
printf("\n");
hpoly.print();
hpoly2.print();
printf("\n");

printf("%d\n", hpoly.getOrder());

Complex* roots = hpoly.getRoots();
for (int I=0; I<hpoly.getOrder(); I++) {
roots[I].print();
}*/

double A1 = 1;
double A2 = 0.2;
double A3 = 0.1;
double A4 = 1.5;





double x = 0.5;

HighPrecisionComplexPolynom anaCon(2,1000);

for (int N=2; N<30; N+=2) {
anaCon = getAnalyticalContinuationOfPropagatoFit(A1, A2, A3, A4, N, 1000);

Complex res = anaCon.evaluateAt(Complex(x,0));
double exact = (x*x+A2*A2+A3*x*x/(x*x+A4))/A1;
printf("exact: %1.15f\n", exact);

printf("%1.15f\n", res.x);
printf("diff: %1.15f\n\n", res.x-exact);

}

FILE* file=fopen("RootsOfFit.dat","w");
Complex* roots = anaCon.getRoots();
for (int I=0; I<anaCon.getOrder(); I++) {
  fprintf(file, "%1.15e %1.15e\n", roots[I].x, roots[I].y);
}

fclose(file);



/*(poly.calcDerivative()).print();
printf("\n");
(hpoly.calcDerivative()).print();*/

/*

//double q = sqrt(0.472894247);
double q = sqrt(1.10);
double h = 1E-5;
  double res = calcLuescherPhiFunction(q-0.5*h, 1E-5);
  double res2 = calcLuescherPhiFunction(q+0.5*h, 1E-5);
  double d = calcLuescherPhiDerivative(q, 1E-5);
  
  printf("q=%1.5f num derive = %1.10f ana derive = %1.10f\n", q, (res2-res)/h , d);
*/
}

