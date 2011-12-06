#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <iostream>


#include "Complex.h"
#include "Quat.h"
#include "Global.h"
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include "Tools.h"
#include "NeubergerMatrix.h"
#include "WilsonMatrix.h"
#include "FermionMatrixOperations.h"


//Types
struct NeighboursType {
  int neighbourIndex[8];
};


//Variables
int Parameter_N = 0;
int Parameter_Nf = 0;
double Parameter_RHO = 0;
double Parameter_R = 0;
double Parameter_Lambda = 0;
double Parameter_Y = 0;
double Parameter_Kappa = 0;
double Parameter_PhiChangeRadius = 0;
double Parameter_TurnAngelMax = 0;
int Parameter_nhit = 0;
int Parameter_nSweeps = 0;
int Parameter_thermalSweeps = 0;
char* Parameter_filenamePrefix;
int Parameter_FLAG_CalcNeubergerDetWithXi = 0;
int Parameter_FLAG_CalcWilsonDetWithoutXi = 0;
int Parameter_FLAG_CalcNeubergerDetWithoutXi = 0;
int Parameter_FLAG_CalcWilsonDetWithXi = 0;
int Parameter_FLAG_ImprovedPhiDynamics = 0;
NeighboursType* neighbours = NULL;
vector4D* phiField = NULL;
FermionMatrixOperations* fermiOps = NULL;
ComplexMatrix calcMatrix(1);


// X-Werte duerfen negativ werden bis -Parameter_N. Nach oben keine Schranke.
#define LPOS(x0, x1, x2, x3) (((x0+Parameter_N) % Parameter_N)*Parameter_N*Parameter_N*Parameter_N + ((x1+Parameter_N) % Parameter_N)*Parameter_N*Parameter_N + ((x2+Parameter_N) % Parameter_N)*Parameter_N + ((x3+Parameter_N) % Parameter_N))


void loadParameters(int MultiProcessNr, bool ll) {
  FILE* file;
  if (LogLevel>0) printf("Trying to read parameter-file 'SimulationParameters.txt'...\n");
  file = fopen("SimulationParameters.txt","r");
  char *Comment = new char[500];
  bool error = false;
  Parameter_filenamePrefix = new char[500];

  if (fscanf(file,"%s",Parameter_filenamePrefix)==1) {
    if (LogLevel>0) printf(" -> Filename-Prefix: %s\n",Parameter_filenamePrefix);
  } else error = true;
  fgets(Comment, 500, file);
  
  if (fscanf(file,"%d",&Parameter_N)==1) {
    if (LogLevel>0) printf(" -> N: %d\n",Parameter_N);
  } else error = true;
  fgets(Comment, 500, file);
  
  if (fscanf(file,"%d",&Parameter_Nf)==1) {
    if (LogLevel>0) printf(" -> Nf: %d\n",Parameter_Nf);
  } else error = true;
  fgets(Comment, 500, file);
  
  if (fscanf(file,"%lf",&Parameter_RHO)==1) {
    if (LogLevel>0) printf(" -> Rho: %f\n",Parameter_RHO);
  } else error = true;
  fgets(Comment, 500, file);
  
  if (fscanf(file,"%lf",&Parameter_R)==1) {
    if (LogLevel>0) printf(" -> R: %f\n",Parameter_R);
  } else error = true;
  fgets(Comment, 500, file);
  
  if (fscanf(file,"%lf",&Parameter_Lambda)==1) {
    if (LogLevel>0) printf(" -> Lambda: %f\n",Parameter_Lambda);
  } else error = true;
  fgets(Comment, 500, file);

  double yMin, yMax;
  int yTics;
  if (fscanf(file,"%lf",&yMin)==1) {
    if (LogLevel>0) printf(" -> yMin: %f\n",yMin);
  } else error = true;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf",&yMax)==1) {
    if (LogLevel>0) printf(" -> yMax: %f\n",yMax);
  } else error = true;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d",&yTics)==1) {
    if (LogLevel>0) printf(" -> yTics: %d\n",yTics);
  } else error = true;
  fgets(Comment, 500, file);
  
  double kappaMin, kappaMax;
  int kappaTics;
  if (fscanf(file,"%lf",&kappaMin)==1) {
    if (LogLevel>0) printf(" -> kappaMin: %f\n",kappaMin);
  } else error = true;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf",&kappaMax)==1) {
    if (LogLevel>0) printf(" -> kappaMax: %f\n",kappaMax);
  } else error = true;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d",&kappaTics)==1) {
    if (LogLevel>0) printf(" -> kappaMTics: %d\n",kappaTics);
  } else error = true;
  fgets(Comment, 500, file);
  
  MultiProcessNr = MultiProcessNr % (yTics*kappaTics);
  if (kappaTics <= 1) {
    Parameter_Kappa = kappaMin;
  } else {
    Parameter_Kappa = kappaMin+(kappaMax-kappaMin)*(MultiProcessNr % kappaTics)/(kappaTics-1);
  }
  if (yTics <= 1) {
    Parameter_Y = yMin;
  } else {
    Parameter_Y = yMin+(yMax-yMin)*(MultiProcessNr / kappaTics)/(yTics-1);
  }
  if (LogLevel>0) printf(" -> Y (selected): %f\n",Parameter_Y);
  if (LogLevel>0) printf(" -> Kappa (selected): %f\n",Parameter_Kappa);
  
  if (fscanf(file,"%lf",&Parameter_PhiChangeRadius)==1) {
    if (LogLevel>0) printf(" -> Phi change radius: %f\n",Parameter_PhiChangeRadius);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%lf",&Parameter_TurnAngelMax)==1) {
    if (LogLevel>0) printf(" -> Phi change angel max: %f\n",Parameter_TurnAngelMax);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_nhit)==1) {
    if (LogLevel>0) printf(" -> Number of Phi hits per sweep: %d\n",Parameter_nhit);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_FLAG_CalcNeubergerDetWithXi)==1) {
    if (LogLevel>0) printf(" -> Flag - Calculate Neuberger Determinant with Xi-fields: %d\n",Parameter_FLAG_CalcNeubergerDetWithXi);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_FLAG_CalcWilsonDetWithoutXi)==1) {
    if (LogLevel>0) printf(" -> Flag - Calculate Wilson Determinant without Xi-fields: %d\n",Parameter_FLAG_CalcWilsonDetWithoutXi);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_FLAG_CalcNeubergerDetWithoutXi)==1) {
    if (LogLevel>0) printf(" -> Flag - Calculate Neuberger Determinant without Xi-fields: %d\n",Parameter_FLAG_CalcNeubergerDetWithoutXi);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_FLAG_CalcWilsonDetWithXi)==1) {
    if (LogLevel>0) printf(" -> Flag - Calculate Wilson Determinant with Xi-fields: %d\n",Parameter_FLAG_CalcWilsonDetWithXi);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_FLAG_ImprovedPhiDynamics)==1) {
    if (LogLevel>0) printf(" -> Flag - Use improved Phi-Dynamics: %d\n",Parameter_FLAG_ImprovedPhiDynamics);
  } else error = true;
  fgets(Comment, 500, file);

  if ((Parameter_FLAG_ImprovedPhiDynamics) && (Parameter_Y==0)) {
    printf("Improved Phi Dynamics can not be used for yN = 0!!!\n");
    exit(0);
  }

  if (fscanf(file,"%d",&Parameter_nSweeps)==1) {
    if (LogLevel>0) printf(" -> Number of sweeps before measurement: %d\n",Parameter_nSweeps);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_thermalSweeps)==1) {
    if (LogLevel>0) printf(" -> Number of thermalizing sweeps: %d\n",Parameter_thermalSweeps);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&LogLevel)==1) {
    if ((LogLevel>0) && (ll)) printf(" -> Log - Level: %d\n",LogLevel);
  } else error = true;
  fgets(Comment, 500, file);

  fclose(file);
  if (error) {
    printf("... ERROR!\n");
    exit(0);
  }  
  delete [] Comment;
  if ((LogLevel>0) && (ll)) printf("...sucessfully.\n");
}


void iniFields() {
  int x0,x1,x2,x3,I;
  int count = 0;
  
  if (LogLevel>0) printf("Initializing fields...");
  neighbours = new NeighboursType[Parameter_N*Parameter_N*Parameter_N*Parameter_N];
  
  for (x0=0; x0<Parameter_N; x0++) {
    for (x1=0; x1<Parameter_N; x1++) {
      for (x2=0; x2<Parameter_N; x2++) {
        for (x3=0; x3<Parameter_N; x3++) {
	  neighbours[count].neighbourIndex[0] = LPOS(x0+1,x1,x2,x3);
	  neighbours[count].neighbourIndex[1] = LPOS(x0-1,x1,x2,x3);
	  neighbours[count].neighbourIndex[2] = LPOS(x0,x1+1,x2,x3);
	  neighbours[count].neighbourIndex[3] = LPOS(x0,x1-1,x2,x3);
	  neighbours[count].neighbourIndex[4] = LPOS(x0,x1,x2+1,x3);
	  neighbours[count].neighbourIndex[5] = LPOS(x0,x1,x2-1,x3);
	  neighbours[count].neighbourIndex[6] = LPOS(x0,x1,x2,x3+1);
	  neighbours[count].neighbourIndex[7] = LPOS(x0,x1,x2,x3-1);
	  count++;
	}
      }
    }
  }
  
  phiField = new vector4D[Parameter_N*Parameter_N*Parameter_N*Parameter_N];  
  if (!(Parameter_Lambda==Parameter_Lambda)) {
    //Lambda == Infinity (i.e. NaN)
    for (I=0; I<Parameter_N*Parameter_N*Parameter_N*Parameter_N; I++) {
      phiField[I][0] = sqrt(Parameter_Nf);
      phiField[I][1] = 0;
      phiField[I][2] = 0;
      phiField[I][3] = 0;
    }
  } else {
    //Lambda < Infinity
    for (I=0; I<Parameter_N*Parameter_N*Parameter_N*Parameter_N; I++) {
      phiField[I][0] = sqrt(Parameter_Nf);
      phiField[I][1] = 0;
      phiField[I][2] = 0;
      phiField[I][3] = 0;
    }
  }
  if (LogLevel>0) printf("sucessfully.\n");    
}


void desini() {
  if (LogLevel>0) printf("Desinitializing...");
  delete [] neighbours;
  if (LogLevel>0) printf("sucessfully.\n");      
}


double MonteCarloSweep() {
  int pos, mu, nu, n;
  double deltaPhi, deltaS, Phi;
  int changedCount = 0;
  bool changed;
  double neighbourSum;
  vector4D neighbourSumVec4D;
  double phiSqr;
  double Nf2Lam = 2*Parameter_Nf*Parameter_Lambda;
  double Lam2 = 2*Parameter_Lambda;
  double Kap2 = 2*Parameter_Kappa;
  double deltaPhiPhi, deltaPhideltaPhi;
  double acceptRate;
  vector4D phiSum4D;
  double phiSumNormSquared;
  double deltaPhiSumNormSquared;


  if (LogLevel>5) printf("Performing Monte Carlo Sweep...");
  
  phiSum4D[0] = 0;
  phiSum4D[1] = 0;
  phiSum4D[2] = 0;
  phiSum4D[3] = 0;
  for (pos=0; pos<Parameter_N*Parameter_N*Parameter_N*Parameter_N; pos++) {
    phiSum4D[0] += phiField[pos][0];
    phiSum4D[1] += phiField[pos][1];
    phiSum4D[2] += phiField[pos][2];
    phiSum4D[3] += phiField[pos][3];
  }
  phiSumNormSquared = phiSum4D[0] * phiSum4D[0] 
                    + phiSum4D[1] * phiSum4D[1] 
                    + phiSum4D[2] * phiSum4D[2] 
                    + phiSum4D[3] * phiSum4D[3];
  
  if (!(Parameter_Lambda==Parameter_Lambda)) {
    //Sweep for Lambda == Infinity (i.e. NaN here)
    if (Parameter_FLAG_ImprovedPhiDynamics) {
      printf("Monte Carlo Step for Improved Phi Dynamics not implemented yet!!!\n");
      exit(0);
    }
    for (pos=0; pos<Parameter_N*Parameter_N*Parameter_N*Parameter_N; pos++) {
  
      for (mu=0; mu<4; mu++) {
        neighbourSumVec4D[mu] = phiField[neighbours[pos].neighbourIndex[0]][mu] + phiField[neighbours[pos].neighbourIndex[1]][mu]
                              + phiField[neighbours[pos].neighbourIndex[2]][mu] + phiField[neighbours[pos].neighbourIndex[3]][mu]
                              + phiField[neighbours[pos].neighbourIndex[4]][mu] + phiField[neighbours[pos].neighbourIndex[5]][mu]
                              + phiField[neighbours[pos].neighbourIndex[6]][mu] + phiField[neighbours[pos].neighbourIndex[7]][mu];
      }
      for (mu=0; mu<4; mu++) {
        for (nu=mu+1; nu<4; nu++) {
          changed = false;
	  
	  double SOld = phiField[pos][mu]*neighbourSumVec4D[mu] + phiField[pos][nu]*neighbourSumVec4D[nu];
          for (n=0; n<Parameter_nhit; n++) {
  	    double turnAngel = 2*(AdvancedZufall(AdvancedSeed)-0.5)*Parameter_TurnAngelMax;
	    
	    double s = sin(turnAngel);
	    double c = cos(turnAngel);
	    
	    double Tmu = c*phiField[pos][mu] + s*phiField[pos][nu];
	    double Tnu = c*phiField[pos][nu] - s*phiField[pos][mu];

  	    double SNew = Tmu*neighbourSumVec4D[mu] + Tnu*neighbourSumVec4D[nu];
	    deltaS = Kap2*(SOld - SNew);  //wg. -Kappa in Wirkung
      	    if (AdvancedZufall(AdvancedSeed) < exp(-deltaS)) {
	      phiField[pos][mu] = Tmu;
	      phiField[pos][nu] = Tnu;
	      SOld = SNew;
	      changed = true;
	    }
	  }
          if (changed) changedCount++;
	}
      }
      double norm = phiField[pos][0]*phiField[pos][0]
                  + phiField[pos][1]*phiField[pos][1]
                  + phiField[pos][2]*phiField[pos][2]
                  + phiField[pos][3]*phiField[pos][3];
      norm = sqrt(norm) / Parameter_Nf;
      phiField[pos][0] /= norm;
      phiField[pos][1] /= norm;
      phiField[pos][2] /= norm;
      phiField[pos][3] /= norm;
    }
    acceptRate = (1.0*changedCount)/(6*Parameter_N*Parameter_N*Parameter_N*Parameter_N);
  } else {
    //Sweep for Lambda < Infinity
    for (pos=0; pos<Parameter_N*Parameter_N*Parameter_N*Parameter_N; pos++) {
      phiSqr = phiField[pos][0]*phiField[pos][0] 
             + phiField[pos][1]*phiField[pos][1] 
             + phiField[pos][2]*phiField[pos][2] 
             + phiField[pos][3]*phiField[pos][3];
  
      for (mu=0; mu<4; mu++) {
        changed = false;
        Phi = phiField[pos][mu];
        neighbourSum = phiField[neighbours[pos].neighbourIndex[0]][mu] + phiField[neighbours[pos].neighbourIndex[1]][mu]
                     + phiField[neighbours[pos].neighbourIndex[2]][mu] + phiField[neighbours[pos].neighbourIndex[3]][mu]
                     + phiField[neighbours[pos].neighbourIndex[4]][mu] + phiField[neighbours[pos].neighbourIndex[5]][mu]
                     + phiField[neighbours[pos].neighbourIndex[6]][mu] + phiField[neighbours[pos].neighbourIndex[7]][mu];
  		    
        for (n=0; n<Parameter_nhit; n++) {
          deltaPhi = Parameter_PhiChangeRadius*2*(AdvancedZufall(AdvancedSeed)-0.5);
	  deltaPhiPhi = deltaPhi * Phi;
  	  deltaPhideltaPhi = deltaPhi*deltaPhi;

          deltaPhiSumNormSquared = 2*deltaPhi*phiSum4D[mu] + deltaPhideltaPhi;
         
          deltaS = -Kap2*deltaPhi*neighbourSum
	         +  2*deltaPhiPhi*(1-Nf2Lam+Lam2*phiSqr+Lam2*deltaPhi*deltaPhi)
	         +  deltaPhideltaPhi*(1-Nf2Lam+Lam2*phiSqr)
	         +  Parameter_Lambda*(4*deltaPhiPhi*deltaPhiPhi + deltaPhideltaPhi*deltaPhideltaPhi);

          if (Parameter_FLAG_ImprovedPhiDynamics) {
	    deltaS -= 4 * Parameter_Nf * log(1 + (deltaPhiSumNormSquared/phiSumNormSquared));
	  }

    	  if (AdvancedZufall(AdvancedSeed) < exp(-deltaS)) {
	    phiSqr -= phiField[pos][mu]*phiField[pos][mu];
  	    phiField[pos][mu] += deltaPhi;
	    phiSqr += phiField[pos][mu]*phiField[pos][mu];	  
            Phi = phiField[pos][mu];
	    phiSumNormSquared += deltaPhiSumNormSquared;
	    phiSum4D[mu] += deltaPhi;
	    changed = true;
  	  }
        }
        if (changed) changedCount++;
      }
    }
    acceptRate = (1.0*changedCount)/(4*Parameter_N*Parameter_N*Parameter_N*Parameter_N);
  }
  if (LogLevel>5) printf("ready. Update quote: %1.2f\n",acceptRate);
  return acceptRate;
}



void calculateDeterminants(double& logDetNeubergerWithXiFull, double& logDetNeubergerWithXiNoZeros,double& logDetWilsonWithOutXiFull, double& logDetWilsonWithOutXiNoZeros,
                           double& logDetNeubergerWithOutXiFull, double& logDetNeubergerWithOutXiNoZeros, double& logDetWilsonWithXiFull, double& logDetWilsonWithXiNoZeros) {
  logDetNeubergerWithXiFull = 0.0;
  logDetNeubergerWithXiNoZeros = 0.0;  
  logDetWilsonWithOutXiFull = 0.0;
  logDetWilsonWithOutXiNoZeros = 0.0;  
  logDetNeubergerWithOutXiFull = 0.0;
  logDetNeubergerWithOutXiNoZeros = 0.0;
  logDetWilsonWithXiFull = 0.0;
  logDetWilsonWithXiNoZeros = 0.0;

  if (Parameter_Y == 0.0) return;

  if (Parameter_FLAG_CalcNeubergerDetWithXi) {
    fermiOps->constructNeubergerWithXiFermionMatrix(calcMatrix, phiField);
  
    if (LogLevel>2) printf("Calculating Neuberger Determinant with Xi-fields...\n");
    logDetNeubergerWithXiFull = Parameter_Nf*calcLogDetScaledAbsNorm(calcMatrix, 0, 3.2*Parameter_RHO, true);
    logDetNeubergerWithXiNoZeros = Parameter_Nf*calcLogDetScaledAbsNorm(calcMatrix, 8, 3.2*Parameter_RHO, false);
    if (LogLevel>3) printEigenvalues("TempNeubergerWithXi",calcMatrix,ComplexUnity);
    if (LogLevel>2) printf("ready.\n");
    
    double pow = logDetNeubergerWithXiFull / log(10);
    if (LogLevel>1) printf("Full (normed) Neuberger determinant with Xi-fields approx. 10^%1.2f\n",pow);
    pow = logDetNeubergerWithXiNoZeros / log(10);
    if (LogLevel>1) printf("Non-zero-mode (normed) Neuberger determinant with Xi-fields approx. 10^%1.2f\n",pow);
  } else {
    logDetNeubergerWithXiFull = NaN;
    logDetNeubergerWithXiNoZeros = NaN;  
  }
  
  if (Parameter_FLAG_CalcWilsonDetWithoutXi) {
    fermiOps->constructWilsonWithOutXiFermionMatrix(calcMatrix, phiField);
  
    if (LogLevel>2) printf("Calculating Wilson Determinant without Xi-fields...\n");
    logDetWilsonWithOutXiFull = Parameter_Nf*calcLogDetScaledAbsNorm(calcMatrix, 0, 6.5*Parameter_R, true);
    logDetWilsonWithOutXiNoZeros = Parameter_Nf*calcLogDetScaledAbsNorm(calcMatrix, 8, 6.5*Parameter_R, false);
    if (LogLevel>3) printEigenvalues("TempWilsonWithOutXi",calcMatrix,ComplexUnity);
    if (LogLevel>2) printf("ready.\n");

    double pow = logDetWilsonWithOutXiFull / log(10);
    if (LogLevel>1) printf("Full (normed) Wilson determinant without Xi-fields approx. 10^%1.2f\n",pow);
    pow = logDetWilsonWithOutXiNoZeros / log(10);
    if (LogLevel>1) printf("Non-zero-mode (normed) Wilson determinant without Xi-fields approx. 10^%1.2f\n",pow);
  } else {
    logDetWilsonWithOutXiFull = NaN;
    logDetWilsonWithOutXiNoZeros = NaN;  
  }

  if (Parameter_FLAG_CalcNeubergerDetWithoutXi) {
    fermiOps->constructNeubergerWithOutXiFermionMatrix(calcMatrix, phiField);
  
    if (LogLevel>2) printf("Calculating Neuberger Determinant without Xi-fields...\n");
    logDetNeubergerWithOutXiFull = Parameter_Nf*calcLogDetScaledAbsNorm(calcMatrix, 0, 6.5*Parameter_R, true);
    logDetNeubergerWithOutXiNoZeros = Parameter_Nf*calcLogDetScaledAbsNorm(calcMatrix, 8, 6.5*Parameter_R, false);
    if (LogLevel>3) printEigenvalues("TempNeubergerWithOutXi",calcMatrix,ComplexUnity);
    if (LogLevel>2) printf("ready.\n");
    
    double pow = logDetNeubergerWithOutXiFull / log(10);
    if (LogLevel>1) printf("Full (normed) Neuberger determinant without Xi-fields approx. 10^%1.2f\n",pow);
    pow = logDetNeubergerWithOutXiNoZeros / log(10);
    if (LogLevel>1) printf("Non-zero-mode (normed) Neuberger determinant without Xi-fields approx. 10^%1.2f\n",pow);
  } else {
    logDetNeubergerWithOutXiFull = NaN;
    logDetNeubergerWithOutXiNoZeros = NaN;  
  }
  
  if (Parameter_FLAG_CalcWilsonDetWithXi) {
    fermiOps->constructWilsonWithXiFermionMatrix(calcMatrix, phiField);
  
    if (LogLevel>2) printf("Calculating Wilson Determinant with Xi-fields...\n");
    logDetWilsonWithXiFull = Parameter_Nf*calcLogDetScaledAbsNorm(calcMatrix, 0, 3.2*Parameter_RHO, true);
    logDetWilsonWithXiNoZeros = Parameter_Nf*calcLogDetScaledAbsNorm(calcMatrix, 8, 3.2*Parameter_RHO, false);
    if (LogLevel>3) printEigenvalues("TempWilsonWithXi",calcMatrix,ComplexUnity);
    if (LogLevel>2) printf("ready.\n");

    double pow = logDetWilsonWithXiFull / log(10);
    if (LogLevel>1) printf("Full (normed) Wilson determinant with Xi-fields approx. 10^%1.2f\n",pow);
    pow = logDetWilsonWithXiNoZeros / log(10);
    if (LogLevel>1) printf("Non-zero-mode (normed) Wilson determinant with Xi-fields approx. 10^%1.2f\n",pow);
  } else {
    logDetWilsonWithXiFull = NaN;
    logDetWilsonWithXiNoZeros = NaN;  
  }
  
  
  logDetNeubergerWithOutXiFull = NaN;
  logDetNeubergerWithOutXiNoZeros = NaN;
  logDetWilsonWithXiFull = NaN;
  logDetWilsonWithXiNoZeros = NaN;
}


void measure() {
  int I1,I2,I3,I4;
  int count = 0;
  vector4D measurePhi, measureStaggeredPhi;
  double measurePhiNorm, measureStaggeredPhiNorm;
  double avgNorm, sigmaNorm;
  double dummy;
  
  measurePhi[0] = 0;
  measurePhi[1] = 0;
  measurePhi[2] = 0;
  measurePhi[3] = 0;
  measureStaggeredPhi[0] = 0;
  measureStaggeredPhi[1] = 0;
  measureStaggeredPhi[2] = 0;
  measureStaggeredPhi[3] = 0;
  avgNorm = 0;
  sigmaNorm = 0;
  
  count = 0;
  for (I1=0; I1<Parameter_N; I1++) {
    for (I2=0; I2<Parameter_N; I2++) {
      for (I3=0; I3<Parameter_N; I3++) {
        for (I4=0; I4<Parameter_N; I4++) {
          measurePhi[0] += phiField[count][0];
          measurePhi[1] += phiField[count][1];
          measurePhi[2] += phiField[count][2];
          measurePhi[3] += phiField[count][3];
	  
	  avgNorm += dummy = sqrt( phiField[count][0]*phiField[count][0] 
	                          +phiField[count][1]*phiField[count][1]
	                          +phiField[count][2]*phiField[count][2]
	                          +phiField[count][3]*phiField[count][3]);
          sigmaNorm += dummy*dummy;
	  
	  double stagFac = 1;
	  if (((I1+I2+I3+I4) % 2) == 1) stagFac = -1;
          measureStaggeredPhi[0] += stagFac*phiField[count][0];
          measureStaggeredPhi[1] += stagFac*phiField[count][1];
          measureStaggeredPhi[2] += stagFac*phiField[count][2];
          measureStaggeredPhi[3] += stagFac*phiField[count][3];
	  count++;
	}
      }
    }
  }
  double norm = Parameter_N*Parameter_N*Parameter_N*Parameter_N;
  measurePhi[0] = measurePhi[0]/norm;
  measurePhi[1] = measurePhi[1]/norm;
  measurePhi[2] = measurePhi[2]/norm;
  measurePhi[3] = measurePhi[3]/norm;
  measureStaggeredPhi[0] = measureStaggeredPhi[0]/norm;
  measureStaggeredPhi[1] = measureStaggeredPhi[1]/norm;
  measureStaggeredPhi[2] = measureStaggeredPhi[2]/norm;
  measureStaggeredPhi[3] = measureStaggeredPhi[3]/norm;
  
  avgNorm /= norm;
  sigmaNorm = sqrt(sigmaNorm/norm - avgNorm*avgNorm);
  
  measurePhiNorm = measurePhi[0]*measurePhi[0]
                 + measurePhi[1]*measurePhi[1]
                 + measurePhi[2]*measurePhi[2]
                 + measurePhi[3]*measurePhi[3];
  measureStaggeredPhiNorm = measureStaggeredPhi[0]*measureStaggeredPhi[0]
                          + measureStaggeredPhi[1]*measureStaggeredPhi[1]
                          + measureStaggeredPhi[2]*measureStaggeredPhi[2]
                          + measureStaggeredPhi[3]*measureStaggeredPhi[3];
		 
  measurePhiNorm = sqrt(measurePhiNorm);
  measureStaggeredPhiNorm = sqrt(measureStaggeredPhiNorm);
  
  double logDetNeubergerWithXiFull;
  double logDetNeubergerWithXiNoZeros;
  double logDetWilsonWithOutXiFull;
  double logDetWilsonWithOutXiNoZeros;
  double logDetNeubergerWithOutXiFull;
  double logDetNeubergerWithOutXiNoZeros;
  double logDetWilsonWithXiFull;
  double logDetWilsonWithXiNoZeros;
  
  calculateDeterminants(logDetNeubergerWithXiFull, logDetNeubergerWithXiNoZeros, logDetWilsonWithOutXiFull, logDetWilsonWithOutXiNoZeros,
                        logDetNeubergerWithOutXiFull, logDetNeubergerWithOutXiNoZeros, logDetWilsonWithXiFull, logDetWilsonWithXiNoZeros);

  char* fileName = new char[300];

  if (Parameter_FLAG_ImprovedPhiDynamics) {
    double ImprovedPhiDynamicsDetQuotient = sqr(sqr(sqr(norm*measurePhiNorm)));
    double ImprovedPhiDynamicsLogDetSummand = -Parameter_Nf*log(ImprovedPhiDynamicsDetQuotient);
    
    logDetNeubergerWithXiFull += ImprovedPhiDynamicsLogDetSummand;
    logDetNeubergerWithXiNoZeros = NaN;
    logDetWilsonWithOutXiFull += ImprovedPhiDynamicsLogDetSummand;
    logDetWilsonWithOutXiNoZeros = NaN;
    logDetNeubergerWithOutXiFull += ImprovedPhiDynamicsLogDetSummand;
    logDetNeubergerWithOutXiNoZeros = NaN;
    logDetWilsonWithXiFull += ImprovedPhiDynamicsLogDetSummand;
    logDetWilsonWithXiNoZeros = NaN;
    
    snprintf(fileName,300,"%s/data/results/metroImprovedPhiDynamics/%sN%dNf%dKap%1.3fLam%1.3fY%1.3fRho%1.3fR%1.3f.dat",DataBaseDirectory,Parameter_filenamePrefix,Parameter_N,Parameter_Nf,Parameter_Kappa,Parameter_Lambda,Parameter_Y,Parameter_RHO,Parameter_R);	 
  } else {
    snprintf(fileName,300,"%s/data/results/metro/%sN%dNf%dKap%1.3fLam%1.3fY%1.3fRho%1.3fR%1.3f.dat",DataBaseDirectory,Parameter_filenamePrefix,Parameter_N,Parameter_Nf,Parameter_Kappa,Parameter_Lambda,Parameter_Y,Parameter_RHO,Parameter_R);	 
  }
  
  FILE* file;
  file = fopen(fileName,"a");
  fprintf(file,"%f %f %lf %lf %lf %lf %lf %lf %lf %lf\n",
          measurePhiNorm,measureStaggeredPhiNorm,
	  logDetNeubergerWithXiFull, logDetNeubergerWithXiNoZeros, logDetWilsonWithOutXiFull, logDetWilsonWithOutXiNoZeros,
	  logDetNeubergerWithOutXiFull, logDetNeubergerWithOutXiNoZeros, logDetWilsonWithXiFull, logDetWilsonWithXiNoZeros);
  fclose(file);  
  delete [] fileName;

  if (LogLevel>1) {
    print(measurePhiNorm);
    print(measureStaggeredPhiNorm);
  }
  if (LogLevel>1) printf("\n");
}


int main(int argc,char **argv) {
  int SweepCount;
  double acceptRate;
  int NConf = 0;

  loadParameters(0,false);


  if (LogLevel>1) printf("Number of arguments = %d\n",argc);
  int I;
  for (I=0; I<argc; I++) {
    if (LogLevel>1) printf("Argument %d: %s\n",I+1,argv[I]);  
  }
  
  int JobNr = -1;
  if (argc>=2) {
    if (sscanf(argv[1],"%d",&JobNr)!=1)  JobNr = -1;
  }
  int MultiProcessNr = -1;
  if (argc>=3) {
    if (sscanf(argv[2],"%d",&MultiProcessNr)!=1)  MultiProcessNr= -1;
  }
  
  if (JobNr>0) {
    if (LogLevel>1) printf("Use Job-Nr: %d for rand seed generation.\n",JobNr);
  } else {
    if (LogLevel>1) printf("Start rand seed from system time only.\n");
    JobNr = 1;
  }
  if (MultiProcessNr>0) {
    if (LogLevel>1) printf("Use multi-process number: %d for simulation parameter determination.\n",JobNr);
  } else {
    if (LogLevel>1) printf("Single process run.\n");
    MultiProcessNr = 0;
  }
  if (LogLevel>0) printf("\n");
  
  iniTools(JobNr);
  loadParameters(MultiProcessNr,true);
  iniFields();
  
  if (Parameter_Y>0) {
    fermiOps = new FermionMatrixOperations(Parameter_N, Parameter_N, Parameter_N, Parameter_N, Parameter_RHO, Parameter_R, 0);

    if (Parameter_FLAG_CalcNeubergerDetWithXi || Parameter_FLAG_CalcNeubergerDetWithoutXi) {
      fermiOps->constructNeubergerWithXiFermionMatrix(calcMatrix, phiField);
      if (LogLevel>3) {
        printf("Calculating (Neuberger) eigenvalues...");
        calcMatrix.calcEigenvalues();
        printf("ready.\n");
        printEigenvalues("DiracNeuberger",calcMatrix,ComplexUnity*(-1/(2*Parameter_RHO)));
      }
    }

    if (Parameter_FLAG_CalcWilsonDetWithoutXi || Parameter_FLAG_CalcWilsonDetWithXi) {
      fermiOps->constructWilsonWithOutXiFermionMatrix(calcMatrix, phiField);
      if (LogLevel>3)  {
        printf("Calculating (Wilson) eigenvalues...");
        calcMatrix.calcEigenvalues();
        printf("ready.\n");
        printEigenvalues("DiracWilson",calcMatrix,ComplexUnity*(-1/(2*Parameter_RHO)));
      }
    }
    
    fermiOps->setYukawaCoupling(Parameter_Y);    
  }
  
  if (LogLevel>1) printf("Performing %d thermalizing Monte Carlo Sweeps...\n",Parameter_thermalSweeps);
  acceptRate = 0;
  for (SweepCount=0; SweepCount<Parameter_thermalSweeps; SweepCount++) {
    acceptRate += MonteCarloSweep() / Parameter_thermalSweeps;
  }
  if (LogLevel>1) printf("...ready. Update quote: %1.2f\n",acceptRate);
  
  for (NConf=0; NConf<10000; NConf++) {
    if (LogLevel>2) printf("Performing %d Monte Carlo Sweeps...\n",Parameter_nSweeps);
    acceptRate = 0;
    for (SweepCount=0; SweepCount<Parameter_nSweeps; SweepCount++) {
      acceptRate += MonteCarloSweep() / Parameter_nSweeps;
    }
    if (LogLevel>2) printf("...ready. Update quote: %1.2f\n",acceptRate);
    if (acceptRate<0.20) {
      printf("Accept rate is smaller than 20 percent on average!!!\n");
      exit(0);
    }
   
    measure(); 
  }


  desini();
}
