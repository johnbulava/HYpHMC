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



//Variables
int Parameter_L0 = 0;
int Parameter_L1 = 0;
int Parameter_L2 = 0;
int Parameter_L3 = 0;
double Parameter_RHO = 0;
double Parameter_R = 0;
double Parameter_Y = 0;
double Parameter_MassSplit = 0;
char* Parameter_filenamePrefix;
int Parameter_FLAG_CalcNeubergerDetWithXi = 0;
int Parameter_FLAG_CalcWilsonDetWithoutXi = 0;
int Parameter_FLAG_CalcNeubergerDetWithoutXi = 0;
int Parameter_FLAG_CalcWilsonDetWithXi = 0;
vector4D* phiField = NULL;
FermionMatrixOperations* fermiOps = NULL;
ComplexMatrix calcMatrix(1);
int RunCount = 0;



void iniFields() {
  if (LogLevel>0) printf("Initializing fields...");
  phiField = new vector4D[Parameter_L0*Parameter_L1*Parameter_L2*Parameter_L3];  
  int I;
  for (I=0; I<Parameter_L0*Parameter_L1*Parameter_L2*Parameter_L3; I++) {
    phiField[I][0] = 0;
    phiField[I][1] = 0;
    phiField[I][2] = 0;
    phiField[I][3] = 0;
  }
  if (LogLevel>0) printf("sucessfully.\n");    
}


void randomPhi() {
  if (LogLevel>0) printf("Randomizing fields...\n");
  fermiOps->fillGaussRandomVector((Complex*) phiField, 2*Parameter_L0*Parameter_L1*Parameter_L2*Parameter_L3);
  for (int I=0; I<Parameter_L0*Parameter_L1*Parameter_L2*Parameter_L3; I++) {
    printf("%f %f %f %f\n", phiField[I][0], phiField[I][1], phiField[I][2], phiField[I][3]);
  }
  if (LogLevel>0) printf("sucessfully.\n");    
}



void calculateDeterminants() {
  int I;
  Complex normedDet;
  double logDet;
  FILE* file1;
  FILE* file2;
  char* fileName1 = new char[1000];
  char* fileName2 = new char[1000];
  

  if (Parameter_FLAG_CalcNeubergerDetWithXi) {
    if (LogLevel>2) printf("Calculating Neuberger Operator with Xi-fields...\n");
    fermiOps->constructNeubergerWithXiFermionMatrix(calcMatrix, phiField, Parameter_MassSplit);
  
    if (LogLevel>2) printf("Calculating Neuberger Determinant with Xi-fields...\n");
    calcMatrix.calcEigenvalues();
   
    snprintf(fileName1,1000,"%s/data/results/spectrum/%sNeubergerXiEigenValuesL%dx%dx%dx%dY%1.3fSplit%1.3fRho%1.3fR%1.3fRun%d.dat", DataBaseDirectory,Parameter_filenamePrefix,Parameter_L0,Parameter_L1,Parameter_L2,Parameter_L3,Parameter_Y,Parameter_MassSplit,Parameter_RHO,Parameter_R,RunCount);
    snprintf(fileName2,1000,"%s/data/results/spectrum/%sNeubergerXiDeterminantL%dx%dx%dx%dY%1.3fSplit%1.3fRho%1.3fR%1.3f.dat", DataBaseDirectory,Parameter_filenamePrefix,Parameter_L0,Parameter_L1,Parameter_L2,Parameter_L3,Parameter_Y,Parameter_MassSplit,Parameter_RHO,Parameter_R);
    file1 = fopen(fileName1,"w");
    file2 = fopen(fileName2,"a");
    
    normedDet = Complex(1.0, 0);
    logDet = 0;
    for (I=0; I<calcMatrix.matrixSize; I++) {
      double norm = sqrt(sqr(calcMatrix.eigenvalues[I].x) + sqr(calcMatrix.eigenvalues[I].y));
      if (norm > 1E-6) {      
        normedDet = (normedDet * calcMatrix.eigenvalues[I]) / norm;
        logDet += log(norm);
      } else {
        logDet = NaN;
      }
      fprintf(file1,"%1.15f %1.15f\n",calcMatrix.eigenvalues[I].x, calcMatrix.eigenvalues[I].y);
    }
    
    fprintf(file2,"%1.15f %1.15f %1.15f\n",normedDet.x, normedDet.y, logDet);
    
    fclose(file1);
    fclose(file2);
    if (LogLevel>2) printf("...successfully!\n");
  }
  
  if (Parameter_FLAG_CalcWilsonDetWithoutXi) {
    fermiOps->constructWilsonWithOutXiFermionMatrix(calcMatrix, phiField);
  
  }

  if (Parameter_FLAG_CalcNeubergerDetWithoutXi) {
    fermiOps->constructNeubergerWithOutXiFermionMatrix(calcMatrix, phiField);
  
  }
  
  if (Parameter_FLAG_CalcWilsonDetWithXi) {
    fermiOps->constructWilsonWithXiFermionMatrix(calcMatrix, phiField);
  
  } 
  
  delete[] fileName1;
  delete[] fileName2;
}


void spacetimeReversePhi(int dir) {
  int ind[4];
  int L[4];
  L[0] = Parameter_L0;
  L[1] = Parameter_L1;
  L[2] = Parameter_L2;
  L[3] = Parameter_L3;
  
  for (ind[0]=0; ind[0]<Parameter_L0; ind[0]++) {
    for (ind[1]=0; ind[1]<Parameter_L1; ind[1]++) {
      for (ind[2]=0; ind[2]<Parameter_L2; ind[2]++) {
        for (ind[3]=0; ind[3]<Parameter_L3; ind[3]++) {
	  if (ind[dir]<L[dir]/2) {
            int index = ind[3] + ind[2]*Parameter_L3 + ind[1]*Parameter_L2*Parameter_L3 + ind[0]*Parameter_L1*Parameter_L2*Parameter_L3;
            ind[dir] = L[dir]-1 - ind[dir];
            int revIndex = ind[3] + ind[2]*Parameter_L3 + ind[1]*Parameter_L2*Parameter_L3 + ind[0]*Parameter_L1*Parameter_L2*Parameter_L3;
            ind[dir] = L[dir]-1 - ind[dir];	
	  
	    for (int i=0; i<4; i++) {
	      double dummy = phiField[index][i];
	       phiField[index][i] = phiField[revIndex][i];
	       phiField[revIndex][i] = dummy;
	    }
	  }
	}
      }
    }
  }
  if (LogLevel>0) printf("Reversed field in direction %d\n",dir);
  for (int I=0; I<Parameter_L0*Parameter_L1*Parameter_L2*Parameter_L3; I++) {
    printf("%f %f %f %f\n", phiField[I][0], phiField[I][1], phiField[I][2], phiField[I][3]);
  }
}


void addConstantToPhi(double c) {
  for (int I=0; I<Parameter_L0*Parameter_L1*Parameter_L2*Parameter_L3; I++) {
    phiField[I][0] += c;
    phiField[I][1] += c;
    phiField[I][2] += c;
    phiField[I][3] += c;
  }
}


int main(int argc,char **argv) {
  LogLevel = 3;

  Parameter_L0 = 4;
  Parameter_L1 = 4;
  Parameter_L2 = 4;
  Parameter_L3 = 8;
  Parameter_RHO = 1.0;
  Parameter_R = 0.5;
  Parameter_Y = 2.000;
  Parameter_MassSplit = 0.024;
  Parameter_filenamePrefix = "spectrum";
  Parameter_FLAG_CalcNeubergerDetWithXi = 1;
  Parameter_FLAG_CalcWilsonDetWithoutXi = 0;
  Parameter_FLAG_CalcNeubergerDetWithoutXi = 0;
  Parameter_FLAG_CalcWilsonDetWithXi = 0;
  
  iniTools(-1);
  iniFields();
  
  fermiOps = new FermionMatrixOperations(Parameter_L0, Parameter_L1,Parameter_L2,Parameter_L3,Parameter_RHO, Parameter_R, 0);
  fermiOps->setYukawaCoupling(Parameter_Y);    
  
  RunCount = 0;
  while (true) {
    randomPhi();
//addConstantToPhi(1);
    calculateDeterminants();
    
/*    spacetimeReversePhi(3);
    
    calculateDeterminants();*/
    
    RunCount++;
    if (RunCount>100) break;
  }
}
