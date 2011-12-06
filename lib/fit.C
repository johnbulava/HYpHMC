#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include <pthread.h>

#include "Tools.C"


char* fileName = NULL;
int column1 = 0;
int column2 = 0;
int column3 = 0;
int NrOfTerms = 0;
int FitIterations = 0;
char** terms = NULL;
int dataCount = 0;
double** data;
double* fitConst = NULL;
double* fitConstError = NULL;
double fitRedChiSqr = 0;
double* dataX = NULL;
double* dataY = NULL;
double* dataYErr = NULL;






void invalidParameters() {
  printf("Parameters must be: filename \n");
  exit(0);
}


void loadCommandLineParameters(int argc,char **argv) {
  bool error = false;
  if (argc<6) invalidParameters();
  fileName = new char[2000];
  snprintf(fileName, 2000, "%s", argv[1]);
  if (sscanf(argv[2],"%d",&column1)!=1) error = true;
  if (sscanf(argv[3],"%d",&column2)!=1) error = true;
  if (sscanf(argv[4],"%d",&column3)!=1) error = true;
  if (sscanf(argv[5],"%d",&FitIterations)!=1) error = true;
  NrOfTerms = argc - 6;
  
  terms = new char*[NrOfTerms];
  for (int I=0; I<NrOfTerms; I++) {
    terms[I] = new char[2000];
    if (sscanf(argv[6+I],"%s",terms[I])!=1) error = true;    
  }
  
  if (error) invalidParameters();
  
  printf("Read Parameters: \n");
  printf("  -> Column1       : %d\n", column1);
  printf("  -> Column2       : %d\n", column2);
  printf("  -> Column3       : %d\n", column3);
  printf("  -> FitIterations : %d\n", FitIterations);
  printf("  -> Nr of Terms   : %d\n", NrOfTerms);
  for (int I=0;I<NrOfTerms; I++) {
    printf("  -> term %d        : %s\n", I, terms[I]);
  }
}


void loadData() {
  printf("Loading file %s\n", fileName);
  int colMax = column1;
  if (column2>colMax) colMax = column2;
  if (column3>colMax) colMax = column3;  
  
  FILE* file = fopen(fileName, "r");
  double dummy = 0;
  dataCount = 0;
  bool error = false;
  char* restLine = new char[3000];
  while (!error) {
    for (int I=0; I<colMax; I++) {
      if (fscanf(file, "%lf", &dummy) != 1) {
        error = true;
        break;
      }
    }
    dataCount++;
    fgets(restLine, 2500, file);
  }
  dataCount--;  
  fclose(file);

  data = new double*[dataCount];
  file = fopen(fileName, "r");
  dataCount = 0;
  error = false;
  while (!error) {
    data[dataCount] = new double[colMax];
    for (int I=0; I<colMax; I++) {
      if (fscanf(file, "%lf", &(data[dataCount][I])) != 1) {
        error = true;
        break;
      }
    }
    dataCount++;
    fgets(restLine, 2500, file);
  }
  dataCount--;
  fclose(file);
  
  delete[] restLine;
  
  printf("Read %d data lines\n", dataCount);
  dataX = new double[dataCount];
  dataY = new double[dataCount];
  dataYErr = new double[dataCount];
  for (int I=0; I<dataCount; I++) {
    dataX[I] = data[I][column1-1];
    dataY[I] = data[I][column2-1];
    dataYErr[I] = data[I][column3-1];
    
    printf(" **> %1.3e %1.3e %1.3e\n",dataX[I],dataY[I],dataYErr[I]);
  }
}


void performFit() {
  fitConst = new double[NrOfTerms];
  fitConstError = new double[NrOfTerms];
  fitRedChiSqr = 0;
  char* functionBody = new char[2000];
  char* dummy = new char[2000];
  

  printf("Finding initial fit parameters...\n");
  snprintf(functionBody,2000," ");
  for (int I=0; I<NrOfTerms; I++) {
    snprintf(dummy, 2000, "%s", functionBody);
    snprintf(functionBody,2000,"%s", terms[I]);
     
    fitConst[I] = 1;
    printf("...with fit function %s\n", functionBody);
    performGnuplotFit(functionBody, dataX, dataY, dataYErr, dataCount, I+1, fitConst, fitConstError, fitRedChiSqr);
  }
  
  printf("\nInitial fit parameters are:\n");
  for (int I=0; I<NrOfTerms; I++) {
    printf("  A%d = %1.3e +- %1.3e\n", I+1, fitConst[I], fitConstError[I]);
  }

  delay(2);
  
  performGnuplotFitWithErrorEstimateFromResampling(functionBody,dataX, dataY, dataYErr, dataCount, NrOfTerms, fitConst, fitConstError, fitRedChiSqr, FitIterations); 

  printf("\nFinal fit parameters are:\n");
  for (int I=0; I<NrOfTerms; I++) {
    printf("  A%d = %1.3e +- %1.3e (%1.1f %%)\n", I+1, fitConst[I], fitConstError[I], 100*fitConstError[I]/ abs(fitConst[I]));
  }

  
  delete[] functionBody;
  delete[] dummy;
}


int main(int argc, char **argv) {
  loadCommandLineParameters(argc, argv);
  LogLevel = 0;
  loadData();  
  iniTools(5517);

  
  performFit();
  
  
  
  delete[] fileName;
  for (int I=0; I<NrOfTerms; I++) {
    delete[] terms[I];
  }
  delete[] terms;
  delete[] fitConst;
  delete[] fitConstError;
  delete[] dataX;
  delete[] dataY;
  delete[] dataYErr;
}
