#ifndef AnalyzerPhiFieldConfiguration_included
#define AnalyzerPhiFieldConfiguration_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "Global.h"
#include "Tools.h"
#include "FermionMatrixOperations.h"



class AnalyzerPhiFieldConfiguration {
private:
  FermionMatrixOperations* fermiOps;
  double* phiField;
  double** phiFieldCopies;
  int phiFieldCopyCount;
  double weight;
  int errorState;
  bool weightAvail;
  double magnetizationM;
  double magnetizationS;
  double phiFieldAvgVectorLength;
  double phiFieldAvgVectorLengthVariation;  
  vector4D avgPhiFieldVector;
  
  void loadPhiconfiguration(char* fileName);
  void getHiggsFieldDirection(double* pField, vector4D dir);
  void rotateHiggsField(double* pField, int ind1, int ind2, double w);
  void measureMagnetizations();

public:
  AnalyzerPhiFieldConfiguration(char* fileName, FermionMatrixOperations* fOps); 
  ~AnalyzerPhiFieldConfiguration();
  
  double* getPhiFieldCopy();  //Only delete by calling clearPhiFieldCopies()
  void clearPhiFieldCopies();

  int getErrorState();
  bool isWeightAvail();
  double getWeight();
  double getMagnetizationM();
  double getMagnetizationS();
  double getPhiFieldAvgVectorLength();
  double getPhiFieldAvgVectorLengthVariation();
  double getAvgPhiFieldVectorComponent(int comp);
  double getAvgPhiFieldVectorLength();
  
  
  void multiplyHiggsFieldWithConst(double* pField, double fac);
  void alignHiggsFieldDirection(double* pField);
  void randomGaugeRotation(double* pField);  
  void randomizeHiggsField(double* pField);
  void randomizeGaussHiggsField(double* pField);

  
  Complex* performFourierTransform(double* pField, bool forward, int howmanyComponents);  
};


#endif
