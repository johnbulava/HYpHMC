#ifndef AnalyzerObservableDetSign_included
#define AnalyzerObservableDetSign_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "AnalyzerIOControl.h"
#include "FermionMatrixOperations.h"
#include "StateDescriptorReader.h"
#include "AnalyzerObservable.h"
#include "AnalyzerPhiFieldConfiguration.h"


#define LARGEST_REAL_VALUE 0		
#define SMALLEST_REAL_VALUE 1
#define LARGEST_ABSOLUT_VALUE 4
#define SMALLEST_ABSOLUT_VALUE 5

class AnalyzerObservableDetSign : public AnalyzerObservable {

private:  
	va_list args;
	char* str;
	int max_iterations;
	int curr_eigenvalue_count;
	int eigenvalue_increment;
	int logLevel;
	double eigenvalue_prec;
	bool usePreconditionerP;
	
	int getDeterminant_old (FermionMatrixOperations& fermiOps, double* phiField, int initialEigenValueCount, double tol, Complex& det);
	int getSignOfDeterminant (FermionMatrixOperations& fermiOps, double* phiField, int initialEigenValueCount, double tol, Complex& det);
	void println(int level, char * format, ...);
	
public:    
  double* phiField_2save; // only for testing: this and the following line shoule be deleted for later use
  bool neg_det;
  double* extPhiField;
  
  AnalyzerObservableDetSign(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
  ~AnalyzerObservableDetSign();
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();      
};


#include "AnalyzerObservableDetSign.C"

#endif
