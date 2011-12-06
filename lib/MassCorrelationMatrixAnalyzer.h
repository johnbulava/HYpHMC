#ifndef MassCorrelationMatrixAnalyzer_included
#define MassCorrelationMatrixAnalyzer_included

#include <math.h>
#include <fstream>
#include "Complex.h"
#include "Quat.h"
#include "Global.h"
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include "Tools.h"
#include "ControlLogger.h"


class MassCorrelationMatrixAnalyzer {
protected:
  int dimension;
  int timeExtent;
char* uniqueAnalyzerFileNameExtension;
  bool operatorCorrelationMode;
  bool generalEWproblem;
  bool iterativeEigenVectorAlignment;
  bool subtractVacuumExpectation;  
  bool fitVacuumExpectation;    
  int dataMax;
  int dataCount;
  int* dataID;
  double* dataWeights;
  Complex*** OperatorData;
  Complex**** OperatorCorrelationData;
  ComplexMatrix** MassCorrelationMatrices;
  ComplexMatrix* GeneralEWNormalizationMatrix;

  Complex** MassCorrelationEigenvalues;
    Complex** MassCorrelationEigenvaluesError;
  double** effectiveMasses;
    double** effectiveMassesError;
  double* asymptoticEffectiveMasses;
double* asymptoticEffectiveMassesError;  
  double* fittedVacuumExpectationValue;
    double* fittedVacuumExpectationValueError;
  int fittedMassCount;
  double** fittedMasses;
    double** fittedMassesError;
  double** fittedCoefficients;
    double** fittedCoefficientsError;
  double* fitChiSquare;
    double* fitChiSquareError;
  int fittingRangeReduction;


  void sortData();
  void calcMassCorrelationMatrices(int ignoreStart, int ignoreEnd, bool boot);
  void calcOptimalEVassignment(ComplexVector** v1, ComplexVector** v2, int* assignment);
  bool calcMassCorrelationEigenvalues(int ignoreStart, int ignoreEnd, bool boot, Complex** eigenvalues, ComplexVector** sortOrder);
  void calcEigenvaluesAndErrors(bool boot, int iterations);
  void calcMassesAndErrors(bool boot, int iterations, int massCount);
  double calcErrorRescaleFactor(bool boot, int iterations);
  int calcDecentJackKnifeIterations();
  
  
  

public:
  MassCorrelationMatrixAnalyzer(int dim, int timeEx, bool opCorMode, bool genEWprob, bool iterEValign, bool SubVac, bool FitVac, int fitRangeRed, int datMx, char* uniqueID);
  ~MassCorrelationMatrixAnalyzer();

  void clearData();
  void addOperatorData(int ID, double weight, Complex** opData);
  void addOperatorCorrelationData(int ID, double weight, Complex*** opCorData);
  void addOperatorCorrelationData(int ID, double weight, Complex*** opCorData, Complex** opDagData, Complex** opData);  
  int getDataCount();

  void calcEigenvaluesWithBootStrapAnalysis(int iterations);
  void calcEigenvaluesWithJackKnifeAnalysis(int iterations);
  void calcEigenvaluesWithJackKnifeAnalysis();
  void calcEigenvaluesAndMassesWithBootStrapAnalysis(int iterations, int massCount);
  void calcEigenvaluesAndMassesWithJackKnifeAnalysis(int iterations, int massCount);
  void calcEigenvaluesAndMassesWithJackKnifeAnalysis(int massCount);
    
void plotEigenvalues(bool logScale);
void plotEffectiveMasses();
void plotEigenvalues(ControlLogger* logger, bool logScale);
void plotEffectiveMasses(ControlLogger* logger);
  
  int getFittedMassCount();
  int getDimension();
  int getTimeExtent();
  double getFittedMass(int operatorIndex, int massIndex);
  double getFittedMassError(int operatorIndex, int massIndex);
  double getFittedChiSquare(int operatorIndex);
  double getFittedChiSquareError(int operatorIndex);
  double getFittedCoefficient(int operatorIndex, int massIndex);
  double getFittedCoefficientError(int operatorIndex, int massIndex);
  double getFittedVacuumExpectationValue(int operatorIndex);
  double getFittedVacuumExpectationValueError(int operatorIndex);
  double getFittedEffectiveMasses(int timeIndex, int operatorIndex);
  double getFittedEffectiveMassesError(int timeIndex, int operatorIndex);
  double getFittedAsymptoticEffectiveMasses(int operatorIndex);
  double getFittedAsymptoticEffectiveMassesError(int operatorIndex);  
  double getMassCorrelationEigenvalue(int timeIndex, int operatorIndex);
  double getMassCorrelationEigenvalueError(int timeIndex, int operatorIndex);
  
  
void plotEigenvalues();
};

#endif
