#ifndef AnalyzerTop_included
#define AnalyzerTop_included

#include <stdlib.h>

#include "Global.h"
#include "Complex.h"
#include "Quat.h"
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include FFTWIncludeFile
#include "AutoCorrelation.h"
#include "FermionMatrixOperations.h"
#include "ControlLogger.h"
#include "Tools.C"

#define AnalyzerTopDataMAX 1000000

class AnalyzerTop {
private:
  void ini(int l0, int l1, int l2, int l3, double y, ControlLogger* log, bool FLAG_GFit, char* fileNameExt, double tol);
  void desini();
  void sampleFermioncScanVector(Complex* v, int t, int FermionIndex);
  void calcPsiPsiBarMatrixForPhiField(vector4D* phiField, Complex**** PsiPsiBarMatrix);
  int confInPoolIndex(int confNr);
  
  ControlLogger*  MassControlLog;
  FermionMatrixOperations* fermiOps;
  Complex* LeftVector;
  Complex* RightVector;
  Complex* SolutionVector;
  int L0,L1,L2,L3;
  double Parameter_Y;
  double Parameter_R;
  double Parameter_RHO;
  int timeDirection;
  int LargestL;
  double* sinPSqr;
  AutoCorrelation** TopTimeSliceCorrelator;
  Complex** TopTimeSliceCorrelatorData;
  int* TopTimeSliceCorrelatorConfNr;
  int* confPoolIndOrder;
  int TopTimeSliceCorrelatorDataCount;
  double* weightData;
  int totalN;
  bool FLAG_GnuplotFit;
  AutoCorrelation* TopMass;
  double* TopMassData;
  char* confFileNameExtension;
  double TOL;
  
  MassCorrelationMatrixAnalyzer* TopMassAnalyzer;
  
public:
  AnalyzerTop(int l0, int l1, int l2, int l3, double y, ControlLogger* log, bool FLAG_GFit, char* fileNameExt, double tol);
  ~AnalyzerTop();
  
  double LatticeResult_PhysicalTopMass;
  double LatticeResult_PhysicalTopMassError;
  
  int AnalyzerTop::getTotalN();
  void analyzeHiggsField(vector4D* phiField, double weight, int confNr);

  void calcTopTimeSliceCorrelator();
  void plotTopTimeSliceCorrelator();
  void loadTopTimeSliceCorrelator();
  void saveTopTimeSliceCorrelator();
void plotTopMasses();

};

#include "AnalyzerTop.C"

#endif
