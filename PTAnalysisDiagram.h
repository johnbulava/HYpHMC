#ifndef PTAnalysisDiagram_included
#define PTAnalysisDiagram_included

#include <math.h>

#include "Global.h"
#include "Tools.h"
#include "Complex.h"
#include "ComplexVector.h"
#include "ComplexMatrix.h"


class PTAnalysisDiagram {
private:
  ComplexMatrix* lastValid_Result;


protected:
  int L0;
  int L1;
  int L2;
  int L3;
  int TypeOfOperator;        // 0: Neuberger-Operator, 1: Wilson-Operator

  vector4D lastValid_p;
  double lastValid_mF0;
  double lastValid_mH0;
  double lastValid_mG0;
  double lastValid_y0;
  double lastValid_lambda0;
  double lastValid_massSplit;
  
  virtual bool isLastValidResultStillValid(vector4D p, double mF0, double mH0, double mG0, double y0, double lambda0, double massSplit) = 0;
  virtual ComplexMatrix calcDiagramContributionToPropagator(vector4D p, double mF0, double mH0, double mG0, double y0, double lambda0, double massSplit) = 0;
  
public:
  PTAnalysisDiagram(int l0, int l1, int l2, int l3, int typeOfOp); 
  virtual ~PTAnalysisDiagram();
  
  ComplexMatrix getDiagramContributionToPropagator(vector4D p, double mF0, double mH0, double mG0, double y0, double lambda0, double massSplit);

};


#endif
