#include "PTAnalysisDiagram.h"

PTAnalysisDiagram::PTAnalysisDiagram(int l0, int l1, int l2, int l3, int typeOfOp) { 
  L0 = l0;
  L1 = l1;
  L2 = l2;
  L3 = l3;
  TypeOfOperator = typeOfOp;
  lastValid_Result = new ComplexMatrix(1);
  
  lastValid_p[0] = 0;
  lastValid_p[1] = 0;
  lastValid_p[2] = 0;
  lastValid_p[3] = 0;
  lastValid_mF0 = -1;
  lastValid_mH0 = -1;
  lastValid_mG0 = -1;
  lastValid_y0 = -1;
  lastValid_lambda0 = -1;
  lastValid_massSplit = -1;  
}


PTAnalysisDiagram::~PTAnalysisDiagram() { 
  delete lastValid_Result;
}


ComplexMatrix PTAnalysisDiagram::getDiagramContributionToPropagator(vector4D p, double mF0, double mH0, double mG0, double y0, double lambda0, double massSplit) {
  if (isLastValidResultStillValid(p, mF0, mH0, mG0, y0, lambda0, massSplit)) {
    ComplexMatrix res(*lastValid_Result);
    return res;
  } 
  
  ComplexMatrix res = calcDiagramContributionToPropagator(p, mF0, mH0, mG0, y0, lambda0, massSplit);
  *lastValid_Result = res;
  lastValid_p[0] = p[0];
  lastValid_p[1] = p[1];
  lastValid_p[2] = p[2];
  lastValid_p[3] = p[3];
  lastValid_mF0 = mF0;
  lastValid_mH0 = mH0;
  lastValid_mG0 = mG0;
  lastValid_y0 = y0;
  lastValid_lambda0 = lambda0;
  lastValid_massSplit = massSplit;  
  return res;
}
