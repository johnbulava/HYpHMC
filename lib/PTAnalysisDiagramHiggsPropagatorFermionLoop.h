#ifndef PTAnalysisDiagramHiggsPropagatorFermionLoop_included
#define PTAnalysisDiagramHiggsPropagatorFermionLoop_included

#include <math.h>

#include "Global.h"
#include "Tools.h"
#include "Complex.h"
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include "PTAnalysisDiagram.h"
#include "NeubergerMatrix.h"


class PTAnalysisDiagramHiggsPropagatorFermionLoop : PTAnalysisDiagram {
private:
	double rho;
	ComplexMatrix gamma5;
	ComplexMatrix rProjector;
	ComplexMatrix lProjector;
	
	NeubergerMatrix* diracOp;
	ComplexVector v[4];		
	ComplexMatrix getZetaTimesPropagator(vector4D k, double mF);

public:	
	PTAnalysisDiagramHiggsPropagatorFermionLoop(int l0, int l1, int l2, int l3, int typeOfOp); 
	~PTAnalysisDiagramHiggsPropagatorFermionLoop();
	bool isLastValidResultStillValid(vector4D p, double mF0, double mH0, double mG0, double y0, double lambda0, double massSplit);
	ComplexMatrix calcDiagramContributionToPropagator(vector4D p, double mF0, double mH0, double mG0, double y0, double lambda0, double massSplit);
};


#endif
