#include "PTAnalysisDiagramHiggsPropagatorFermionLoop.h"

PTAnalysisDiagramHiggsPropagatorFermionLoop::PTAnalysisDiagramHiggsPropagatorFermionLoop(int l0, int l1, int l2, int l3, int typeOfOp) : PTAnalysisDiagram(l0, l1, l2, l3, typeOfOp) {
	rho = 1.0;
	diracOp = new NeubergerMatrix(rho, 0.5, 1, 1, 1, 1, 2);	
	v[0].resize(4);
	v[1].resize(4);
	v[2].resize(4);
	v[3].resize(4);
		 
	gamma5.resize(4);
	gamma5.matrixElements[0].x = 1.0;
	gamma5.matrixElements[5].x = 1.0;
	gamma5.matrixElements[10].x = -1.0;
	gamma5.matrixElements[15].x = -1.0;
	
	rProjector.resize(4);
	rProjector.matrixElements[0].x = 1.0;
	rProjector.matrixElements[5].x = 1.0;	
	
	lProjector.resize(4);
	lProjector.matrixElements[10].x = 1.0;
	lProjector.matrixElements[15].x = 1.0;	
}


PTAnalysisDiagramHiggsPropagatorFermionLoop::~PTAnalysisDiagramHiggsPropagatorFermionLoop() {
	delete diracOp; 
}

ComplexMatrix PTAnalysisDiagramHiggsPropagatorFermionLoop::getZetaTimesPropagator(vector4D k, double mF) {
	diracOp->analyticalEigenvectors(k, v);
	Complex ew_massless = diracOp->analyticalEigenvalue(k);
	ComplexMatrix m1(v[0]);  
	ComplexMatrix m2(v[1]);  
	ComplexMatrix m3(v[2]);  
	ComplexMatrix m4(v[3]);	
	
	ComplexMatrix propagator(4);
	Complex propEW(0.0, 0.0);
	Complex propEWInv(0.0, 0.0);
	Complex propEWInvAdj(0.0, 0.0);
	
	ComplexMatrix zeta(4);
	Complex zetaEW(0.0, 0.0);
	Complex zetaEWAdj(0.0, 0.0);
	
	propEW.x = ew_massless.x + mF - mF*ew_massless.x/(2.0*rho);  
	propEW.y = ew_massless.y - mF*ew_massless.y/(2.0*rho);
	 
	propEWInv.x = propEW.x/(propEW.x*propEW.x + propEW.y*propEW.y);
	propEWInv.y = propEW.y/(propEW.x*propEW.x + propEW.y*propEW.y);		
	propEWInvAdj.x = propEWInv.x;
	propEWInvAdj.y = -1.0*propEWInv.y;
	
	for (int j = 0; j < 16; j++) {
		propagator.matrixElements[j] = propEWInv * m1.matrixElements[j] + propEWInv * m2.matrixElements[j] + propEWInvAdj * m3.matrixElements[j] + propEWInvAdj * m4.matrixElements[j];
	}
	
	zetaEW.x = 1.0 - ew_massless.x/(2.0*rho);
	zetaEW.y = -1.0 * ew_massless.y/(2.0*rho);
	zetaEWAdj.x = zetaEW.x;
	zetaEWAdj.y = -1.0 * zetaEW.y;
	
	for (int j = 0; j < 16; j++) {
		zeta.matrixElements[j] = zetaEW * m1.matrixElements[j] + zetaEW * m2.matrixElements[j] + zetaEWAdj * m3.matrixElements[j] + zetaEWAdj * m4.matrixElements[j];
	}
	
	return propagator * zeta; 
}


bool PTAnalysisDiagramHiggsPropagatorFermionLoop::isLastValidResultStillValid(vector4D p, double mF0, double mH0, double mG0, double y0, double lambda0, double massSplit) {
  return false;
}


ComplexMatrix PTAnalysisDiagramHiggsPropagatorFermionLoop::calcDiagramContributionToPropagator(vector4D p, double mF0, double mH0, double mG0, double y0, double lambda0, double massSplit) {
  ComplexMatrix res(1);
  
  vector4D k;  
  vector4D pPlusk;  
  
  Complex tmp (0,0);
  Complex trace(1,0);
  
  for (int i0=0; i0<L0; i0++) {
	  k[0] = 2*pi*i0 / L0;
	  pPlusk[0] = p[0] + k[0];
	  for (int i1=0; i1<L1; i1++) {
		  k[1] = 2*pi*i1 / L1;
		  pPlusk[1] = p[1] + k[1];
		  for (int i2=0; i2<L2; i2++) {
			  k[2] = 2*pi*i2 / L2;
			  pPlusk[2] = p[2] + k[2];
			  for (int i3=0; i3<L3; i3++) {
				  k[3] = 2*pi*i3 / L3;
				  pPlusk[3] = p[3] + k[3];
	  
				  ComplexMatrix mat1 = getZetaTimesPropagator(pPlusk, mF0);
				  ComplexMatrix mat2 = getZetaTimesPropagator(k, mF0);
				  mat1 = mat1 * mat2;
				  trace = mat1.matrixElements[0] + mat1.matrixElements[5] + mat1.matrixElements[10] + mat1.matrixElements[15];
				  tmp = tmp + trace;
			  }
		  }
	  }
  }
  tmp = (2.0/(L0*L1*L2*L3)) * tmp;  // factor of 2 due to top and bottom contributions
  res.matrixElements[0] = -1.0*tmp;
  
  return res;
}
