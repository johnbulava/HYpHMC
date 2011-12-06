// Date: December 2007
// Institut fuer theoretische Physik
// Humboldt Universitaet zu Berlin
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <iostream>

#include "PTAnalysis.h"
using namespace std;

inline double getPSqr(vector4D p) {
	double res = 0.0;
	double tmp = 0;
	
	tmp = 2* sin(0.5*p[0]);
	tmp *= tmp;	
	res += tmp; 
	
	tmp = 2* sin(0.5*p[1]);
	tmp *= tmp;	
	res += tmp;
	
	tmp = 2* sin(0.5*p[2]);
	tmp *= tmp;	
	res += tmp;
	
	tmp = 2* sin(0.5*p[3]);
	tmp *= tmp;	
	res += tmp;
		
	return res;
}

void createSelfEnergyData() {
	int l0,l1,l2,l3;
	l0=4;
	l1=4;
	l2=4;
	l3=8;	
	
	int momCuttOff = 0;
	
	double kappa = 0.12390;
	double yN = 0.1;
	double vev = 1.217;
	double field_renormalization = 0.99924;
	double lambda0 = 0.0;
	double massSplit = 1.0;
	double mH0 = 0.0;
	double mG0 = 0.0;
	
	double y0 = yN/sqrt(2 * kappa);
	double mF0 = vev * yN;		
	
	PTAnalysisDiagramHiggsPropagatorFermionLoop higgsFermiLoop(l0,l1,l2,l3,0); //l0, l1, l2, l3, typeOfOp.: 0 = neuberger, 1 = Wilson	
			
	vector4D p;
	p[0] = 0.0;
	p[1] = 0.0;
	p[2] = 0.0;
	p[3] = 0.0;
	
	ComplexMatrix mat(1);
	double selfEnergy[2][l3*l2*l1*l0];
	int counter = 0;	
	
	char* fileName = new char[1024];
	sprintf(fileName, "./PTAnalysisData/SelfEnergy_%dx%d_vev%1.5f_kappa%1.5f_yN%1.5f_mF%1.5f.dat", l0,l3,vev, kappa, yN, mF0);
	FILE* pFile = fopen(fileName, "w");	
		
	for (int i_3 = 0; i_3 < l3 - momCuttOff; i_3++) {
		p[3] = 2*pi*i_3 / l3;
		for (int i_2 = 0; i_2 < l2 - momCuttOff; i_2++) {
			p[2] = 2*pi*i_2 / l2;
			for (int i_1 = 0; i_1 < l1 - momCuttOff; i_1++) {
				p[1] = 2*pi*i_1 / l0;
				for(int i_0 = i_1; i_0 < l0 - momCuttOff; i_0++) {
					p[0] = 2*pi*i_0 / l0;
					mat = higgsFermiLoop.calcDiagramContributionToPropagator(p, mF0, mH0, mG0, y0, lambda0, massSplit);
					selfEnergy[0][counter] = getPSqr(p);
					selfEnergy[1][counter] = mat.matrixElements[0].x;
					fprintf(pFile, "%1.15f \t%1.15f \t%1.15f\n", selfEnergy[0][counter], selfEnergy[1][counter], y0*y0*selfEnergy[1][counter]/field_renormalization); 
					counter++;			
				}	
				fflush(pFile);			
			}
		}
	}	
	fclose(pFile);
	
	// rewrite the file and sort the data points:
	pFile = fopen(fileName, "w");		
	delete[] fileName;
	double precision = 1E-10;
	double x_value = 0.0;
	double y_value = 0.0;
	int occurances = 1;
			
	for(int i = 0; i < counter; i++) {
		if (selfEnergy[0][i] != -1) {
			x_value = selfEnergy[0][i];
			y_value = selfEnergy[1][i];
			occurances = 1;
			
			for (int j = i+1; j < counter; j++) {
				if (fabs(x_value - selfEnergy[0][j]) <= precision) {
					y_value += selfEnergy[1][j];
					selfEnergy[0][j] = -1.0;
					occurances++; 
				}
			}
			y_value /= occurances;
			fprintf(pFile, "%1.15f \t%1.15f \t%1.15f\n", x_value, y_value, y0*y0*y_value/field_renormalization);
		}
	}
	
	fclose(pFile);
}

void printMatrix(ComplexMatrix& m) {
	int size = 4;
	printf("Contents of ComplexMatrix:\n");
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			printf("%+1.6f %+1.6f i\t", m.matrix[i][j].x, m.matrix[i][j].y);
		}
		printf("\n");
	}
}

int main(void) {
	cout << "Pertubation Theory Analysis." << endl;			
	createSelfEnergyData();		
	return 0;
}
