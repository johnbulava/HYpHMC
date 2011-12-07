#include "EvaluateObservableDetSign.h"

EvaluateObservableDetSign::EvaluateObservableDetSign(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, double relStart, double relEnd) : EvaluateObservable(aIOcon, sdr, "DetSign", "det", relStart, relEnd) { 
  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
  neg_det_count = 0;
}


EvaluateObservableDetSign::~EvaluateObservableDetSign() {
}


int EvaluateObservableDetSign::getAnalyzerResultsCount() {
  return 3;
}

void EvaluateObservableDetSign::defineObsDependencies() { 
}


bool EvaluateObservableDetSign::evaluate() {
	bool res = false;
	printf("Evaluation ...\n");
	if (dataAvailCount > 0 && dataAvail != NULL) {
		for (int I=0; I<dataAvailCount; I++) {
			printf("ID:%d  Sign is %+1.0f ", dataAvailID[I], dataAvail[I][0]);
			if (dataAvail[I][0] < 0.0)
				neg_det_count++;
		}
		res = true;
	}
				
	return res;
}


void EvaluateObservableDetSign::generateLatexAndPlotsAndXML() {
	if (LAPsystem != NULL) {
		LAPsystemPlot* plot = LAPsystem->createNewPlot("plot_det_sgn");		
		int dataLineCount = dataAvailCount;
		int dataColumnCount = getAnalyzerResultsCount()+1;

		double** plotData = new double*[dataLineCount];
		for (int i = 0; i < dataLineCount; i++) {
			plotData[i] = new double[dataColumnCount];
			plotData[i][0] = dataAvailID[i];
			for (int j = 1; j < dataColumnCount; j++) {				
				plotData[i][j] = dataAvail[i][j-1];
			}
		}		
		
		plot->setPlotData(dataLineCount, dataColumnCount, plotData);
				
		for (int i = 0; i < dataLineCount; i++) {
			delete[] plotData[i];
			plotData[i] = NULL;
		}
		delete[] plotData;
		plotData = NULL;
		
		char* str = new char[1000000];						
		plot->setXLabel("Conf ID");		
		plot->setYLabel("Sign Det M");				 
		sprintf(str, "Sign (Det M) evaluated on %d configurations, Lattice: $%d \\times %d \\times %d \\times %d, negDetCount: %d$", dataLineCount, SDReader->getL0(), SDReader->getL1(), SDReader->getL2(), SDReader->getL3(), neg_det_count);				
		plot->setCaption(str);
		plot->setPlotTitle("Sign of determinant of M");
		plot->setXLogScale(false);
		plot->setYLogScale(false);
		plot->setXRange(0.0, INFINITY);
		plot->setYRange(-1.5, 1.5);
		plot->setXErrorBars(false);
		plot->setYErrorBars(false);  
		plot->setSize(0.7, 0.7);
  		plot->setTitle("Det");
  		
		sprintf(str, "1:2");
		plot->plotData(str);						
		LAPsystem->addPlot(plot);
				
		startLatexOutputSummaryTable();

		addXML_And_LatexOutputSummaryTableLine("negDetCount", "Number of negative determinants", "Sgn( det(M) )", neg_det_count, 0, NULL, "%+1.0f");
// 		addXML_And_LatexOutputSummaryTableLine("autoCorVEV", "Auto-correlation time of vev", "$\\tau_{int}^{vev}$", autoCorMtime, 0, NULL, "%1.1f");  

		endLatexOutputSummaryTable();
		
		delete[] str;		
	}
}
