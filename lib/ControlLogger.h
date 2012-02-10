#ifndef ControlLogger_included
#define ControlLogger_included

#include <stdlib.h>

#include "Global.h"
#include "Complex.h"
#include "Quat.h"
#include "ComplexVector.h"
#include "ComplexMatrix.h"



class ControlLogger {
protected:
  void ini();
  void desini();
  
  double PLOT_PointSize;
  double PLOT_FrameSizeX;
  double PLOT_FrameSizeY;
  
  char* latexFileName;
  int plotCounter;
  bool loggingActive;
  void makePlotWithErrors(char* &fName, int N, double* x, double* y, double* err, char* xlabel, char* ylabel, char* title, char* fitCommand);
  void makePlotWithErrorsAndAdditionalLinePlots(char* &fName, int N, double* x, double* y1, double* err1, double* y2, double* err2, double* l1, double* l2, char* xlabel, char* ylabel, char* title1, char* title2, char* fitCommand);
  
  void addPlotTOTEX(char* eps1, char* caption, int angle1, double w1);
  void addDoublePlotTOTEX(char* eps1, char* eps2, char* caption, int angle1, double w1, int angle2, double w2);
  
  
public:
  ControlLogger();
  ~ControlLogger();
  
  void generateLatex();
  void setLatexFileName(char* fname);

  void addSection(const char* caption);
  void clearPage();
  void newPage();
  void addPlot(const char* title1, const char* title2, const char* caption, const char* xlabel, const char* ylabel, double* x, double* y1, double* err1, double* y2, double* err2, int N, const char* fitCommand);
  void addPlotWithAdditionalLinePlots(char* title1, char* title2, char* caption, char* xlabel, char* ylabel, double* x, double* y1, double* err1, double* y2, double* err2, double* l1, double* l2, int N, char* fitCommand1, char* fitCommand2);
  void addGnuplotTableLinePlot(char* title1, char* title2, char* caption, char* xlabel, char* ylabel, double* x, double* y1, double* y2, int N, char* fitCommand1, char* fitCommand2);
  void addTablePlot(char* title1, char* title2, char* caption, char* xlabel1, char* xlabel2, char* ylabel1, char* ylabel2, double* x1, double* x2, double* y1, double* y2, double* yErr1, double* yErr2, int N1, int N2, char* fitCommand1, char* fitCommand2);
  void setLogging(bool on);
  void addDirectText(char* text);
  void setPLOT_PointSize(double size);
  void setPLOT_FrameSize(double x, double y);
};


#endif
