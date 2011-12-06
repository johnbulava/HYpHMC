#ifndef LAPsystemPlot_included
#define LAPsystemPlot_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>
#include <cstring>

#include "Global.h"


class LAPsystemPlot {
private:  
	char* PlotName;
	char* BaseName;
	double** dataLines;
	int dataLineCount;
	int dataColCount;
	char* plotCommand;
	char* plotSubCommand;
	char* plotTitle;
	
	bool xLogScale, yLogScale;
	char* xLabel;
	char* yLabel;
	char* caption;
	double xSize, ySize;
	double pointSize;
	int lineWidth;
	double xRangeMin, xRangeMax;
	double yRangeMin, yRangeMax;
	bool xErrorBars, yErrorBars;
	int pointType, lineType;
	char* title;
	char* usingStr;
	int plotCmdCount; // to distinguish between plot and replot
	
	bool writeDataToFile(char* fileName);
	void appendEndLine(char* str);
	void str_replace(char* str, const char* token, const char* replacement);	

public:    
  LAPsystemPlot(const char* baseName, const char* plName); 
  ~LAPsystemPlot();

  //Called from LatexAndPlottingSystem
  void writePlotFilesToDisk(char* baseDirectoryName);
  
  //Called from user
  void setPlotData(int datLinCnt, int datColCnt, double** dat);
    
  void setCaption(const char* label);
  void setXLabel(const char* label);
  void setYLabel(const char* label);
  void setTitle(const char* title);    
  void setPlotTitle(const char* title);
  void setXLogScale(bool on);
  void setYLogScale(bool on);
  void setXErrorBars(bool on);
  void setYErrorBars(bool on);
  void setSize(double x, double y);
  void setPointSize(double size);  
  void setXRange(double min, double max);
  void setYRange(double min, double max);      
  void setPointType(int type);
  void setLineType(int type);
  void setLineWidth(int width);
  
  void plotData(const char* usingStr);
  void plotDirect(const char* plotCmd);
  
  char* getPlotScriptFileNamePrefix();
  char* getCaption();
};

#endif
