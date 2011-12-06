#ifndef LatexAndPlottingSystem_included
#define LatexAndPlottingSystem_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>

#include "Global.h"
#include "LAPsystemPlot.h"


class LatexAndPlottingSystem {
private:  
	char* BaseDirectoryName;
	char* LAPsystemBaseName;
	LAPsystemPlot** plotList;
	int plotCount;		
	char* latexStr;		
	
	char* createLatexBodyFile(); // returns the latex body file name PREFIX!
	char* createLatexHeaderFile(char* latexBodyFileName); // returns the latex header file name PREFIX!
	void createMakefile(char* latexHeaderFileName); 
	void createPlotsMakefile(LAPsystemPlot** plots, int plotCount);
	bool fileExists(std::string strFilename);
		
public:	    
  LatexAndPlottingSystem(const char* baseDirName, const char* lapBaseName); 
  ~LatexAndPlottingSystem();  // Destructor creates Makefiles (1. plots, 2. latex-file) and writes Latex-Body and writes Latex-Wrapper and so on.

  LAPsystemPlot* createNewPlot(const char* plotScriptPrefix);
  void addPlot(LAPsystemPlot* plot);
  void addSection(char* caption);
  void clearPage();
  void newPage();
  void addDirectText(char* text);
  void processContent();
};


#include "LatexAndPlottingSystem.C"

#endif
