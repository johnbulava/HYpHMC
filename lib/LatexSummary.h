#ifndef LatexSummary_included
#define LatexSummary_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>
#include <string>
#include <vector>

#include "Global.h"


class LatexSummary {
private:
	char* latexBaseFilename;
	char* baseDirName;
	int subDirCount;
	char** subDirNames;
	char** localLatexBodyFileNames;
	char* preIncludeText;
	char* postIncludeText;
	
	bool fileExists(std::string strFilename);
	char* createLatexBodyFile();	
	char* createLatexHeaderFile(char* latexBodyFilePrefix);  
	void createMakefile(char* latexHeadFilePrefix);

public:    
  LatexSummary(char* latexBaseFilename, char* baseDirName, int subDirCount, char** subDirNames, char** localLatexBodyFileNames); 
  ~LatexSummary();  // Body- und Head-Files werden geschrieben, wenn Destructor aufgerufen wird.
  
  void addDirectTextBeforeIncludes(char* text);
  void addDirectTextAfterIncludes(char* text);
};


#include "LatexSummary.C"

#endif
