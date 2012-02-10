#include "LatexAndPlottingSystem.h"

LatexAndPlottingSystem::LatexAndPlottingSystem(const char* baseDirName, const char* lapBaseName) { 
	if (LogLevel>1) printf("Initializing LatexAndPlottingSystem with BaseDirectory=%s and LAPsystem BaseName=%s ...\n",baseDirName,lapBaseName);
	BaseDirectoryName = new char[1+strlen(baseDirName)];
	snprintf(BaseDirectoryName, 1+strlen(baseDirName), "%s", baseDirName);
	LAPsystemBaseName = new char[1+strlen(lapBaseName)];
	snprintf(LAPsystemBaseName, 1+strlen(lapBaseName), "%s", lapBaseName);
	
	plotCount = 0;
	plotList = new LAPsystemPlot*[100000];
	for (int i = 0; i < 100000; i++) {
   	plotList[i] = NULL;
	}
		
	latexStr = new char[1000000];
	sprintf(latexStr, "%s", "");
}


LatexAndPlottingSystem::~LatexAndPlottingSystem() {	
	processContent();
	delete[] BaseDirectoryName;
	delete[] LAPsystemBaseName;
	
	
	for (int i = 0; i < 100000; i++) {
		delete plotList[i];
		plotList[i] = NULL;
	}
	delete[] plotList;
	plotList = NULL;	
	delete[] latexStr;
	latexStr = NULL;
}

LAPsystemPlot* LatexAndPlottingSystem::createNewPlot(const char* plotScriptPrefix) {
	LAPsystemPlot* plot = NULL;		
	
	if (plotScriptPrefix!= NULL && strlen(plotScriptPrefix) != 0) {
		plotList[plotCount] = new LAPsystemPlot(LAPsystemBaseName, plotScriptPrefix);
		plot = plotList[plotCount];
		plotCount++; 
	}	
	return plot;
}

void LatexAndPlottingSystem::addPlot(LAPsystemPlot* plot) {
	char* strBuffer = new char[100000];		
	char* captionStr = NULL;
	char* epsFileNamePrefix = NULL;
	
	epsFileNamePrefix = plot->getPlotScriptFileNamePrefix();
	captionStr = plot->getCaption();
	strcat(latexStr, "\n\\begin{figure}[h]\n");
	strcat(latexStr, "\\centering\n");			
		sprintf(strBuffer, "\\includegraphics{%s.eps}\n", epsFileNamePrefix);
	strcat(latexStr, strBuffer);
		sprintf(strBuffer, "\\caption{%s}\n", captionStr);
	strcat(latexStr, strBuffer);
	strcat(latexStr, "\\end{figure}\n\n");
			
	delete[] epsFileNamePrefix;
	delete[] strBuffer;
	delete[] captionStr;
	captionStr = NULL;
	epsFileNamePrefix = NULL;		
}	

void LatexAndPlottingSystem::addSection(char* caption) {
	char* str = new char[strlen(caption)+100];
	sprintf(str, "\n\\section{%s}\n", caption);
	strcat(latexStr, str);
	delete[] str;
}

void LatexAndPlottingSystem::clearPage() {

}

void LatexAndPlottingSystem::newPage() {
	strcat(latexStr, "\\newpage\n");
}

void LatexAndPlottingSystem::addDirectText(const char* text) {
	strcat(latexStr, text);
}

void LatexAndPlottingSystem::processContent() {
//	1.) Write plot data and GNU-Plot scripts for each plot object
// 2.) Create latex body for all plots in latexPlotList
// 	 Create latex header file which needs the body file
//		 Create corresponding Makefile

//	1.)
	for (int i = 0; i < plotCount; i++) {
		plotList[i]->writePlotFilesToDisk(BaseDirectoryName);
	}
	
//	2.)
	createPlotsMakefile(plotList, plotCount);
	char* latexBodyFileName = createLatexBodyFile();
	char* latexHeaderFileName = createLatexHeaderFile(latexBodyFileName);
	createMakefile(latexHeaderFileName);
	
	delete[] latexBodyFileName;
	delete[] latexHeaderFileName;
}

char* LatexAndPlottingSystem::createLatexBodyFile() {
	char* latexBodyFileName = NULL;
	char* strBuffer = new char[1000000];	
	
	latexBodyFileName = new char[100000];
	sprintf(latexBodyFileName, "%s_body", LAPsystemBaseName);
	sprintf(strBuffer, "%s/%s.tex", BaseDirectoryName, latexBodyFileName);
	
	if (!fileExists(strBuffer) || strlen(latexStr) > 0) {
		FILE* latexBodyFile = fopen(strBuffer, "w");
		fprintf(latexBodyFile, "%s", latexStr);		
		fclose(latexBodyFile);			
	}
	
	delete[] strBuffer;
	strBuffer = NULL;
	return latexBodyFileName;
}

char* LatexAndPlottingSystem::createLatexHeaderFile(char* latexBodyFileName) {
	char* latexHeadFileName = NULL;
	char* strBuffer = new char[2048];
	char* latexHead = new char[1000000];
		
	sprintf(latexHead, "\\documentclass[11pt, a4paper, oneside]{article}\n");
	strcat(latexHead, "\\usepackage[latin1]{inputenc}\n");
	strcat(latexHead, "\\usepackage{german,amsfonts}\n");
	strcat(latexHead, "\\usepackage{graphicx, graphics}\n");
	strcat(latexHead, "\\graphicspath{{plots/}}\n");
	strcat(latexHead, "\\usepackage{latexsym}\n");
	strcat(latexHead, "\\usepackage{mathrsfs}\n\n");
	strcat(latexHead, "\\tolerance 1414\n");
	strcat(latexHead, "\\hbadness 1414\n");
	strcat(latexHead, "\\emergencystretch 1.5em\n");
	strcat(latexHead, "\\hfuzz 0.3pt\n");
	strcat(latexHead, "\\widowpenalty=10000\n");
	strcat(latexHead, "\\vfuzz \\hfuzz\n");
	strcat(latexHead, "\\raggedbottom\n");
	strcat(latexHead, "\\pagestyle{plain}\n\n");
	strcat(latexHead, "\\begin{document}\n");
		sprintf(strBuffer, "\\input{%s.tex}\n", latexBodyFileName);
		strcat(latexHead, strBuffer);		
	strcat(latexHead, "\\end{document}\n");	
		
	latexHeadFileName = new char[2048];
	sprintf(latexHeadFileName, "%s", LAPsystemBaseName);
	sprintf(strBuffer, "%s/%s.tex", BaseDirectoryName, latexHeadFileName);
	
	if (!fileExists(strBuffer) || latexBodyFileName != NULL) {
		FILE* latexHeadFile = fopen(strBuffer, "w");
		fprintf(latexHeadFile, "%s", latexHead);		
		fclose(latexHeadFile);
	}		
				
	delete[] strBuffer;
	strBuffer = NULL;		
	delete[] latexHead;
	latexHead = NULL;

	return latexHeadFileName;
}

void LatexAndPlottingSystem::createMakefile(char* latexHeaderFileName) {
	char* strBuffer = new char[2048];
	char* makeFileStr = new char[1000000];
	
	sprintf(makeFileStr, "# Makefile for latex document\n# This file was automatically created by the class LatexAndPlottingSystem\n\n");
		sprintf(strBuffer, "MAIN_FILE=%s\n", latexHeaderFileName);
	strcat(makeFileStr, strBuffer);
	strcat(makeFileStr, "OUTPUT=.\n");
	strcat(makeFileStr, "FORMAT=dvi\n\n");
	strcat(makeFileStr, "all: $(MAIN_FILE).tex\n");
	strcat(makeFileStr, "\tmake images\n");
	strcat(makeFileStr, "\tlatex -interaction=nonstopmode $(MAIN_FILE).tex\n");
	strcat(makeFileStr, "\tdvipdf $(OUTPUT)/$(MAIN_FILE).dvi\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).aux\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).log\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).dvi\n");
	strcat(makeFileStr, "images:\n");
	strcat(makeFileStr, "\tcd plots;make\n\n");
	strcat(makeFileStr, ".PHONY: clean\n");
	strcat(makeFileStr, "clean:\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).aux\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).dvi\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).log\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).ps\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).pdf\n");
	strcat(makeFileStr, "\t-rm ./plots/*.eps\n\n");
	strcat(makeFileStr, "clean_all:\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).aux\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).dvi\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).log\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).ps\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).pdf\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/*.tex\n");
	strcat(makeFileStr, "\t-rm ./plots/*.eps\n");
	strcat(makeFileStr, "\t-rm ./plots/*.dat\n");
	strcat(makeFileStr, "\t-rm ./plots/*.gnu\n");
	strcat(makeFileStr, "\t-rm ./plots/Makefile\n\n");
	
	sprintf(strBuffer, "%s/Makefile", BaseDirectoryName);
	if (!fileExists(strBuffer) || latexHeaderFileName != NULL) {		
		FILE* makefile = fopen(strBuffer, "w");
		fprintf(makefile, "%s", makeFileStr);
		fclose(makefile);
	}
	
	delete[] strBuffer;
	strBuffer = NULL;
	delete[] makeFileStr;
	makeFileStr = NULL;
}

void LatexAndPlottingSystem::createPlotsMakefile(LAPsystemPlot** plots, int plotCount) {
	char* str_1 = new char[1024];
	char* str_2 = new char[1024];
	char* plotMakefileStr = new char[1000000];		
	
	sprintf(plotMakefileStr, "# makefile for images\n");		
	strcat(plotMakefileStr, "# GNU script file names without suffix .gnu\n");
	sprintf(str_2, "# Objects: image files in eps format\nOBJS=");
	char* tmp = NULL;
	for(int i = 0; i < plotCount; i++) {
		tmp = plots[i]->getPlotScriptFileNamePrefix();
		sprintf(str_1, "SCRIPT_FILE_%d=%s\n", i+1, tmp);
		delete[] tmp;
		
		strcat(plotMakefileStr, str_1);
		sprintf(str_1, "$(SCRIPT_FILE_%d).eps ", i+1);
		strcat(str_2, str_1);
	}
	strcat(str_2,"\n\n");		
	strcat(plotMakefileStr, str_2);				
	strcat(plotMakefileStr, "all: $(OBJS)\n\n");
	for(int i = 0; i < plotCount; i++) {
		sprintf(str_1, "$(SCRIPT_FILE_%d).eps: $(SCRIPT_FILE_%d).gnu\n", i+1, i+1);
		strcat(plotMakefileStr, str_1);
		sprintf(str_1, "\t./makePlot -bw $(SCRIPT_FILE_%d)\n", i+1);
		strcat(plotMakefileStr, str_1);				
	}
	strcat(plotMakefileStr, "\n\n");
	strcat(plotMakefileStr, ".PHONY: clean\n");
	strcat(plotMakefileStr, "clean:\n");
	strcat(plotMakefileStr, "\t-rm $(OBJS)\n");
	
	sprintf(str_1, "%s/plots/Makefile", BaseDirectoryName);
	if (!fileExists(str_1) || plotCount > 0) {		
		FILE* plotMakefile = fopen(str_1, "w");
		fprintf(plotMakefile, "%s", plotMakefileStr);
		fclose(plotMakefile);
	}
	
	delete[] str_1;
	str_1 = NULL;
	delete[] str_2;
	str_2 = NULL;
	delete[] plotMakefileStr;
	plotMakefileStr = NULL;
}

bool LatexAndPlottingSystem::fileExists(std::string strFilename) {
	struct stat stFileInfo;
	bool blnReturn;
	int intStat;

  // Attempt to get the file attributes
	intStat = stat(strFilename.c_str(),&stFileInfo);
	if(intStat == 0) {
    // We were able to get the file attributes
    // so the file obviously exists.
		blnReturn = true;
	} else {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
		blnReturn = false;
	}
  
	return(blnReturn);
}
