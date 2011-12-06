LatexSummary::LatexSummary(char* latexBaseFilename, char* baseDirName, int subDirCount, char** subDirNames, char** localLatexBodyFileNames) {
	printf("\n\n\n    ---   LATEX SUMMARY   --- \n\n");
	printf("latexBaseFilename: %s\n",latexBaseFilename);
	printf("baseDirName: %s\n",baseDirName);
	printf("Sub directories:\n");
	for (int I=0; I<subDirCount; I++) {
		printf("%d. %s ... %s\n", I, subDirNames[I], localLatexBodyFileNames[I]);  
	}
	this->latexBaseFilename = new char[strlen(latexBaseFilename)+1];
		strcpy(this->latexBaseFilename, latexBaseFilename);
	this->baseDirName = new char[strlen(baseDirName)+1];
		strcpy(this->baseDirName, baseDirName);
	this->subDirCount = subDirCount;
	this->subDirNames = new char*[subDirCount];
	this->localLatexBodyFileNames = new char*[subDirCount];
	for (int i = 0; i < subDirCount; i++) {
		this->subDirNames[i] = new char[strlen(subDirNames[i])+1];
		strcpy(this->subDirNames[i], subDirNames[i]);
					
		this->localLatexBodyFileNames[i] = new char[strlen(localLatexBodyFileNames[i])+1];
		strcpy(this->localLatexBodyFileNames[i], localLatexBodyFileNames[i]);
	}	
	preIncludeText = NULL;
	postIncludeText = NULL;
	// ACHTUNG: Möglicherweise ist body-datei nicht existent  ==> checken
}


LatexSummary::~LatexSummary() {
	char* bodyFilePrefix = createLatexBodyFile();
	char* headFilePrefix = createLatexHeaderFile(bodyFilePrefix);
	createMakefile(headFilePrefix);
	delete[] bodyFilePrefix;
	delete[] headFilePrefix;
	
	delete[] latexBaseFilename;
	delete[] baseDirName;
	for (int i = 0; i < subDirCount; i++) {
		delete[] subDirNames[i];
		delete[] localLatexBodyFileNames[i];
	}
	delete[] subDirNames;
	delete[] localLatexBodyFileNames;
	delete[] preIncludeText;
	delete[] postIncludeText; 
}


void LatexSummary::addDirectTextBeforeIncludes(char* text) {
	printf("addDirectBefore: %s\n",text);
	if(text != NULL && strlen(text) != 0) {
		if (preIncludeText != NULL && strlen(preIncludeText) != 0) {
			char* tmp = new char[strlen(preIncludeText)+1];
			strcpy(tmp, preIncludeText);
			delete[] preIncludeText;
			preIncludeText = new char[strlen(text) + strlen(tmp) + 1];
			strcpy(preIncludeText, tmp);
			delete[] tmp;
			strcat(preIncludeText, text);
		}
		else {
			preIncludeText = new char[strlen(text) + 1];
			strcpy(preIncludeText, text);
		}			
	}	
}


void LatexSummary::addDirectTextAfterIncludes(char* text) {
  printf("addDirectAfter: %s\n",text);  
  if(text != NULL && strlen(text) != 0) {
	  if (postIncludeText != NULL && strlen(postIncludeText) != 0) {
		  char* tmp = new char[strlen(postIncludeText)+1];
		  strcpy(tmp, postIncludeText);
		  delete[] postIncludeText;
		  postIncludeText = new char[strlen(text) + strlen(tmp) + 1];
		  strcpy(postIncludeText, tmp);
		  delete[] tmp;
		  strcat(postIncludeText, text);
	  }
	  else {
		  postIncludeText = new char[strlen(text) + 1];
		  strcpy(postIncludeText, text);
	  }			
  }  
}

bool LatexSummary::fileExists(std::string strFilename) {
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

char* LatexSummary::createLatexBodyFile() {
	char* subBodyFileName = new char[1000000];
	char* latexBodyContent = new char[10000000];
	char* strBuffer = new char[100000];
	char* bodyFileName = new char[100000];
	
	bool newContent = false;
	
	if(preIncludeText != NULL && strlen(preIncludeText) != 0) {
		newContent = true;
		sprintf(latexBodyContent, "%s\n\n", preIncludeText);
	}
	else 
		sprintf(latexBodyContent, "%s", "");	 
	
	for (int i = 0; i < subDirCount; i++) {
		sprintf(subBodyFileName, "%s/%s/%s", baseDirName, subDirNames[i], localLatexBodyFileNames[i]);
		if (fileExists(subBodyFileName)) {
			newContent = true;
			sprintf(strBuffer, "\\input{%s/%s}\n", subDirNames[i], localLatexBodyFileNames[i]);
			strcat(latexBodyContent, strBuffer);
			strcat(latexBodyContent, "\\clearpage\n");
		}
	}
	if(postIncludeText != NULL && strlen(postIncludeText) != 0) {
		newContent = true;
		strcat(latexBodyContent, postIncludeText);
	}
	
	
	sprintf(bodyFileName, "%s_body", latexBaseFilename);
	sprintf(strBuffer, "%s/%s.tex", baseDirName, bodyFileName);
	if (!fileExists(strBuffer) || newContent) {		
		FILE* latexBodyFile = fopen(strBuffer, "w");
		fprintf(latexBodyFile, "%s", latexBodyContent);		
		fclose(latexBodyFile);
	}	
	
	delete[] subBodyFileName;
	delete[] latexBodyContent;
	delete[] strBuffer;
	return bodyFileName;
}

char* LatexSummary::createLatexHeaderFile(char* latexBodyFilePrefix) {
	char* strBuffer = new char[1000000];
	char* latexHead = new char[10000000];
	char* tmpCharArr = new char[100000];
	char* headerFilePrefix = new char[100000];
		
	sprintf(latexHead, "\\documentclass[11pt, a4paper, oneside]{article}\n");
	strcat(latexHead, "\\usepackage[latin1]{inputenc}\n");
	strcat(latexHead, "\\usepackage{german,amsfonts}\n");
	strcat(latexHead, "\\usepackage{graphicx, graphics}\n");
	
	strcat(latexHead, "\\graphicspath{");
	for (int i = 0; i < subDirCount; i++) {
		sprintf(tmpCharArr, "%s/%s/%s", baseDirName, subDirNames[i], localLatexBodyFileNames[i]);
		if (fileExists(tmpCharArr)) {
			sprintf(strBuffer, "{%s/plots/} ", subDirNames[i]);
			strcat(latexHead, strBuffer);
		}
	}
	strcat(latexHead, "}\n");
	delete[] tmpCharArr;
	
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
		
	sprintf(strBuffer, "\\input{%s.tex}\n", latexBodyFilePrefix);
	strcat(latexHead, strBuffer);
	strcat(latexHead, "\\clearpage\n");
	
	strcat(latexHead, "\\end{document}\n");	
	if(latexHead != NULL && strlen(latexHead) != 0) {	
		sprintf(headerFilePrefix, "%s", latexBaseFilename);
		sprintf(strBuffer, "%s/%s.tex", baseDirName, headerFilePrefix);		
		FILE* latexHeadFile = fopen(strBuffer, "w");
		fprintf(latexHeadFile, "%s", latexHead);		
		fclose(latexHeadFile);		
	}
			
	delete[] strBuffer;
	strBuffer = NULL;		
	delete[] latexHead;
	latexHead = NULL;
	return headerFilePrefix;
}

void LatexSummary::createMakefile(char* latexHeadFilePrefix) {
	char* strBuffer = new char[1000000];
	char* makeFileStr = new char[10000000];
	
	sprintf(makeFileStr, "# Makefile for latex document\n# This file was automatically created by the class LatexAndPlottingSystem\n\n");
	sprintf(strBuffer, "MAIN_FILE=%s\n", latexHeadFilePrefix);
	strcat(makeFileStr, strBuffer);
	
	strcat(makeFileStr, "SUB_FILES=");
	
	char* tmpCharArr = new char[100000];
	for(int i = 0; i < subDirCount; i++) {
		sprintf(tmpCharArr, "%s/%s/%s", baseDirName, subDirNames[i], localLatexBodyFileNames[i]);
		if (fileExists(tmpCharArr)) {
			sprintf(strBuffer, "%s/%s ", subDirNames[i], localLatexBodyFileNames[i]);
			strcat(makeFileStr, strBuffer);
		}
	}
	delete[] tmpCharArr;
	strcat(makeFileStr, "\n");
	
	strcat(makeFileStr, "OUTPUT=.\n");
	strcat(makeFileStr, "FORMAT=dvi\n\n");
	strcat(makeFileStr, "all: $(MAIN_FILE).tex $(SUB_FILES)\n");
	strcat(makeFileStr, "\tmake images\n");
	strcat(makeFileStr, "\tlatex -interaction=nonstopmode $(MAIN_FILE).tex\n");
	strcat(makeFileStr, "\tdvipdf $(OUTPUT)/$(MAIN_FILE).dvi\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).aux\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).log\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).dvi\n");
	strcat(makeFileStr, "images:\n");
	for(int i = 0; i < subDirCount; i++) {
		sprintf(strBuffer, "\tcd %s;make\n", subDirNames[i]);
		strcat(makeFileStr, strBuffer);		
	}	
	strcat(makeFileStr, "\n");
			
	strcat(makeFileStr, ".PHONY: clean\n");
	strcat(makeFileStr, "clean:\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).aux\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).dvi\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).log\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).pdf\n");
	for(int i = 0; i < subDirCount; i++) {
		sprintf(strBuffer, "\tcd %s;make clean\n", subDirNames[i]);
		strcat(makeFileStr, strBuffer);		
	}	
	strcat(makeFileStr, "\n");
	
	strcat(makeFileStr, "clean_all:\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).aux\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).dvi\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).log\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).pdf\n");
	strcat(makeFileStr, "\t-rm $(OUTPUT)/$(MAIN_FILE).tex\n");
	for(int i = 0; i < subDirCount; i++) {
		sprintf(strBuffer, "\tcd %s;make clean_all\n", subDirNames[i]);
		strcat(makeFileStr, strBuffer);		
	}	
	strcat(makeFileStr, "\n");
	
	sprintf(strBuffer, "%s/Makefile", baseDirName);		
	FILE* makefile = fopen(strBuffer, "w");
	fprintf(makefile, "%s", makeFileStr);
	fclose(makefile);
	
	delete[] strBuffer;
	strBuffer = NULL;
	delete[] makeFileStr;
	makeFileStr = NULL;
}
