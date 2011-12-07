#ifndef AnalyzerIOControl_included
#define AnalyzerIOControl_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "Global.h"
#include "Tools.h"



class AnalyzerIOControl {
private:
  char* fileNameExtension;
  int JobID;

public:
  AnalyzerIOControl(char* extension, int jID); 
  ~AnalyzerIOControl();
  
  void createWorkFolder();
  void createObsFolder(char* obsName);  
  void createAndPrepareObsPlotsFolder(char* obsName);    
  void markAsInProgress(int ID);
  void unmarkAsInProgress(int ID);
  bool isInProgress(int ID);
  void getConfInProgressIDs(int& count, int* &confInProgIDs);
  void removeDeprecatedInProgressFiles(double maxAgeinHours);
  char* getWorkFolderName(); //creates a new char*-String. Must be deleted by calling function
  char* getObservableFolderName(char* obsName);  
  char* getObservableOutsourcedFolderName(char* obsName);  
  char* getObservablePlotsFolderName(char* obsName);    
  char* getObservableFileName(char* obsName);
  char* getObservableOutsourceFileName(char* obsName, int confID);
  char* getObservableOutsourceFileName(char* obsName, int jID, int confID);  
  char* getObservableXmlFileName(char* obsName);    
  char* getObservableLatexBaseName(char* obsName);
  char* getObservableLatexBodyLocalFileName(char* obsName);
  char* getObservableLatexBodyFileName(char* obsName);  
  char* getLatexSummaryBaseName();
  char* getStateDescriptorFileName();
  char* getMainXmlFileName();  
  char* getPhiConfFileName(int ID);
  int getTotalPhiConfCount();
  void getObservableFileNameList(char* obsName, char** &fileNames, int& fileCount);
  int getJobID();
};


#endif
