#ifndef EvaluateObservable_included
#define EvaluateObservable_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "Global.h"
#include "AnalyzerIOControl.h"
#include "StateDescriptorReader.h"
#include "LatexAndPlottingSystem.h"
#include "SimulationParameterSet.h"


class EvaluateObservable {
private: 
  AnalyzerIOControl* ioControl;
  char* nickName;
  char* obsName;
  void loadData(int& count, int* &IDs, double** &data);  
  void deleteLoadedData(int& count, int* &IDs, double** &data);
  bool detSignAvail;
  bool weightAvail;
  bool physicalScaleIncludesRenormZFactor;  
  char** xmlOutput_Tag;
  char** xmlOutput_Description;
  double* xmlOutput_Value;
  double* xmlOutput_ValueError;  
  EvaluateObservable** otherObs;
  int otherObscount;
  EvaluateObservable** dependOnObs;
  int dependOnObsCount;
  int xmlOutput_Count;
  int bareObsDataSetsAvailCount;
  int bareObsDataSetsGapCount;
  double relativeEvaluationStartIndex;
  double relativeEvaluationEndIndex;
  void addXmlOutput(char* tag, char* des, double val, double err);
  

protected:
  StateDescriptorReader* SDReader;
  LatexAndPlottingSystem* LAPsystem;
  double physicalScaleInGEV;
  double physicalScaleErrorInGEV;
  int dataAvailCount;
  double** dataAvail;
  int* dataAvailID;
  double* dataAvailWeightAndSign;
  SimulationParameterSet* SimParaSet;
  
  void ini(int resCount, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign);
  void addDependOnObsByName(char* obsName);
  EvaluateObservable* getDependOnObsByName(char* name);
  double* getDataForID(int ID);
  void startLatexOutputSummaryTable();
  void addXMLonly(char* xmltag, char* des, double val, double error);
  void addXML_And_LatexOutputSummaryTableLine(char* xmltag, char* des, char* shortCut, double val, double error, char* unit, char* valFormat);
  void endLatexOutputSummaryTable();
  
 
public:
  EvaluateObservable(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, char* oName, char* nick, double relStart, double relEnd); 
  virtual ~EvaluateObservable();
  
  bool isNick(char* nick);
  bool isObs(char* obs);  
  char* getObsName();
  char* getNickName();
  void setPhysicalScale(double physScale, double physScaleErr, bool ZFacDet);
  void getWeights(int &count, int* &IDs, double* &weigths);
  void getDetSigns(int &count, int* &IDs, double* &detSigns);
  void getData(int &count, int* &IDs, int &anaResCount, double** &data);
  void getDataIDs(int &count, int* &IDs);
  void writeXmlOutput();
  bool evaluateWrapper();
  void generateLatexAndPlotsAndXMLWrapper();  
  void printStatusSummaryToScreen();
  void addOtherObsPointer(EvaluateObservable* obs);
  bool reduceDataDueToDependencies();

  
  virtual bool evaluate() = 0;
  virtual void generateLatexAndPlotsAndXML() = 0;  
  virtual void defineObsDependencies() = 0;    
  virtual int getAnalyzerResultsCount() = 0;
};


#include "EvaluateObservable.C"

#endif
