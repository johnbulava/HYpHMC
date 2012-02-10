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
  void addXmlOutput(const char* tag, const char* des, double val, double err);
  

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
  void addDependOnObsByName(const char* obsName);
  EvaluateObservable* getDependOnObsByName(const char* name);
  double* getDataForID(int ID);
  void startLatexOutputSummaryTable();
  void addXMLonly(const char* xmltag, const char* des, double val, double error);
  void addXML_And_LatexOutputSummaryTableLine(const char* xmltag, const char* des, const char* shortCut, double val, double error, const char* unit, const char* valFormat);
  void endLatexOutputSummaryTable();
  
 
public:
  EvaluateObservable(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, const char* oName, const char* nick, double relStart, double relEnd); 
  virtual ~EvaluateObservable();
  
  bool isNick(char* nick);
  bool isObs(const char* obs);  
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


#endif
