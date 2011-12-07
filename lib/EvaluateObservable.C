#include "EvaluateObservable.h"

EvaluateObservable::EvaluateObservable(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, char* oName, char* nick, double relStart, double relEnd) { 
  if (LogLevel>1) printf("EvaluateObservable %s (%s) initializing with rel. eval-indices (%1.2f/%1.2f)...\n",oName,nick, relStart, relEnd);
  ioControl = aIOcon;
  SDReader = sdr;
  nickName = new char[1+strlen(nick)];
  snprintf(nickName, 1+strlen(nick), "%s", nick);
  obsName = new char[1+strlen(oName)];
  snprintf(obsName, 1+strlen(oName), "%s", oName);  
  relativeEvaluationStartIndex = relStart;
  relativeEvaluationEndIndex = relEnd;

  physicalScaleInGEV = NaN;
  physicalScaleErrorInGEV = NaN;
  dataAvailCount = 0;
  dataAvail = NULL;
  dataAvailID = NULL;
  dataAvailWeightAndSign = NULL;
  detSignAvail = false;
  weightAvail = false;
  physicalScaleIncludesRenormZFactor = false;
  bareObsDataSetsAvailCount = 0;
  bareObsDataSetsGapCount = 0;

  otherObs = new EvaluateObservable*[10000];
  otherObscount = 0;
  dependOnObs = new EvaluateObservable*[10000];
  dependOnObsCount = 0;

  SimParaSet = new SimulationParameterSet(SDReader, SimulationParameterSet_NfNotation);
  
  xmlOutput_Tag = new char*[10000];
  xmlOutput_Description = new char*[10000];
  xmlOutput_Value = new double[10000];
  xmlOutput_ValueError = new double[10000]; 
  xmlOutput_Count = 0;
  
  char* baseDir = ioControl->getObservableFolderName(obsName);
  char* fileNameExtension = SDReader->getFileNameExtension();
  char* latexFileName = ioControl->getObservableLatexBaseName(obsName);
  LAPsystem = new LatexAndPlottingSystem(baseDir, latexFileName);
  delete[] baseDir;
  delete[] fileNameExtension;
  delete[] latexFileName;
}


EvaluateObservable::~EvaluateObservable() {
  if (LogLevel>1) printf("Desinitializing Observable %s...",obsName);
  
  deleteLoadedData(dataAvailCount, dataAvailID, dataAvail);
  delete[] dataAvailWeightAndSign;
  delete LAPsystem;
  delete[] nickName;
  delete[] obsName;  
  for (int I=0; I<xmlOutput_Count; I++) {
    delete[] xmlOutput_Tag[I];
    delete[] xmlOutput_Description[I];
  }
  delete[] xmlOutput_Value;
  delete[] xmlOutput_ValueError;
  delete[] xmlOutput_Tag;
  delete[] xmlOutput_Description;
  xmlOutput_Count = 0;
  delete[] otherObs;
  delete[] dependOnObs;
  delete SimParaSet;
  
  if (LogLevel>1) printf("sucessfully\n");
}


void EvaluateObservable::ini(int resCount, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign) {
  if (LogLevel>2) printf("Loading Data with %d entries per data line...\n", resCount);

  int count = 0;
  int* IDs = NULL;
  double** data = NULL;
  loadData(count, IDs, data);
  if (LogLevel>2) printf("...%d data lines read for observable %s.\n",count,obsName);
  bareObsDataSetsAvailCount = count;
  bareObsDataSetsGapCount = 0;
  if (count>0) {
    bareObsDataSetsGapCount += IDs[0]-1;
    for (int I=1; I<count; I++) {
      bareObsDataSetsGapCount += IDs[I]-IDs[I-1]-1;  
    }
  }

  int obsWeightCount=0;
  int* obsWeightIDs = NULL;
  double* obsWeightData = NULL;
  if (obsWeight != NULL) obsWeight->getWeights(obsWeightCount, obsWeightIDs, obsWeightData);
  weightAvail = (obsWeight != NULL);
  if (LogLevel>2) {
    if (obsWeight!=NULL) {
      printf("...%d data lines read for observable Weight.\n",obsWeightCount);
    } else {
      printf("No reference to observable Weight specified. Ignoring weigths...\n");    
    }
  }
  
  int obsDetSignCount=0;
  int* obsDetSignIDs = NULL;
  double* obsDetSignData = NULL;
  if (obsDetSign != NULL) obsDetSign->getDetSigns(obsDetSignCount, obsDetSignIDs, obsDetSignData);
  detSignAvail = (obsDetSign != NULL);
  if (LogLevel>2) {
    if (obsDetSign!=NULL) {
      printf("...%d data lines read for observable DetSign.\n",obsDetSignCount);
    } else {
      printf("No reference to observable DetSign specified. Ignoring DetSigns...\n");    
    }
  }
  
  double* localWeightData = new double[count];
  for (int I=0; I<count; I++) localWeightData[I] = NaN;
  int startI2 = 0;
  for (int I=0; I<count; I++) {
    int id1 = IDs[I];
    for (int I2=startI2; I2<obsWeightCount; I2++) {
      int id2 = obsWeightIDs[I2];
      if (id1==id2) {
        localWeightData[I] = obsWeightData[I2];
	startI2 = I2+1;
	break;
      }
      if (id2>id1) break;
    }
  }
  
  double* localDetSignData = new double[count];
  for (int I=0; I<count; I++) localDetSignData[I] = NaN;
  startI2 = 0;
  for (int I=0; I<count; I++) {
    int id1 = IDs[I];
    for (int I2=startI2; I2<obsDetSignCount; I2++) {
      int id2 = obsDetSignIDs[I2];
      if (id1==id2) {
        localDetSignData[I] = obsDetSignData[I2];
	startI2 = I2+1;
	break;
      }
      if (id2>id1) break;
    }
  }

  dataAvailCount = 0;
  for (int I=0; I<count; I++) {
    if (obsWeightData==NULL) {
      if (obsDetSignData==NULL) {
        dataAvailCount++;
      } else {
        if (!isNaN(localDetSignData[I])) dataAvailCount++;
      }
    } else {
      if (obsDetSignData==NULL) {
        if (!isNaN(localWeightData[I])) dataAvailCount++;      
      } else {
        if ((!isNaN(localWeightData[I])) && (!isNaN(localDetSignData[I]))) dataAvailCount++;            
      }
    }
  }
  if (LogLevel>2) printf("...%d data lines ready for evaluation of observable %s.\n",count,obsName);
  
  dataAvail = new double*[dataAvailCount];
  dataAvailID = new int[dataAvailCount];
  dataAvailWeightAndSign = new double[dataAvailCount];

  dataAvailCount = 0;
  for (int I=0; I<count; I++) {
    if (obsWeightData==NULL) {
      if (obsDetSignData==NULL) {
        dataAvail[dataAvailCount] = new double[resCount];
        for (int I2=0; I2<resCount; I2++) dataAvail[dataAvailCount][I2] = data[I][I2];
	dataAvailID[dataAvailCount] = IDs[I];
	dataAvailWeightAndSign[dataAvailCount] = 1.0;
      
        dataAvailCount++;
      } else {
        if (!isNaN(localDetSignData[I])) { 
          dataAvail[dataAvailCount] = new double[resCount];
          for (int I2=0; I2<resCount; I2++) dataAvail[dataAvailCount][I2] = data[I][I2];
	  dataAvailID[dataAvailCount] = IDs[I];
	  dataAvailWeightAndSign[dataAvailCount] = localDetSignData[I];
      
	  dataAvailCount++;
	}
      }
    } else {
      if (obsDetSignData==NULL) {
        if (!isNaN(localWeightData[I])) {
          dataAvail[dataAvailCount] = new double[resCount];
          for (int I2=0; I2<resCount; I2++) dataAvail[dataAvailCount][I2] = data[I][I2];
	  dataAvailID[dataAvailCount] = IDs[I];
	  dataAvailWeightAndSign[dataAvailCount] = localWeightData[I];
	
	  dataAvailCount++;      
	}
      } else {
        if ((!isNaN(localWeightData[I])) && (!isNaN(localDetSignData[I]))) {
          dataAvail[dataAvailCount] = new double[resCount];
          for (int I2=0; I2<resCount; I2++) dataAvail[dataAvailCount][I2] = data[I][I2];
	  dataAvailID[dataAvailCount] = IDs[I];
	  dataAvailWeightAndSign[dataAvailCount] = localWeightData[I] * localDetSignData[I];
	
	  dataAvailCount++;  
	}
      }
    }
  }

  deleteLoadedData(count, IDs, data);  
  delete[] obsWeightIDs;
  delete[] obsWeightData;
  delete[] obsDetSignIDs;
  delete[] obsDetSignData;
  delete[] localWeightData;
  delete[] localDetSignData;
}


void EvaluateObservable::printStatusSummaryToScreen() {
  printf("==> %s\n",obsName);
  printf("   -> Evaluated data sets:     %d\n",dataAvailCount);
  printf("   -> Available data sets:     %d\n",bareObsDataSetsAvailCount);
  printf("   -> Gaps between data sets:  %d\n",bareObsDataSetsGapCount);
  printf("   -> FLAGS (weight, det, Z):  %d,%d,%d\n",weightAvail,detSignAvail,physicalScaleIncludesRenormZFactor);


  printf("\n");
}


bool EvaluateObservable::evaluateWrapper() {
  if (dataAvailCount > 0) {
    return evaluate();
  } else {
    return true;
  }
}


void EvaluateObservable::generateLatexAndPlotsAndXMLWrapper() {
  char* sectionName = new char[1000];
  snprintf(sectionName,1000,"Observable %s:",obsName);
  LAPsystem->addSection(sectionName);
  delete[] sectionName;

  if (dataAvailCount > 0) {
    bool det = isObs("DetSign");
    bool wei = isObs("Weight");
    bool gpr = isObs("GoldstonePropagator");
//    bool mag = isObs("Magnetizations");

    if ((!weightAvail) && (!wei)) {
      LAPsystem->addDirectText("\\noindent WARNING: Weight of configurations not considered\n\n");
    }
    if ((!detSignAvail) && (!det) && (!wei)) {
      LAPsystem->addDirectText("\\noindent WARNING: Sign of fermionic determinant not considered\n\n");
    }
    if ((!physicalScaleIncludesRenormZFactor) && (!gpr) && (!det) && (!wei)) {
      LAPsystem->addDirectText("\\noindent WARNING: Goldstone Renormailzation Factor Z not considered\n\n");
    }
    char* text = new char[1000]; 
    snprintf(text,1000,"\\noindent Evaluated data sets: %d\n\n",dataAvailCount);
    LAPsystem->addDirectText(text);
    if ((relativeEvaluationStartIndex>0) || (relativeEvaluationEndIndex<1)) {
      snprintf(text,1000,"\\noindent WARNING: Only part of available data evaluated (%1.2f/%1.2f)\n\n", relativeEvaluationStartIndex, relativeEvaluationEndIndex);
      LAPsystem->addDirectText(text);
    }
    delete[] text;

    generateLatexAndPlotsAndXML();
  }

  addXmlOutput("evalDataSets", "Evaluated data sets", dataAvailCount, 0);
  addXmlOutput("evalRelStartIndex", "Relative start index for evaluation", relativeEvaluationStartIndex, 0);
  addXmlOutput("evalRelEndIndex", "Relative end index for evaluation", relativeEvaluationEndIndex, 0);  
  addXmlOutput("availDataSets", "Available data sets", bareObsDataSetsAvailCount, 0);
  addXmlOutput("gapsInDataSets", "Gaps between data sets", bareObsDataSetsGapCount, 0);
  addXmlOutput("weightAvail", "Weight of configurations considered?", weightAvail, 0);
  addXmlOutput("detSignAvail", "Sign of fermionic determinant considered?", detSignAvail, 0);
  addXmlOutput("physicalScaleIncludesRenormZFactor", "Goldstone Renormailzation Factor Z considered?", physicalScaleIncludesRenormZFactor, 0);
}


void EvaluateObservable::setPhysicalScale(double physScale, double physScaleErr, bool ZFacDet) {
  physicalScaleInGEV = physScale;
  physicalScaleErrorInGEV = physScaleErr;
  physicalScaleIncludesRenormZFactor = ZFacDet;
}


bool EvaluateObservable::isNick(char* nick) {
  if (strlen(nick) != strlen(nickName)) return false;
  
  for (int I=0; I<(int)strlen(nick); I++) {
    if (nick[I] != nickName[I]) return false;
  }

  return true;
}


bool EvaluateObservable::isObs(char* obs) {
  if (strlen(obs) != strlen(obsName)) return false;
  
  for (int I=0; I<(int)strlen(obs); I++) {
    if (obs[I] != obsName[I]) return false;
  }

  return true;
}


void EvaluateObservable::getWeights(int &count, int* &IDs, double* &weights) {
  count = 0;
  IDs = NULL;
  weights = NULL;
  if (!isObs("Weight")) return;
  
  count = dataAvailCount;
  if (count <= 0) return;
  IDs = new int[count];
  weights = new double[count];
  
  for (int I=0; I<count; I++) {
    IDs[I] = dataAvailID[I];
    weights[I] = dataAvail[I][0];
  }
}


void EvaluateObservable::getDetSigns(int &count, int* &IDs, double* &detSigns) {
  count = 0;
  IDs = NULL;
  detSigns = NULL;
  if (!isObs("DetSign")) return;

  count = dataAvailCount;
  if (count <= 0) return;
  IDs = new int[count];
  detSigns = new double[count];
  for (int I=0; I<count; I++) {
    IDs[I] = dataAvailID[I];
    detSigns[I] = dataAvail[I][0];
  }
}


void EvaluateObservable::getData(int &count, int* &IDs, int &anaResCount, double** &data) {
  count = 0;
  IDs = NULL;
  data = NULL;

  count = dataAvailCount;
  if (count <= 0) return;
  IDs = new int[count];
  data = new double*[count];

  anaResCount = getAnalyzerResultsCount();
  for (int I=0; I<count; I++) {
    IDs[I] = dataAvailID[I];
    data[I] = new double[anaResCount];
    for (int I2=0; I2<anaResCount; I2++) {
      data[I][I2] = dataAvail[I][I2];    
    }
  }
}


void EvaluateObservable::getDataIDs(int &count, int* &IDs) {
  count = 0;
  IDs = NULL;

  count = dataAvailCount;
  if (count <= 0) return;
  IDs = new int[count];
  for (int I=0; I<count; I++) {
    IDs[I] = dataAvailID[I];
  }
}


void EvaluateObservable::loadData(int& count, int* &IDs, double** &data) {
  count = 0;
  IDs = NULL;
  data = NULL;
  char** fileNames = NULL;
  int fileCount = 0;
  int interimDataMAX = 1000000;
  double** interimData = new double*[interimDataMAX];
  int dataSetLength = getAnalyzerResultsCount();

  for (int I=0; I<interimDataMAX; I++) interimData[I] = NULL;
  ioControl->getObservableFileNameList(obsName, fileNames, fileCount);
  
  
  for (int I=0; I<fileCount; I++) {
    FILE* file = fopen(fileNames[I],"r");  
    int id = -1;
    while (fscanf(file,"%d",&id)==1) {
      if (id>=interimDataMAX) {
        printf("ERROR: EvaluateObservable::loadData: ID=%d higher than %d\n",id,interimDataMAX);
	exit(0);
      }
      if (interimData[id]==NULL) {
        interimData[id] = new double[dataSetLength];
      } else {
        if (LogLevel>1) printf("EvaluateObservable::loadData: Obs %s, ConfIF=%d exists twice (at least).\n",obsName, id);
      }
      if (dataSetLength<100) {
        for (int I2=0; I2<dataSetLength; I2++) {
          if (fscanf(file,"%lf",&(interimData[id][I2]))!=1) {
            printf("ERROR: EvaluateObservable::loadData: Could not read datum: ConfIF=%d, dataNr=%d, confFile=%s\n",id,I2,fileNames[I]);
	    exit(0);
   	  }
        } 
      } else {
        int jID = -1;
        if (fscanf(file,"%d",&jID)!=1) {
          printf("ERROR: EvaluateObservable::loadData: Could not read JobID for outsourced data set: ConfIF=%d, confFile=%s\n",id,fileNames[I]);
          exit(0);
	}
        char* outFileName = ioControl->getObservableOutsourceFileName(obsName, jID, id);
        std::fstream outFile;
        outFile.open(outFileName, std::ios::in);
        if (!outFile.good()) {
          printf("ERROR: EvaluateObservable::loadData: Could open outsourced data set: ConfIF=%d, outsourcedFile=%s\n",id,outFileName);
          exit(0);
        }

        outFile.read((char*)(interimData[id]), 8*dataSetLength);
        if (outFile.eof()) {
          printf("ERROR: EvaluateObservable::loadData: Could load outsourced data set: ConfIF=%d, outsourcedFile=%s\n",id,outFileName);
          exit(0);
        }

        outFile.close();

	delete[] outFileName;
      }
    }
    fclose(file);
  }
  deleteFileNameList(fileNames, fileCount);

  count = 0;
  for (int I=1; I<interimDataMAX; I++) {
    if (interimData[I] != NULL) count++;  
  }
  int totalCount = count;

  IDs = new int[count];
  data = new double*[count];
  count = 0;
  int dummyCount = 0;
  for (int I=1; I<interimDataMAX; I++) {
    if (interimData[I] != NULL) {
      if ((dummyCount>=relativeEvaluationStartIndex*totalCount) && (dummyCount<=relativeEvaluationEndIndex*totalCount)) {
        IDs[count] = I;
        data[count] = interimData[I];
        count++;
      } else {
        delete[] interimData[I];
	interimData[I] = NULL;
      }
      dummyCount++;
    }
  }
  
  delete[] interimData;
}


void EvaluateObservable::deleteLoadedData(int& count, int* &IDs, double** &data) {
  for (int I=0; I<count; I++) {
    delete[] data[I];
  }
  delete[] IDs;
  delete[] data;
  IDs = NULL;
  data = NULL;
  count = 0;
}


void EvaluateObservable::addXmlOutput(char* tag, char* des, double val, double err) {
  xmlOutput_Tag[xmlOutput_Count] = cloneString(tag);
  xmlOutput_Description[xmlOutput_Count] = cloneString(des);
  xmlOutput_Value[xmlOutput_Count] = val;
  xmlOutput_ValueError[xmlOutput_Count] = err;
  xmlOutput_Count++;
}
 
 
void EvaluateObservable::writeXmlOutput() {
  if (LogLevel>2) printf("Writing XML-File to disk for observable %s...",obsName);
  char* fileName = ioControl->getObservableXmlFileName(obsName);    
  FILE* file = fopen(fileName,"w");
  delete[] fileName;
  fprintf(file,"<?xml version=\'1.0\'?>\n");
  fprintf(file,"<%s>\n",obsName);
  
  for (int I=0; I<xmlOutput_Count; I++) {
    fprintf(file,"  <%s>\n",xmlOutput_Tag[I]);
    if ((!isInteger(xmlOutput_Value[I])) || (isNaN(xmlOutput_ValueError[I])) || (xmlOutput_ValueError[I]!=0)) {
      fprintf(file,"    <value>%1.15e</value>\n",xmlOutput_Value[I]);
      if ((isNaN(xmlOutput_ValueError[I])) || (xmlOutput_ValueError[I]!=0)) {
        fprintf(file,"    <error>%1.15e</error>\n",xmlOutput_ValueError[I]);    
      }
    } else {
      fprintf(file,"    <value>%d</value>\n",roundToInt(xmlOutput_Value[I]));    
    }
    fprintf(file,"    <description>%s</description>\n",xmlOutput_Description[I]);
    fprintf(file,"  </%s>\n",xmlOutput_Tag[I]);    
  }
  fprintf(file,"</%s>\n",obsName);
  fclose(file);
  if (LogLevel>2) printf("ready!\n");
}

  
char* EvaluateObservable::getObsName() {
  return obsName;
}


char* EvaluateObservable::getNickName() {
  return nickName;
}




void EvaluateObservable::addDependOnObsByName(char* obsName) {
  bool found = false;
  for (int I=0; I<otherObscount; I++) {
    if (otherObs[I]->isObs(obsName)) {
      dependOnObs[dependOnObsCount] = otherObs[I];
      dependOnObsCount++;      
      found = true;
      break;
    }
  }
  if (!found) {
    printf("EvaluateObservable::addDependOnObsByName: Cannot add observable %s\n", obsName);
    exit(0);
  }
}


EvaluateObservable* EvaluateObservable::getDependOnObsByName(char* name) {
  for (int I=0; I<dependOnObsCount; I++) {
    if (dependOnObs[I]->isObs(name)) {
      return dependOnObs[I];
    }
  }
  printf("EvaluateObservable::getDependOnObsByName: Cannot find observable %s\n", name);
  exit(0);
}


double* EvaluateObservable::getDataForID(int ID) {
  double* res = NULL;
  
  for (int I=0; I<dataAvailCount; I++) {
    if (dataAvailID[I]==ID) {
      res = new double[getAnalyzerResultsCount()];
      for (int I2=0; I2<getAnalyzerResultsCount(); I2++) {
        res[I2] = dataAvail[I][I2];
      }    
      return res;
    }  
  }
  printf("ERROR: EvaluateObservable::getDataForID could not deliver results for ID=%d\n",ID);
  exit(0);
}


void EvaluateObservable::addOtherObsPointer(EvaluateObservable* obs) {
  otherObs[otherObscount] = obs;
  otherObscount++;  
}


bool EvaluateObservable::reduceDataDueToDependencies() {
  bool dataReduced = false;

  if (dependOnObsCount>0) {
    if (LogLevel>1) printf("Reducing data of observable %s ...\n", obsName);
    int totalAnalyzerResultsCount = getAnalyzerResultsCount();
    for (int I=0; I<dependOnObsCount; I++) {
      if (LogLevel>1) printf("  depending on observable %s ==>", dependOnObs[I]->getObsName());
      int* depIDs = NULL;
      int depCount = 0;
      dependOnObs[I]->getDataIDs(depCount, depIDs);
      totalAnalyzerResultsCount += dependOnObs[I]->getAnalyzerResultsCount();
      int mismatchCount = 0;
    
      int p = 0;
      for (int I2=0; I2<dataAvailCount; I2++) {
        bool match = false;
        for (int I3=p; I3<depCount; I3++) {
          if (depIDs[I3]>dataAvailID[I2]) break;	
          if (depIDs[I3] == dataAvailID[I2]) {
  	    match = true;
  	    p = I3+1;
 	    break;
  	  }
        }
        if (!match) {
          //delete entry;  
  	  delete[] dataAvail[I2];
	  dataAvail[I2] = NULL;
	  dataReduced = true;
 	  mismatchCount++;
        }
      }
      if (LogLevel>1) printf("  %d mismatches\n", mismatchCount);      
      delete[] depIDs;
    }

    int dataRemoved = 0;  
    for (int I=0; I<dataAvailCount; I++) {
      if (dataAvail[I] == NULL) {
        for (int I2=I; I2<dataAvailCount-1; I2++) {
          dataAvail[I2] = dataAvail[I2+1];
          dataAvailID[I2] = dataAvailID[I2+1];
          dataAvailWeightAndSign[I2] = dataAvailWeightAndSign[I2+1];
        }
        I--;
        dataAvailCount--;
        dataRemoved++;
      }
    }
  
    if (LogLevel>1) printf("--> removed %d data sets in total\n", dataRemoved);
    
    //extend space for data
    int anaResCount = getAnalyzerResultsCount();    
    for (int I2=0; I2<dataAvailCount; I2++) {
      double* interim = new double[totalAnalyzerResultsCount];
      for (int I3=0; I3<anaResCount; I3++) {
        interim[I3] = dataAvail[I2][I3];      
      }
      delete[] dataAvail[I2];
      dataAvail[I2] = interim;
    }
    
    //fill new space with depending oberservable data
    int fillInPos = getAnalyzerResultsCount();
    for (int I=0; I<dependOnObsCount; I++) {
      int* depIDs;
      int depCount = 0;
      int depResCount = 0;
      double** depData = NULL;
    
      dependOnObs[I]->getData(depCount, depIDs, depResCount, depData);
    
      int p = 0;
      for (int I2=0; I2<dataAvailCount; I2++) {
        bool match = false;
        for (int I3=p; I3<depCount; I3++) {
          if (depIDs[I3]>dataAvailID[I2]) break;	
          if (depIDs[I3] == dataAvailID[I2]) {
  	    match = true;
	    for (int I4=0; I4<depResCount; I4++) {
	      dataAvail[I2][I4+fillInPos] = depData[I3][I4];
	    }
  	    p = I3+1;
 	    break;
  	  }
        }
        if (!match) {
	  printf("ERROR in EvaluateObservable::reduceDataDueToDependencies: No data match found!!!\n");
	  exit(0);
        }
      }
      delete[] depIDs;
      for (int I2=0; I2<depCount; I2++) {
        delete[] depData[I2];      
      }
      delete[] depData;
      fillInPos += depResCount;
    }
    if (LogLevel>1) printf("--> dataSets available %d with data-chunk size %d\n", dataAvailCount, totalAnalyzerResultsCount);
  }

  return dataReduced;
}


void EvaluateObservable::startLatexOutputSummaryTable() {
  LAPsystem->addDirectText("\\begin{center}\n");
  LAPsystem->addDirectText("\\begin{tabular}{lcr}\n");
}


void EvaluateObservable::addXMLonly(char* xmltag, char* des, double val, double error) {
  addXmlOutput(xmltag, des, val, error);
}


void EvaluateObservable::addXML_And_LatexOutputSummaryTableLine(char* xmltag, char* des, char* shortCut, double val, double error, char* unit, char* valFormat) {
  char* line = new char[2000];
  char* formatStr = new char[1000];
  
  if ((error!=0) && (unit!=NULL)) {
    snprintf(formatStr,1000,"%s & %s & ($%s \\pm %s$)\\,%s \\\\ \n",des,shortCut,valFormat,valFormat, unit);
    snprintf(line,2000,formatStr,val,error);
  }
  if ((error==0) && (unit!=NULL)) {
    snprintf(formatStr,1000,"%s & %s & %s \\,%s \\\\ \n",des,shortCut,valFormat, unit);
    snprintf(line,2000,formatStr,val,error);
  }
  if ((error!=0) && (unit==NULL)) {
    snprintf(formatStr,1000,"%s & %s & $%s \\pm %s$ \\\\ \n",des,shortCut,valFormat,valFormat);
    snprintf(line,2000,formatStr,val,error);
  }
  if ((error==0) && (unit==NULL)) {
    snprintf(formatStr,1000,"%s & %s & %s \\\\ \n",des,shortCut,valFormat);
    snprintf(line,2000,formatStr,val,error);
  }
  
  LAPsystem->addDirectText(line);
  
  addXmlOutput(xmltag, des, val, error);
  
  delete[] formatStr;
  delete[] line;
}


void EvaluateObservable::endLatexOutputSummaryTable() {
  LAPsystem->addDirectText("\\end{tabular}\n");
  LAPsystem->addDirectText("\\end{center}\n");
}
