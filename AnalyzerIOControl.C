AnalyzerIOControl::AnalyzerIOControl(char* extension, int jID) {
  if (LogLevel>1) printf("AnalyzerIOControl initialized with extension = %s, JobID = %d\n",extension,jID);
  JobID = jID;
  fileNameExtension = new char[1000];
  snprintf(fileNameExtension,1000,"%s",extension);
}


AnalyzerIOControl::~AnalyzerIOControl() { 
  delete[] fileNameExtension;
}


void AnalyzerIOControl::createWorkFolder() {
  char* baseDir = new char[1000];
  snprintf(baseDir,1000,"%s/data/results/pHMC/analysis",DataBaseDirectory);
  char* command = new char[1000];
  snprintf(command,1000,"cd %s\n mkdir subFolder%s",baseDir, fileNameExtension);

  if (LogLevel>2) printf("Creating Work-Folder with command:\n %s\n",command);
  system(command);

  delete[] baseDir;
  delete[] command;  
}


void AnalyzerIOControl::createObsFolder(char* obsName) {
  char* baseDir = getWorkFolderName();
  char* command = new char[1000];
  snprintf(command,1000,"cd %s\n mkdir %s\n cd %s\n mkdir outsourced",baseDir, obsName, obsName);

  if (LogLevel>2) printf("Creating Obs-Folder and outsourced-Folder with command:\n %s\n",command);
  system(command);

  delete[] baseDir;
  delete[] command;  
}


void AnalyzerIOControl::createAndPrepareObsPlotsFolder(char* obsName) {
  char* baseDir = getWorkFolderName();
  char* command = new char[1000];
  snprintf(command,1000,"cd %s\n mkdir %s\n cd %s\n mkdir plots",baseDir, obsName, obsName);

  if (LogLevel>2) printf("Creating Obs-Folder and plots folder with command:\n %s\n",command);
  system(command);
  
  snprintf(command,1000,"cp makePlot %s/%s/plots\n cp fixbb %s/%s/plots",baseDir, obsName, baseDir, obsName);

  if (LogLevel>2) printf("Copying makePlot and fixbb to plots folder with command:\n %s\n",command);
  system(command);

  delete[] baseDir;
  delete[] command;  
}


void AnalyzerIOControl::markAsInProgress(int ID) {
  if (LogLevel>2) printf("Mark confID = %d as in progress.\n", ID);
  char* fileName = new char[2000];
  char* workDir = getWorkFolderName();
  long int timeStamp = (long int) zeitwert();
  snprintf(fileName, 2000, "%s/inProgress%s_ConfID%d_JobID%d_Time%ld.info",workDir,fileNameExtension,ID,JobID,timeStamp);
  delete[] workDir;
  
  FILE* file = fopen(fileName,"w");  
  fprintf(file,"Conf %d in Progress on Job %d from time %f\n",ID,JobID,zeitwert());
  fclose(file);
    
  delete[] fileName;
}


void AnalyzerIOControl::unmarkAsInProgress(int ID) {
  if (LogLevel>2) printf("Unmark confID = %d.\n", ID);
  char* searchDir = getWorkFolderName();
  char** fileNames; 
  int fileCount = 0;
  char* searchString = new char[2000];
  snprintf(searchString, 2000, "inProgress%s_ConfID%d_JobID%d",fileNameExtension,ID,JobID);  
  getFileNameList(searchDir, searchString, fileNames, fileCount);
  
  if (fileCount == 0) {
    printf("AnalyzerIOControl::unmarkAsInProgress: InProgress-File for ID=%d and jobID=%d not found\n",ID,JobID);
  }
  if (fileCount > 1) {
    printf("ERROR: AnalyzerIOControl::unmarkAsInProgress: Several InProgress-Files for ID=%d and jobID=%d found\n",ID,JobID);
    for (int I=0; I<fileCount; I++) {
      printf("%s\n",fileNames[I]);
    }
    exit(0);
  }
  if (fileCount == 1) {
    char* command = new char[2000];
    snprintf(command,2000,"rm %s",fileNames[0]);    
    system(command);
    delete[] command;
  }  
  
  deleteFileNameList(fileNames, fileCount);  
  delete[] searchDir;
  delete[] searchString;  
}


bool AnalyzerIOControl::isInProgress(int ID) {
  char* searchDir = getWorkFolderName();
  char** fileNames; 
  int fileCount = 0;
  char* searchString = new char[2000];
  snprintf(searchString, 2000, "inProgress%s_ConfID%d",fileNameExtension,ID);  
  getFileNameList(searchDir, searchString, fileNames, fileCount);
  
  bool inProg = (fileCount > 0);
  
  deleteFileNameList(fileNames, fileCount);  
  delete[] searchDir;
  delete[] searchString;  

  return inProg;
}


void AnalyzerIOControl::getConfInProgressIDs(int& count, int* &confInProgIDs) {
  count = 0;
  confInProgIDs = NULL;  
  char* searchDir = getWorkFolderName();
  char** fileNames; 
  int fileCount = 0;
  char* searchString = new char[2000];
  snprintf(searchString, 2000, "inProgress%s",fileNameExtension);  
  getFileNameList(searchDir, searchString, fileNames, fileCount);
  count = fileCount;
  confInProgIDs = new int[fileCount];
  
  for (int I=0; I<fileCount; I++) {
    char* formatString = new char[2000];
    snprintf(formatString, 2000, "%s/inProgress%s_ConfID%%d_JobID%%d_Time%%ld.info",searchDir,fileNameExtension);
    int cID = -1;
    int jID = -1;
    long int timeStamp = -1;     
    if (sscanf(fileNames[I],formatString,&cID,&jID,&timeStamp) == 3) {
      confInProgIDs[I] = cID;
    } else {
      printf("ERROR: AnalyzerIOControl::getConfInProgressIDs: Could not correctly decode filename %s\n",fileNames[I]);
      exit(0);
    }
    delete[] formatString;
  }
  
  deleteFileNameList(fileNames, fileCount);  
  delete[] searchDir;
  delete[] searchString;    
}


void AnalyzerIOControl::removeDeprecatedInProgressFiles(double maxAgeinHours) {
  if (LogLevel>2) printf("Removing deprecated InProgress-Files older than %1.2f hours...\n",maxAgeinHours);
  char* searchDir = getWorkFolderName();
  char** fileNames; 
  int fileCount = 0;
  char* searchString = new char[2000];
  snprintf(searchString, 2000, "inProgress%s",fileNameExtension);  
  getFileNameList(searchDir, searchString, fileNames, fileCount);
  long int currentTimeStamp = (long int) zeitwert();
  
  for (int I=0; I<fileCount; I++) {
    char* formatString = new char[2000];
    snprintf(formatString, 2000, "%s/inProgress%s_ConfID%%d_JobID%%d_Time%%ld.info",searchDir,fileNameExtension);
    int cID = -1;
    int jID = -1;
    long int timeStamp = -1;         
    if (sscanf(fileNames[I],formatString,&cID,&jID,&timeStamp) == 3) {
      long int diff = currentTimeStamp - timeStamp;
      if (diff<0) {
        printf("ERROR: AnalyzerIOControl::removeDeprecatedInProgressFiles: TimeStamp newer than current time for inProgress-FILE %s\n", fileNames[I]);
	exit(0);
      }
      if (diff>3600*maxAgeinHours) {
        char* command = new char[2000];
        snprintf(command,2000,"rm %s",fileNames[I]);    
        system(command);
	if (LogLevel>2) printf("AnalyzerIOControl::removeDeprecatedInProgressFiles: %s deleted with command %s\n",fileNames[I], command);
        delete[] command;
      }
    }
    delete[] formatString;
  }
  
  deleteFileNameList(fileNames, fileCount);  
  delete[] searchDir;
  delete[] searchString;  
}


char* AnalyzerIOControl::getWorkFolderName() {
  char* fileName = new char[1000];
  snprintf(fileName,1000,"%s/data/results/pHMC/analysis/subFolder%s",DataBaseDirectory,fileNameExtension);
  return fileName;
}


char* AnalyzerIOControl::getObservableFolderName(char* obsName) {
  char* workDir = getWorkFolderName();
  char* folderName = new char[1000];
  snprintf(folderName,1000,"%s/%s",workDir,obsName);
  delete[] workDir;
  return folderName;
}

char* AnalyzerIOControl::getObservableOutsourcedFolderName(char* obsName) {
  char* workDir = getWorkFolderName();
  char* folderName = new char[1000];
  snprintf(folderName,1000,"%s/%s/outsourced",workDir,obsName);
  delete[] workDir;
  return folderName;
}


char* AnalyzerIOControl::getObservablePlotsFolderName(char* obsName) {
  char* workDir = getWorkFolderName();
  char* folderName = new char[1000];
  snprintf(folderName,1000,"%s/%s/plots",workDir,obsName);
  delete[] workDir;
  return folderName;
}


char* AnalyzerIOControl::getObservableFileName(char* obsName) {
  char* workDir = getWorkFolderName();
  char* fileName = new char[1000];
  snprintf(fileName,1000,"%s/%s/%s%s_JID%d.dat",workDir,obsName,obsName,fileNameExtension,JobID);
  delete[] workDir;
  return fileName;
}


char* AnalyzerIOControl::getObservableOutsourceFileName(char* obsName, int confID) {
  return getObservableOutsourceFileName(obsName, JobID, confID); 
}


char* AnalyzerIOControl::getObservableOutsourceFileName(char* obsName, int jID, int confID) {
  char* outDir = getObservableOutsourcedFolderName(obsName);
  char* fileName = new char[1000]; 
  snprintf(fileName,1000,"%s/outsourced_%s%s_JID%d_ConfID%d.dat",outDir,obsName,fileNameExtension,jID, confID);
  delete[] outDir;
  return fileName;
}


char* AnalyzerIOControl::getObservableXmlFileName(char* obsName) {
  char* obsDir = getObservableFolderName(obsName);
  char* fileName = new char[1000]; 
  snprintf(fileName,1000,"%s/%s%s_Results.xml",obsDir,obsName,fileNameExtension);
  delete[] obsDir;
  return fileName;
}


char* AnalyzerIOControl::getMainXmlFileName() {
  char* workDir = getWorkFolderName();
  char* fileName = new char[1000];
  snprintf(fileName,1000,"%s/Summary%s_Results.xml",workDir,fileNameExtension);
  delete[] workDir;
  return fileName;
}


char* AnalyzerIOControl::getObservableLatexBaseName(char* obsName) {
  char* latexBaseName = new char[1000];
  snprintf(latexBaseName,1000,"%s_%s",obsName,fileNameExtension);
  return latexBaseName;
}


char* AnalyzerIOControl::getObservableLatexBodyLocalFileName(char* obsName) {
  char* latexfileName = new char[1000];
  char* latexBaseName = getObservableLatexBaseName(obsName);
  snprintf(latexfileName,1000,"%s_body.tex",latexBaseName);
  delete[] latexBaseName;
  return latexfileName;
}


char* AnalyzerIOControl::getObservableLatexBodyFileName(char* obsName) {
  char* locallatexfileName = getObservableLatexBodyLocalFileName(obsName);
  char* latexfileName = new char[1000];
  char* obsDir = getObservableFolderName(obsName);
  snprintf(latexfileName,1000,"%s/%s",obsDir,locallatexfileName);
  delete[] obsDir;
  delete[] locallatexfileName;
  return latexfileName;
}


char* AnalyzerIOControl::getLatexSummaryBaseName() {
  char* latexBaseName = new char[1000];
  snprintf(latexBaseName,1000,"LatexSummary_%s",fileNameExtension);
  return latexBaseName;
}


char* AnalyzerIOControl::getStateDescriptorFileName() {
  char* fileName = new char[1000];
  snprintf(fileName,1000,"%s/data/results/pHMC/states/StateDescriptor%s.dat",DataBaseDirectory,fileNameExtension);
  return fileName;
}


char* AnalyzerIOControl::getPhiConfFileName(int ID) {
  char* fileName = new char[1000];
  snprintf(fileName,1000,"%s/data/results/pHMC/configurations/subFolder%s/PhiConf%s_%d.dat",DataBaseDirectory,fileNameExtension,fileNameExtension,ID);
  return fileName;
}


int AnalyzerIOControl::getTotalPhiConfCount() {
  int count = 0;
  //Bug-Fix for lost configurations
  char* fileName = getPhiConfFileName(10);
  FILE* file = fopen(fileName,"r");
  if (file != NULL) {
    count = 10;
    fclose(file);
  }
  delete[] fileName;  
  
  
  while (true) {
    char* fileName = getPhiConfFileName(count+1);
    FILE* file = fopen(fileName,"r");
    delete[] fileName;  
    
    if (file == NULL) {
      char* fileName2 = getPhiConfFileName(count+2);
      FILE* file2 = fopen(fileName2,"r");
      delete[] fileName2;  
      
      if (file2 == NULL) return count;
      
      fclose(file2);
    } else {    
      fclose(file);
    }
    
    count++;    
  }
  return count;
}


void AnalyzerIOControl::getObservableFileNameList(char* obsName, char** &fileNames, int& fileCount) {
  char* workDir = getWorkFolderName(); 
  char* searchDir = new char[2000];
  snprintf(searchDir,2000,"%s/%s",workDir,obsName);
  char* searchStr = new char[1000];
  snprintf(searchStr,1000,"%s%s",obsName,fileNameExtension);  
  getFileNameList(searchDir, searchStr, fileNames, fileCount);
  delete[] searchDir;
  delete[] workDir;
  delete[] searchStr;
}


int AnalyzerIOControl::getJobID() {
  return JobID;
}
