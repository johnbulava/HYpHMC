#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include <stddef.h>
#include <sys/types.h>
#include <dirent.h>

#include "Complex.h"
#include "Quat.h"
#include "Global.h"
#include "BootStrapClass.h"
#include "ControlLogger.h"
#include "AutoCorrelation.h"


#define Critical_Nu 0.5
#define Critical_Gamma 1.2  
#define AutoCorrAnalysisLength 4000

class ParameterSelectorType {
public:
  double value[100];
  int valueCount;
  char* readInstruction;
  char* varName;
  
  ParameterSelectorType(char* pre) {
    valueCount = 0;
    readInstruction = new char[100];
    varName = new char[100];
    snprintf(readInstruction,100,"%s=%%s",pre); 
    snprintf(varName,100,"%s",pre);        
  }
  
   ~ParameterSelectorType() {
    delete[] readInstruction;
    delete[] varName;
  }
  
  void readSpecification(char* spec) {
    double v;
    char* s = new char[1000];
    if (sscanf(spec,readInstruction,s)==1) {
      while (strlen(s)>0) {
        if (sscanf(s,"%lf",&v)==1) {
	  value[valueCount] = v;
          printf("Specification read: %s = %f\n",varName,value[valueCount]);
	  valueCount++;
	  
	  int p = -1;
	  int I;
	  for (I=0; I<(int)strlen(s); I++) {
	    if (s[I] == ',') {
	      p = I;
	      break;
	    }
	  }
	  if (p==-1) break;
	  for (I=0; I<(int)strlen(s)-p; I++) {
            s[I] = s[I+p+1];
	  }
	} else {
	  break;
	}
      }
    }
    delete[] s;
  }
  
  bool isSelected(double p) {
    int I;
    if (valueCount == 0) return true;
    for (I=0; I<valueCount; I++) {
      if (!(value[I]==value[I])) {
        if (!(p==p)) return true;
      } else {
        if (value[I] == p) return true;
      }
    }
    return false;    
  }
};
  

//Variablen
ParameterSelectorType Parameter_yN_Selector("yN");
ParameterSelectorType Parameter_rho_Selector("rho");
ParameterSelectorType Parameter_r_Selector("r");
ParameterSelectorType Parameter_kappa_Selector("kappa");
ParameterSelectorType Parameter_lambda_Selector("lambda");
ParameterSelectorType Parameter_L0_Selector("L0");
ParameterSelectorType Parameter_L1_Selector("L1");
ParameterSelectorType Parameter_L2_Selector("L2");
ParameterSelectorType Parameter_L3_Selector("L3");
ParameterSelectorType Parameter_Nf_Selector("Nf");


  
class FileListEntryType {
public:
  char* fileName;
  double kappa, lambda, Y, rho, r;
  int L0,L1,L2,L3, Nf;
  FileListEntryType* next;
  bool readout;
  
  FileListEntryType() {
    fileName = new char[500];
    kappa = 0;
    lambda = 0;
    Y = 0;
    rho = 0;
    r = 0;
    L0 = 0;
    L1 = 0;
    L2 = 0;
    L3 = 0;
    Nf = 0;
    next = NULL;
    readout = false;
  }  
  
  void switchContents(FileListEntryType* x) {
    double d1;
    int d2;
    char* d3;
    bool d4;
    
    d3 = fileName;
    fileName = x->fileName;
    x->fileName = d3;    
    
    d1 = kappa;
    kappa = x->kappa;
    x->kappa = d1;

    d1 = lambda;
    lambda = x->lambda;
    x->lambda = d1;
    
    d1 = Y;
    Y = x->Y;
    x->Y = d1;

    d1 = rho;
    rho = x->rho;
    x->rho = d1;

    d1 = r;
    r = x->r;
    x->r = d1;
    
    d2 = L0;
    L0 = x->L0;
    x->L0 = d2;
    
    d2 = L1;
    L1 = x->L1;
    x->L1 = d2;
    
    d2 = L2;
    L2 = x->L2;
    x->L2 = d2;
    
    d2 = L3;
    L3 = x->L3;
    x->L3 = d2;
    
    d2 = Nf;
    Nf = x->Nf;
    x->Nf = d2;

    d4 = readout;
    readout = x->readout;
    x->readout = d4;
  }
  
  void print() {
    printf("Kappa: %1.3f, Lambda: %1.3f, Y: %1.3f, L: %dx%dx%dx%d, Nf %d, Rho: %1.3f, R: %1.3f, filename = %s\n", kappa,lambda, Y, L0, L1, L2, L3, Nf, rho, r, fileName);
  }
};


class ResultDescriptorType {
public:
  double kappa, lambda, Y, rho, r;
  int L0, L1, L2, L3, Nf;
  double V, Lavg;
  double averagePhiNorm;
  double averagePhiSqrNorm;
  double averagePhiQuartNorm;
  double sigmaPhiNorm;
  double sigmaPhiSqrNorm;
  double sigmaPhiQuartNorm;
  double averageStaggeredPhiNorm;
  double averageStaggeredPhiSqrNorm;
  double averageStaggeredPhiQuartNorm;
  double sigmaStaggeredPhiNorm;
  double sigmaStaggeredPhiSqrNorm;
  double sigmaStaggeredPhiQuartNorm;
  double Susceptibility;
  double sigmaSusceptibility;
  double StaggeredSusceptibility;
  double sigmaStaggeredSusceptibility;
  double BinderCumulant;
  double sigmaBinderCumulant;
  double StaggeredBinderCumulant;
  double sigmaStaggeredBinderCumulant;  
  double DSL;

  int count;
  bool readout;
  
  ResultDescriptorType() {
    count = 0;
    kappa = 0;
    lambda = 0;
    Y = 0;
    rho = 0;
    r = 0;
    L0 = 0;
    L1 = 0;
    L2 = 0;
    L3 = 0;
    V = 0;
    Lavg = 0;
    Nf = 0;
    DSL = 0;
    averagePhiNorm = 0;
    averagePhiSqrNorm = 0;
    averagePhiQuartNorm = 0;
    sigmaPhiNorm = 0;
    sigmaPhiSqrNorm = 0;
    sigmaPhiQuartNorm = 0;
    averageStaggeredPhiNorm = 0;
    averageStaggeredPhiSqrNorm = 0;
    averageStaggeredPhiQuartNorm = 0;
    sigmaStaggeredPhiNorm = 0;
    sigmaStaggeredPhiSqrNorm = 0;
    sigmaStaggeredPhiQuartNorm = 0;
    Susceptibility = 0;
    sigmaSusceptibility = 0;
    StaggeredSusceptibility = 0;
    sigmaStaggeredSusceptibility = 0;
    BinderCumulant = 0;
    sigmaBinderCumulant = 0;
    StaggeredBinderCumulant = 0;
    sigmaStaggeredBinderCumulant = 0;  
    readout = false;
  }

  void print() {
    printf("Kap: %1.3f, Lam: %2.3f, Y: %1.3f, Rho: %1.1f, R: %1.1f, L: %2dx%2dx%2dx%2d, Nf: %2d, Stat: %5d, Phi: %1.4f+-%1.4f, St.Phi: %1.3f+-%1.3f, DSL: %1.0f, Cut: %1.1f GeV\n", kappa,lambda,Y,rho,r,L0,L1,L2,L3,Nf,count,averagePhiNorm,sigmaPhiNorm,averageStaggeredPhiNorm,sigmaStaggeredPhiNorm,DSL/log(10), Physical_VEV_GeV/(sqrt(2.0*kappa)*averagePhiNorm));
  }
  
};
  

//Variables
FileListEntryType fileList;
ResultDescriptorType* MCresults = NULL;
int MCresultCount = 0;
int CommandCount = 0;
char** commands;
bool FLAG_nozeros = false;
bool FLAG_nochi = false;
bool FLAG_wilson = false;
bool FLAG_bootStrap = false;
bool FLAG_ExtrapolateNf = false;
bool FLAG_hybrid = false;
bool FLAG_pHMC = false;
bool FLAG_ImpPhiDyn = false;
bool FLAG_tildeValues = false;
bool FLAG_TexLog = false;
bool FLAG_GnuplotFit = false;
bool FLAG_TransCut = false;
bool FLAG_MCTrajectoryLog = false;
bool FLAG_AutoCorrLog = false;
bool FLAG_FiniteSizeEffect = false;
bool FLAG_Test = false;
int removeThermal = 0;
int bootStrapIter = 0;
int extrapolateToNf = 1;
double transitionCutValue = 0.5;
ResultDescriptorType** FitDataResults = NULL;
int FitDataCount = 0;
ControlLogger FitControlLogger;
ControlLogger MCTrajectoryLogger;
ControlLogger AutoCorrLogger;


double calcScaleFacC(double lambda, double Nf) {
  double p = (1.0-2.0*lambda)/2.0;
  return sqrt(p + sqrt(p*p + 2.0*Nf*lambda));
}


double calcScaleFacCInverse(double lambda_N, double Nf) {
  double p = (1.0-2.0*lambda_N*Nf)/2.0;
  return sqrt(p + sqrt(p*p + 2.0*lambda_N));
}


void loadFileList() {
  DIR *dp = NULL;
  struct dirent *ep;
  FileListEntryType* pointer = &fileList;
  double dummy_kappa, dummy_lambda, dummy_Y, dummy_rho, dummy_r;
  int dummy_Nf, dummy_L0, dummy_L1, dummy_L2, dummy_L3;
     
     
  char* dataDirectory = new char[1000];   
  if ((!FLAG_hybrid) && (!FLAG_ImpPhiDyn) && (!FLAG_pHMC))snprintf(dataDirectory,1000,"%s%s",DataBaseDirectory,"/data/results/metro"); 
  if ((!FLAG_hybrid) && (FLAG_ImpPhiDyn) && (!FLAG_pHMC)) snprintf(dataDirectory,1000,"%s%s",DataBaseDirectory,"/data/results/metroImprovedPhiDynamics"); 
  if ((FLAG_hybrid)  && (!FLAG_ImpPhiDyn) && (!FLAG_pHMC))snprintf(dataDirectory,1000,"%s%s",DataBaseDirectory,"/data/results/hybrid"); 
  if ((FLAG_hybrid)  && (FLAG_ImpPhiDyn) && (!FLAG_pHMC)) snprintf(dataDirectory,1000,"%s%s",DataBaseDirectory,"/data/results/hybridImprovedPhiDynamics"); 
  if ((!FLAG_hybrid)  && (!FLAG_ImpPhiDyn) && (FLAG_pHMC))snprintf(dataDirectory,1000,"%s%s",DataBaseDirectory,"/data/results/pHMC/measure"); 
  
  dp = opendir (dataDirectory);
  
  if (dp != NULL) {
    while ((ep = readdir (dp))) {
      bool fileFound = false;
      if (FLAG_Test) {
        if (sscanf(ep->d_name,"testN%dNf%dKap%lfLam%lfY%lfRho%lfR%lf.dat",&(dummy_L0),&(dummy_Nf),&(dummy_kappa),&(dummy_lambda),&(dummy_Y),&(dummy_rho),&(dummy_r))==7) {
	  dummy_L1 = dummy_L0;
	  dummy_L2 = dummy_L0;
	  dummy_L3 = dummy_L0;
          fileFound = true;	
	}
        if (sscanf(ep->d_name,"testL%dx%dx%dx%dNf%dKap%lfLam%lfY%lfRho%lfR%lf.dat",&(dummy_L0),&(dummy_L1),&(dummy_L2),&(dummy_L3),&(dummy_Nf),&(dummy_kappa),&(dummy_lambda),&(dummy_Y),&(dummy_rho),&(dummy_r))==10) {
          fileFound = true;	
	}
	
      } else {
        if (sscanf(ep->d_name,"measureN%dNf%dKap%lfLam%lfY%lfRho%lfR%lf.dat",&(dummy_L0),&(dummy_Nf),&(dummy_kappa),&(dummy_lambda),&(dummy_Y),&(dummy_rho),&(dummy_r))==7) {
	  dummy_L1 = dummy_L0;
	  dummy_L2 = dummy_L0;
	  dummy_L3 = dummy_L0;
          fileFound = true;	
	}
        if (sscanf(ep->d_name,"measureL%dx%dx%dx%dNf%dKap%lfLam%lfY%lfRho%lfR%lf.dat",&(dummy_L0),&(dummy_L1),&(dummy_L2),&(dummy_L3),&(dummy_Nf),&(dummy_kappa),&(dummy_lambda),&(dummy_Y),&(dummy_rho),&(dummy_r))==10) {
          fileFound = true;	
	}
      }
      if (fileFound) {
        if (Parameter_yN_Selector.isSelected(dummy_Y) && Parameter_rho_Selector.isSelected(dummy_rho) && 
        Parameter_r_Selector.isSelected(dummy_r) && Parameter_kappa_Selector.isSelected(dummy_kappa) &&  
        Parameter_lambda_Selector.isSelected(dummy_lambda) &&  Parameter_L0_Selector.isSelected(dummy_L0) && 
	Parameter_L1_Selector.isSelected(dummy_L1) && Parameter_L2_Selector.isSelected(dummy_L2) && 
	Parameter_L3_Selector.isSelected(dummy_L3) && Parameter_Nf_Selector.isSelected(dummy_Nf)) {
      
          printf("Registering data file: '%s'\n",ep->d_name);	    
          FileListEntryType* entry = new FileListEntryType();
	  entry->kappa = dummy_kappa;
  	  entry->lambda = dummy_lambda;
 	  entry->Y = dummy_Y;
  	  entry->rho = dummy_rho;
	  entry->r = dummy_r;
	  entry->L0 = dummy_L0;
	  entry->L1 = dummy_L1;
	  entry->L2 = dummy_L2;
	  entry->L3 = dummy_L3;
	  entry->Nf = dummy_Nf;
	
          snprintf(entry->fileName,500,"%s/%s",dataDirectory,ep->d_name);

          pointer->next = entry;
          pointer = entry;
	}
      }
    }
    closedir (dp);
  } else {
    printf("Couldn't open the data-directory!!!\n");
    exit(0);
  }
  delete[] dataDirectory;
  
  while (true) {
    pointer = fileList.next;
    bool change = false;
    while (pointer!=NULL) {
      FileListEntryType* pointer2 = pointer->next;
      while (pointer2!=NULL) {
        if (pointer->Y>pointer2->Y) {
          pointer->switchContents(pointer2);
	  change = true;
	}   
        pointer2 = pointer2->next;
      }
      pointer = pointer->next;
    }
    if (!change) break;
  }
}


void printFileList() {
  FileListEntryType* pointer = fileList.next;
  while (pointer!=NULL) {
    pointer->print();
    pointer = pointer->next;
  }
}


ResultDescriptorType calculateMCresult(bool deliverRes, double kappa, double lambda, double Y, double rho, double r, int L0, int L1, int L2, int L3, int Nf) {
  ResultDescriptorType res;
  double dummy_phi, dummy_stagphi;
  double dummy_det_Log = 0.0;
  double dummy_logDet_1, dummy_logDet_2, dummy_logDet_3, dummy_logDet_4;
  double dummy_logDet_5, dummy_logDet_6, dummy_logDet_7, dummy_logDet_8;
  char* restLine = new char[1000];
  double* dataDet = new double[1000000];
  double* dataPhi = new double[1000000];
  double* dataStagPhi = new double[1000000];
  int dataCount = 0;
  double* Traj_X = new double[1000000];
  int RunCount = 0;
  int* RunLengths = new int[1000];
  int I;

  
  res.Y = Y;
  res.kappa = kappa;
  res.lambda = lambda;
  res.Nf = Nf;
  res.rho = rho;
  res.r = r;
  res.L0 = L0;
  res.L1 = L1;
  res.L2 = L2;
  res.L3 = L3;
  res.V = L0*L1*L2*L3;
  res.Lavg = sqrt(sqrt(res.V));
  res.averagePhiNorm = 0;
  res.averagePhiSqrNorm = 0;
  res.averagePhiQuartNorm = 0;
  res.averageStaggeredPhiNorm = 0;
  res.averageStaggeredPhiSqrNorm = 0;
  res.averageStaggeredPhiQuartNorm = 0;
  res.sigmaPhiNorm = 0;
  res.sigmaPhiSqrNorm = 0;
  res.sigmaPhiQuartNorm = 0;
  res.sigmaStaggeredPhiNorm = 0;
  res.sigmaStaggeredPhiSqrNorm = 0;
  res.sigmaStaggeredPhiQuartNorm = 0;
  res.Susceptibility = 0;
  res.sigmaSusceptibility = 0;
  res.StaggeredSusceptibility = 0;
  res.sigmaStaggeredSusceptibility = 0;
  res.BinderCumulant = 0;
  res.sigmaBinderCumulant = 0;
  res.StaggeredBinderCumulant = 0;
  res.sigmaStaggeredBinderCumulant = 0;  
  res.count = 0;
  res.DSL = 0;


  
  double CInverse = 1;
  double C = 1;
  double DetExp = 1;
  if (FLAG_ExtrapolateNf) {
    CInverse = calcScaleFacCInverse(lambda, Nf);

    res.Y = Y / CInverse;
    res.kappa = kappa / sqr(CInverse);
    res.lambda = lambda / sqr(sqr(CInverse));
    res.Nf = 1;

    C = calcScaleFacC(res.lambda, extrapolateToNf);
  
    res.Y = res.Y / C;
    res.kappa = res.kappa / sqr(C);
    res.lambda = res.lambda / sqr(sqr(C));
    res.Nf = extrapolateToNf;
    DetExp = extrapolateToNf;
    DetExp = DetExp / Nf;
  }
  

  //Bring this to the tilde-values
  if (FLAG_tildeValues) {
    res.Y = res.Y * sqrt(res.Nf);
    res.kappa = res.kappa;
    res.lambda = res.lambda * res.Nf;
  }
  
  FileListEntryType* pointer = fileList.next;
  while (pointer!=NULL) {
    if ((pointer->kappa == kappa) && ((pointer->lambda == lambda) || (((!(lambda==lambda))) && (!(pointer->lambda==pointer->lambda)))) && (pointer->Y == Y) && (pointer->rho == rho) && (pointer->r == r) && (pointer->L0 == L0) && (pointer->L1 == L1) && (pointer->L2 == L2) && (pointer->L3 == L3) &&(pointer->Nf == Nf)){
      if (deliverRes) {
        //Evaluate this file
        FILE* file;
        printf("Reading file %s\n",pointer->fileName);
        file = fopen(pointer->fileName,"r");
        RunCount++;
	RunLengths[RunCount-1] = 0;
	int ExtraThermCount = 0;
	
        while (fscanf(file,"%lf %lf %lf	%lf %lf	%lf %lf	%lf %lf	%lf", &dummy_phi, &dummy_stagphi, &dummy_logDet_1, &dummy_logDet_2, &dummy_logDet_3, &dummy_logDet_4, &dummy_logDet_5, &dummy_logDet_6, &dummy_logDet_7,&dummy_logDet_8)==10) {
          fgets(restLine, 1000, file);

          //Choose correct determinant
          if ((FLAG_nozeros==false) && (FLAG_nochi==false) && (FLAG_wilson==false)) dummy_det_Log = dummy_logDet_1;
          if ((FLAG_nozeros==true)  && (FLAG_nochi==false) && (FLAG_wilson==false)) dummy_det_Log = dummy_logDet_2;
          if ((FLAG_nozeros==false) && (FLAG_nochi==true)  && (FLAG_wilson==true)) dummy_det_Log = dummy_logDet_3;
          if ((FLAG_nozeros==true)  && (FLAG_nochi==true)  && (FLAG_wilson==true)) dummy_det_Log = dummy_logDet_4;

          if ((FLAG_nozeros==false) && (FLAG_nochi==true) && (FLAG_wilson==false)) dummy_det_Log = dummy_logDet_5;
          if ((FLAG_nozeros==true)  && (FLAG_nochi==true) && (FLAG_wilson==false)) dummy_det_Log = dummy_logDet_6;
          if ((FLAG_nozeros==false) && (FLAG_nochi==false) && (FLAG_wilson==true)) dummy_det_Log = dummy_logDet_7;
          if ((FLAG_nozeros==true)  && (FLAG_nochi==false) && (FLAG_wilson==true)) dummy_det_Log = dummy_logDet_8;

	  //NaN - Test
	  if (dummy_det_Log == dummy_det_Log) {
	    ExtraThermCount++;
	    if (ExtraThermCount>removeThermal) {
  	      //Rescale Phi - fields
	      dummy_phi *= C*CInverse;// / sqrt(res.Nf);
	      dummy_stagphi *= C*CInverse;// / sqrt(res.Nf);

              //Exponentiate Determinant Adequately
              dataDet[dataCount] = dummy_det_Log * DetExp;
	      dataPhi[dataCount] = dummy_phi;
	      dataStagPhi[dataCount] = dummy_stagphi;

	      Traj_X[dataCount] = RunLengths[RunCount-1];
              RunLengths[RunCount-1]++;
	      dataCount++;
	    }
	  }
	}
        fclose(file);	
      }
      pointer->readout = true;
    }
    pointer = pointer->next;
  }
  delete[] restLine;
  

  if ((dataCount>0) && (deliverRes)) {
    double avgDetLog = 0;
    for (I=0; I<dataCount; I++) avgDetLog += dataDet[I];
    avgDetLog /= dataCount;
    for (I=0; I<dataCount; I++) dataDet[I] -= avgDetLog;
    for (I=0; I<dataCount; I++) dataDet[I] = exp(dataDet[I]);
  
    AutoCorrelation auto1(5,AutoCorrAnalysisLength);
    AutoCorrelation auto2(5,AutoCorrAnalysisLength);
    auto1.loadData(RunCount, RunLengths, dataDet, dataPhi);
    auto2.loadData(RunCount, RunLengths, dataDet, dataStagPhi);

    BootStrapClass boot1;
    BootStrapClass boot2;
    if (FLAG_bootStrap) {
      for (I=0; I<dataCount; I++) {
        boot1.addToDataPool(dataPhi[I], dataDet[I]);
        boot2.addToDataPool(dataStagPhi[I], dataDet[I]);
      }
    }

    double detMin = 0;
    double detMax = 0;
    for (I=0; I<dataCount; I++) {
      dataPhi[I] *= dataDet[I];
      dataStagPhi[I] *= dataDet[I];
      if ((dataDet[I]>detMax) || (detMax==0)) detMax = dataDet[I];
      if ((dataDet[I]<detMin) || (detMin==0)) detMin = dataDet[I];
    }
  
    if (FLAG_MCTrajectoryLog) {
      char* section = new char[1000];
      snprintf(section,1000,"Trajectories for $y_N=%1.3f$, $\\lambda=%1.3f$, $\\kappa=%1.3f$, $L=%dx%dx%dx%d$,$N_f=%d$",Y, lambda, kappa, L0, L1, L2, L3, Nf);
      MCTrajectoryLogger.addSection(section);
      delete[] section;
    
      int index = 0;
      for (I=0; I<RunCount; I++) {
        char* caption = new char[1000];
        snprintf(caption,1000,"Monte Carlo trajectory for $\\Phi$ (left) and staggered $\\Phi$ (right) for $yN=%1.3f$, $\\lambda=%1.3f$, $\\kappa=%1.3f$, $L=%dx%dx%dx%d$, $N_f=%d$",Y,lambda,kappa,L0,L1,L2,L3,Nf);
        double* p1 = &(Traj_X[index]);
        double* p2 = &(dataPhi[index]);
        double* p3 = &(dataStagPhi[index]);
        MCTrajectoryLogger.addGnuplotTableLinePlot("Phi", "Staggered Phi", caption, "Monte Carlo Time", "Monte Carlo Trajectory", p1, p2, p3, RunLengths[I], "", "");
        delete[] caption;
        index += RunLengths[I];
      }
      MCTrajectoryLogger.clearPage();
      MCTrajectoryLogger.newPage();
    }    


    res.count = dataCount;
    res.DSL = (double)(logl(detMax/detMin));
    double VolumeFac = res.V;
    
    if (FLAG_bootStrap) {
      boot1.doBootStrapWithWeights(bootStrapIter, res.averagePhiNorm, res.sigmaPhiNorm, res.Susceptibility, res.sigmaSusceptibility, res.BinderCumulant, res.sigmaBinderCumulant);
      boot2.doBootStrapWithWeights(bootStrapIter, res.averageStaggeredPhiNorm, res.sigmaStaggeredPhiNorm, res.StaggeredSusceptibility, res.sigmaStaggeredSusceptibility, res.StaggeredBinderCumulant , res.sigmaStaggeredBinderCumulant);
      res.Susceptibility *= VolumeFac;
      res.sigmaSusceptibility *= VolumeFac;
      res.StaggeredSusceptibility *= VolumeFac;
      res.sigmaStaggeredSusceptibility *= VolumeFac;
    } else {
      char* fitCommand1 = new char[1000];
      char* fitCommand2 = new char[1000];
      char* caption = new char[1000];
      char* section = new char[1000];
      snprintf(section,1000,"Auto correlation analysis for $y_N=%1.3f$, $\\lambda=%1.3f$, $\\kappa=%1.3f$, $L=%dx%dx%dx%d$, $N_f=%d$",Y, lambda, kappa, L0,L1,L2,L3, Nf);
      AutoCorrLogger.addSection(section);
      delete[] section;

      res.averagePhiNorm = auto1.getAverage(1);
      res.averagePhiSqrNorm = auto1.getAverage(2);
      res.averagePhiQuartNorm = auto1.getAverage(4);
      res.averageStaggeredPhiNorm = auto2.getAverage(1);
      res.averageStaggeredPhiSqrNorm = auto2.getAverage(2);
      res.averageStaggeredPhiQuartNorm = auto2.getAverage(4);
    
      ComplexVector derivatives(5);
      derivatives.setZero();
      derivatives.vectorElements[1].x = 1;
      res.sigmaPhiNorm = auto1.estimateCombinedError(derivatives);
      res.sigmaStaggeredPhiNorm = auto2.estimateCombinedError(derivatives);
      snprintf(fitCommand1,1000,"replot %f notitle\n",auto1.getTotalN()*sqr(res.sigmaPhiNorm)/(auto1.getCombinedCFunction())[0]);
      snprintf(fitCommand2,1000,"replot %f notitle\n",auto2.getTotalN()*sqr(res.sigmaStaggeredPhiNorm)/(auto2.getCombinedCFunction())[0]);
      snprintf(caption,1000,"Auto correlation from $\\Gamma$-strategy for $\\Phi$ (left) and staggered $\\Phi$ (right) with total data count = %d and auto-correlation times: %1.2f (m) and %1.2f (s).", auto1.getTotalN(),auto1.estimateAutoCorrelationTime(),auto2.estimateAutoCorrelationTime());
      AutoCorrLogger.addTablePlot("$\\Phi$", "Staggered $\\Phi$", caption, "$W$", "$W$" , "Auto-correlation $C(W)/\\Gamma(0)$","Auto-correlation $C(W)/\\Gamma(0)$",
       auto1.getWstepArray() , auto2.getWstepArray(), auto1.getCombinedReducedCFunction(), auto2.getCombinedReducedCFunction(),
       auto1.getCombinedReducedCFunctionErrors(), auto2.getCombinedReducedCFunctionErrors(), 
       auto1.getGammaWMax(), auto2.getGammaWMax() , fitCommand1, fitCommand2);
      

      derivatives.setZero();
      derivatives.vectorElements[2].x = 1;
      res.sigmaPhiSqrNorm = auto1.estimateCombinedError(derivatives);
      res.sigmaStaggeredPhiSqrNorm = auto2.estimateCombinedError(derivatives);

      derivatives.setZero();
      derivatives.vectorElements[4].x = 1;
      res.sigmaPhiQuartNorm = auto1.estimateCombinedError(derivatives);
      res.sigmaStaggeredPhiQuartNorm = auto2.estimateCombinedError(derivatives);

      res.Susceptibility = VolumeFac * (res.averagePhiSqrNorm - sqr(res.averagePhiNorm));
      res.StaggeredSusceptibility = VolumeFac * (res.averageStaggeredPhiSqrNorm - sqr(res.averageStaggeredPhiNorm));

      derivatives.setZero();
      derivatives.vectorElements[2].x = VolumeFac;
      derivatives.vectorElements[1].x = -2*VolumeFac*auto1.getAverage(1);
      res.sigmaSusceptibility = auto1.estimateCombinedError(derivatives);
      derivatives.vectorElements[1].x = -2*VolumeFac*auto2.getAverage(1);
      res.sigmaStaggeredSusceptibility = auto2.estimateCombinedError(derivatives);
      snprintf(fitCommand1,1000,"replot %f notitle\n",auto1.getTotalN()*sqr(res.sigmaSusceptibility)/(auto1.getCombinedCFunction())[0]);
      snprintf(fitCommand2,1000,"replot %f notitle\n",auto2.getTotalN()*sqr(res.sigmaStaggeredSusceptibility)/(auto2.getCombinedCFunction())[0]);
      snprintf(caption,1000,"Auto correlation from $\\Gamma$-strategy for susceptibility $\\chi$ (left) and staggered $\\chi$ (right) with total data count = %d.", auto1.getTotalN());
      AutoCorrLogger.addTablePlot("$\\chi$", "Staggered $\\chi$", caption, "$W$", "$W$" , "Auto-correlation $C(W)/\\Gamma(0)$", "Auto-correlation $C(W)/\\Gamma(0)$",
       auto1.getWstepArray() , auto2.getWstepArray(), auto1.getCombinedReducedCFunction(), auto2.getCombinedReducedCFunction(),
       auto1.getCombinedCFunctionErrors(), auto2.getCombinedCFunctionErrors(), 
       auto1.getGammaWMax(), auto2.getGammaWMax() , fitCommand1, fitCommand2);

      res.BinderCumulant = 2 -  res.averagePhiQuartNorm / sqr(res.averagePhiSqrNorm);
      res.StaggeredBinderCumulant = 2 -  res.averageStaggeredPhiQuartNorm / sqr(res.averageStaggeredPhiSqrNorm);
     
      derivatives.setZero();
      derivatives.vectorElements[2].x = 2.0 * res.averagePhiQuartNorm / (sqr(res.averagePhiSqrNorm) * res.averagePhiSqrNorm);
      derivatives.vectorElements[4].x = -1.0 / sqr(res.averagePhiSqrNorm);
      res.sigmaBinderCumulant = auto1.estimateCombinedError(derivatives);
      derivatives.setZero();
      derivatives.vectorElements[2].x = 2.0 * res.averageStaggeredPhiQuartNorm / (sqr(res.averageStaggeredPhiSqrNorm) * res.averageStaggeredPhiSqrNorm);
      derivatives.vectorElements[4].x = -1.0 / sqr(res.averageStaggeredPhiSqrNorm);
      res.sigmaStaggeredBinderCumulant = auto2.estimateCombinedError(derivatives);
      snprintf(fitCommand1,1000,"replot %f notitle\n",auto1.getTotalN()*sqr(res.sigmaBinderCumulant));
      snprintf(fitCommand2,1000,"replot %f notitle\n",auto2.getTotalN()*sqr(res.sigmaStaggeredBinderCumulant));
      snprintf(caption,1000,"Auto correlation ($\\Gamma$-strategy) for binder cumulant $U$ (left) and staggered $U$ (right) with total data count = %d.", auto1.getTotalN());
//      AutoCorrLogger.addTablePlot("$\\chi$", "Staggered $\\chi$", caption, "$W$", "$W$" , "Auto correlation $C(W)/\\Gamma(0)$", "Auto correlation $C(W)/\\Gamma(0)$",
//       auto1.getWstepArray() , auto2.getWstepArray(), auto1.getCombinedCFunction(), auto2.getCombinedCFunction(),
//       auto1.getCombinedCFunctionErrors(), auto2.getCombinedCFunctionErrors(), 
//       auto1.getGammaWMax(), auto2.getGammaWMax() , fitCommand1, fitCommand2);
      
     
      AutoCorrLogger.clearPage();
      AutoCorrLogger.newPage();
      delete[] fitCommand1;
      delete[] fitCommand2;
      delete[] caption;
    }
  }

  delete[] Traj_X;
  delete[] dataDet;
  delete[] dataPhi;
  delete[] dataStagPhi;
  delete[] RunLengths;
  return res;
}


void readMCresults() {
  //Reset readout-flag
  FileListEntryType* pointer = fileList.next;
  while (pointer!=NULL) {
    pointer->readout = false;
    pointer = pointer->next;
  }

  MCresultCount = 0;
  pointer = fileList.next;
  while (pointer!=NULL) {
    if (!pointer->readout) {
      calculateMCresult(false,pointer->kappa,pointer->lambda,pointer->Y,pointer->rho,pointer->r,pointer->L0,pointer->L1,pointer->L2,pointer->L3,pointer->Nf);
      MCresultCount++;
    }
    pointer = pointer->next;
  }  

  MCresults = new ResultDescriptorType[MCresultCount];
  MCresultCount = 0;
  
  //Reset readout-flag
  pointer = fileList.next;
  while (pointer!=NULL) {
    pointer->readout = false;
    pointer = pointer->next;
  }
  
  pointer = fileList.next;
  while (pointer!=NULL) {
    if (!pointer->readout) {
      MCresults[MCresultCount] = calculateMCresult(true,pointer->kappa,pointer->lambda,pointer->Y,pointer->rho,pointer->r, pointer->L0,pointer->L1,pointer->L2,pointer->L3,pointer->Nf);
      MCresultCount++;
    }
    pointer = pointer->next;
  }  
}


void printMCresults() {
  int I;
  printf("\nWriting MC-results to disk...\n");
  if (FLAG_tildeValues) {
    printf("\nWriting MC-results to disk (Tilde-Values)...\n");
  } else {
    printf("\nWriting MC-results to disk...\n");
  }
  FILE* file;
  file = fopen("MCresults.dat","w");
  for (I=0; I<MCresultCount; I++) {
    MCresults[I].print();
    fprintf(file,"%1.15f %1.15f %1.15f %1.15f %1.15f %d %d %d %d %d %d %1.15f %1.15f %1.15f %1.15f %1.15f %1.15f %1.15f %1.15f %1.15f\n", MCresults[I].kappa, MCresults[I].lambda, MCresults[I].Y,
     MCresults[I].rho, MCresults[I].r, MCresults[I].L0, MCresults[I].L1,MCresults[I].L2,MCresults[I].L3,MCresults[I].Nf, MCresults[I].count, MCresults[I].DSL,
     MCresults[I].averagePhiNorm, MCresults[I].sigmaPhiNorm, MCresults[I].averageStaggeredPhiNorm, MCresults[I].sigmaStaggeredPhiNorm,
     MCresults[I].Susceptibility, MCresults[I].sigmaSusceptibility, MCresults[I].StaggeredSusceptibility, MCresults[I].sigmaStaggeredSusceptibility);
  }
  fclose(file);
}


int analyzeKappaDep(bool deliverRes, ResultDescriptorType** dep, double lambda, double Y, double rho, double r, int L0, int L1, int L2, int L3, int Nf) {
  int I;
  int count = 0;
  
  for (I=0; I<MCresultCount; I++) {
    if (((MCresults[I].lambda == lambda) || (((!(lambda==lambda))) && (!(MCresults[I].lambda==MCresults[I].lambda)))) 
     && (MCresults[I].Y == Y) && (MCresults[I].rho == rho) && (MCresults[I].r == r) && (MCresults[I].L0 == L0) 
     && (MCresults[I].L1 == L1) && (MCresults[I].L2 == L2) && (MCresults[I].L3 == L3) && (MCresults[I].Nf == Nf)) {
       
      MCresults[I].readout = true; 
	
      if (deliverRes) {
	dep[count] = &(MCresults[I]);
      }
	
      count++;
    }
  }

  if (deliverRes) {
    int changed = true;
    ResultDescriptorType* dummy;
    
    while (changed) {
      changed = false;
      for (I=0; I<count-1; I++) {
        if (dep[I]->kappa>dep[I+1]->kappa) {
	  dummy = dep[I];
	  dep[I] = dep[I+1];
	  dep[I+1] = dummy;
	  changed = true;
	}
      }
    }  
  }

  return count;
}


bool findPhiPhaseTransition(ResultDescriptorType** data, int count, double thres, ResultDescriptorType& res, ResultDescriptorType& sigmaRes) {
  int I;
  
  if (count<=0) {
    res = ResultDescriptorType();
    sigmaRes = ResultDescriptorType();
    return false;
  }
  if (data[0]->averagePhiNorm >= thres) {
    res = *(data[0]);
    sigmaRes = ResultDescriptorType();   
    return false;
  }
  if (data[count-1]->averagePhiNorm <= thres) {
    res = *(data[count-1]);
    sigmaRes = ResultDescriptorType();   
    return false;   
  }
  
  for (I=0; I<count-1; I++) {
    if ((data[I]->averagePhiNorm <= thres) && (data[I+1]->averagePhiNorm >= thres)) {
      double fac = (thres-data[I]->averagePhiNorm) / (data[I+1]->averagePhiNorm-data[I]->averagePhiNorm);
      double sigmaFac = sqrt(sqr(data[I]->sigmaPhiNorm*(thres-data[I+1]->averagePhiNorm)) 
                           + sqr(data[I+1]->sigmaPhiNorm*(thres-data[I]->averagePhiNorm)) ) 
			   / sqr(data[I+1]->averagePhiNorm - data[I]->averagePhiNorm);
      
      res.kappa = data[I]->kappa + fac * (data[I+1]->kappa - data[I]->kappa);
      res.lambda = data[I]->lambda + fac * (data[I+1]->lambda - data[I]->lambda);
      res.Y = data[I]->Y + fac * (data[I+1]->Y - data[I]->Y);
      res.rho = data[I]->rho + fac * (data[I+1]->rho - data[I]->rho);
      res.r = data[I]->r + fac * (data[I+1]->r - data[I]->r);
      res.L0 = (int) (data[I]->L0 + fac * (data[I+1]->L0 - data[I]->L0));
      res.L1 = (int) (data[I]->L1 + fac * (data[I+1]->L1 - data[I]->L1));
      res.L2 = (int) (data[I]->L2 + fac * (data[I+1]->L2 - data[I]->L2));
      res.L3 = (int) (data[I]->L3 + fac * (data[I+1]->L3 - data[I]->L3));
      res.Nf = (int) (data[I]->Nf + fac * (data[I+1]->Nf - data[I]->Nf));

      sigmaRes.kappa = sigmaFac * (data[I+1]->kappa - data[I]->kappa);
      sigmaRes.lambda = sigmaFac * (data[I+1]->lambda - data[I]->lambda);
      sigmaRes.Y = sigmaFac * (data[I+1]->Y - data[I]->Y);
      sigmaRes.rho = sigmaFac * (data[I+1]->rho - data[I]->rho);
      sigmaRes.r = sigmaFac * (data[I+1]->r - data[I]->r);
      sigmaRes.L0 = (int) (sigmaFac * (data[I+1]->L0 - data[I]->L0));
      sigmaRes.L1 = (int) (sigmaFac * (data[I+1]->L1 - data[I]->L1));
      sigmaRes.L2 = (int) (sigmaFac * (data[I+1]->L2 - data[I]->L2));
      sigmaRes.L3 = (int) (sigmaFac * (data[I+1]->L3 - data[I]->L3));
      sigmaRes.Nf = (int) (sigmaFac * (data[I+1]->Nf - data[I]->Nf));
      
      //Sinnlos...
      res.averagePhiNorm = data[I]->averagePhiNorm + fac * (data[I+1]->averagePhiNorm - data[I]->averagePhiNorm);
      res.averageStaggeredPhiNorm = data[I]->averageStaggeredPhiNorm + fac * (data[I+1]->averageStaggeredPhiNorm - data[I]->averageStaggeredPhiNorm);
      res.DSL = data[I]->DSL + fac * (data[I+1]->DSL - data[I]->DSL);
      res.sigmaPhiNorm = data[I]->sigmaPhiNorm + fac * (data[I+1]->sigmaPhiNorm - data[I]->sigmaPhiNorm);
      res.sigmaStaggeredPhiNorm = data[I]->sigmaStaggeredPhiNorm + fac * (data[I+1]->sigmaStaggeredPhiNorm - data[I]->sigmaStaggeredPhiNorm);
      res.count = (int) (data[I]->count + fac * (data[I+1]->count - data[I]->count));

      return true;
    }
  }
  
  printf("ERROR in findPhiPhaseTransition!!!\n");  
  exit(0);
}


bool findStaggeredPhiPhaseTransition(ResultDescriptorType** data, int count, double thres, ResultDescriptorType& res, ResultDescriptorType& sigmaRes) {
  int I;
  
  if (count<=0) {
    res = ResultDescriptorType();
    sigmaRes = ResultDescriptorType();
    return false;
  }
  if (data[0]->averageStaggeredPhiNorm <= thres) {
    res = *(data[0]);
    sigmaRes = ResultDescriptorType();   
    return false;
  }
  if (data[count-1]->averageStaggeredPhiNorm >= thres) {
    res = *(data[count-1]);
    sigmaRes = ResultDescriptorType();   
    return false;   
  }
  
  for (I=0; I<count-1; I++) {
    if ((data[I]->averageStaggeredPhiNorm >= thres) && (data[I+1]->averageStaggeredPhiNorm <= thres)) {
      double fac = (thres-data[I]->averageStaggeredPhiNorm) / (data[I+1]->averageStaggeredPhiNorm-data[I]->averageStaggeredPhiNorm);
      double sigmaFac = sqrt(sqr(data[I]->sigmaStaggeredPhiNorm*(thres-data[I+1]->averageStaggeredPhiNorm)) 
                           + sqr(data[I+1]->sigmaStaggeredPhiNorm*(thres-data[I]->averageStaggeredPhiNorm)) ) 
			   / sqr(data[I+1]->averageStaggeredPhiNorm - data[I]->averageStaggeredPhiNorm);
      
      res.kappa = data[I]->kappa + fac * (data[I+1]->kappa - data[I]->kappa);
      res.lambda = data[I]->lambda + fac * (data[I+1]->lambda - data[I]->lambda);
      res.Y = data[I]->Y + fac * (data[I+1]->Y - data[I]->Y);
      res.rho = data[I]->rho + fac * (data[I+1]->rho - data[I]->rho);
      res.r = data[I]->r + fac * (data[I+1]->r - data[I]->r);
      res.L0 = (int) (data[I]->L0 + fac * (data[I+1]->L0 - data[I]->L0));
      res.L1 = (int) (data[I]->L1 + fac * (data[I+1]->L1 - data[I]->L1));
      res.L2 = (int) (data[I]->L2 + fac * (data[I+1]->L2 - data[I]->L2));
      res.L3 = (int) (data[I]->L3 + fac * (data[I+1]->L3 - data[I]->L3));
      res.Nf = (int) (data[I]->Nf + fac * (data[I+1]->Nf - data[I]->Nf));

      sigmaRes.kappa = sigmaFac * (data[I+1]->kappa - data[I]->kappa);
      sigmaRes.lambda = sigmaFac * (data[I+1]->lambda - data[I]->lambda);
      sigmaRes.Y = sigmaFac * (data[I+1]->Y - data[I]->Y);
      sigmaRes.rho = sigmaFac * (data[I+1]->rho - data[I]->rho);
      sigmaRes.r = sigmaFac * (data[I+1]->r - data[I]->r);
      sigmaRes.L0 = (int) (sigmaFac * (data[I+1]->L0 - data[I]->L0));
      sigmaRes.L1 = (int) (sigmaFac * (data[I+1]->L1 - data[I]->L1));
      sigmaRes.L2 = (int) (sigmaFac * (data[I+1]->L2 - data[I]->L2));
      sigmaRes.L3 = (int) (sigmaFac * (data[I+1]->L3 - data[I]->L3));
      sigmaRes.Nf = (int) (sigmaFac * (data[I+1]->Nf - data[I]->Nf));
      
      //Sinnlos...
      res.averagePhiNorm = data[I]->averagePhiNorm + fac * (data[I+1]->averagePhiNorm - data[I]->averagePhiNorm);
      res.averageStaggeredPhiNorm = data[I]->averageStaggeredPhiNorm + fac * (data[I+1]->averageStaggeredPhiNorm - data[I]->averageStaggeredPhiNorm);
      res.DSL = data[I]->DSL + fac * (data[I+1]->DSL - data[I]->DSL);
      res.sigmaPhiNorm = data[I]->sigmaPhiNorm + fac * (data[I+1]->sigmaPhiNorm - data[I]->sigmaPhiNorm);
      res.sigmaStaggeredPhiNorm = data[I]->sigmaStaggeredPhiNorm + fac * (data[I+1]->sigmaStaggeredPhiNorm - data[I]->sigmaStaggeredPhiNorm);
      res.count = (int) (data[I]->count + fac * (data[I+1]->count - data[I]->count));

      return true;
    }
  }
  
  printf("ERROR in findStaggeredPhiPhaseTransition!!!\n");  
  exit(0);
}


double FitChiPhi(double* para) {
  int I;
  double chiSqr = 0;
  for (I=0; I<FitDataCount; I++) {
    double f;
    if (FitDataResults[I]->kappa<para[0]) {
      f = para[2]*sqr(FitDataResults[I]->kappa-para[0]);
    } else {
      f = para[3]*sqr(FitDataResults[I]->kappa-para[0]);
    }
    f = para[1] * exp(-0.5*Critical_Gamma*log(f + exp(-2*log(FitDataResults[I]->Lavg)/Critical_Nu)));
    chiSqr += sqr((f-FitDataResults[I]->Susceptibility) / FitDataResults[I]->sigmaSusceptibility);
  }
  return chiSqr;
}


double FitChiStaggeredPhi(double* para) {
  int I;
  double chiSqr = 0;
  for (I=0; I<FitDataCount; I++) {
    double f;
    if (FitDataResults[I]->kappa<para[0]) {
      f = para[2]*sqr(FitDataResults[I]->kappa-para[0]);
    } else {
      f = para[3]*sqr(FitDataResults[I]->kappa-para[0]);    
    }
    f = para[1] * exp(-0.5*Critical_Gamma*log(f + exp(-2*log(FitDataResults[I]->Lavg)/Critical_Nu)));
    chiSqr += sqr((f-FitDataResults[I]->StaggeredSusceptibility) / FitDataResults[I]->sigmaStaggeredSusceptibility);
  }
  return chiSqr;
}


bool findPhiPhaseTransitionFromSusceptibility(double* para, double* error, double& redChiSqr) {
  int I;
  double high = -1E10;
  double M0,M1,M2,M3;

  for (I=0; I<FitDataCount; I++) {
    if (high<FitDataResults[I]->Susceptibility) {
      high = FitDataResults[I]->Susceptibility;
      para[0] = FitDataResults[I]->kappa+1E-6;
    }
  }
  para[2] = para[3] = 20;
  para[1] = high*exp(-(Critical_Gamma/Critical_Nu)*log(FitDataResults[0]->Lavg));
  error[0] = error[1] = error[2] = error[3] = 0;
  
  M0 = para[0];
  M1 = para[1];
  M2 = para[2];
  M3 = para[3];
  
  if (FitDataCount <= 3) return false;
 
  bool b = true;
  int count = 0;
  while (true) {
    b = true;
    para[0] = M0;
    para[1] = M1;
    para[2] = M2;
    para[3] = M3;
    
    if (count<4) b = b && GradientMinimization(&FitChiPhi, 4, 1E-4, 1E-11, 1E-12, para, NULL, NULL,NULL,4,1000);
    if (count<3) b = b && GradientMinimization(&FitChiPhi, 4, 1E-4, 1E-11, 1E-12, para, NULL, NULL,NULL,8,1000);
    if (count<2) b = b && GradientMinimization(&FitChiPhi, 4, 1E-4, 1E-11, 1E-12, para, NULL, NULL,NULL,3,1000);
    if (count<1) b = b && GradientMinimization(&FitChiPhi, 4, 1E-4, 1E-11, 1E-12, para, NULL, NULL,NULL,15,1000);
    if (b) break;
    count++;
  }
  

  if (FLAG_GnuplotFit) {
    double* x = new double[FitDataCount];
    double* y = new double[FitDataCount];
    double* yErr = new double[FitDataCount];
    for (I=0; I<FitDataCount; I++) {
      x[I] = FitDataResults[I]->kappa;
      y[I] = FitDataResults[I]->Susceptibility;
      yErr[I] = FitDataResults[I]->sigmaSusceptibility;
    }
    redChiSqr = 0;
    double NennerAdd = exp(-2*log(FitDataResults[0]->Lavg)/Critical_Nu); 

    char* functionBody = new char[1000];
    snprintf(functionBody,1000,"A2*exp(-0.5*%1.15f*log(0.5*(1-sgn(x-A1))*A3*(x-A1)**2 + 0.5*(1+sgn(x-A1))*A4*(x-A1)**2 + %1.15f))",Critical_Gamma, NennerAdd);
 
    b = /*b &&*/ performGnuplotFit(functionBody, x, y, yErr, FitDataCount, 4, para, error, redChiSqr);

    delete[] x;
    delete[] y;
    delete[] yErr;
    delete[] functionBody;
  } else {
    double* para2 = new double[4];

    for (I=0; I<FitDataCount; I++) {
      double h = 1E-2 * FitDataResults[I]->sigmaSusceptibility;
      FitDataResults[I]->Susceptibility += h;
      para2[0] = para[0];
      para2[1] = para[1];
      para2[2] = para[2];
      para2[3] = para[3];
    
      b = b && GradientMinimization(&FitChiPhi, 4, 1E-4, 1E-11, 1E-12, para2, NULL, NULL,NULL,15,1000);
      
//printf("abl: %f at kappa=%f, err=%f\n",(para2[0]-para[0])/h,FitDataResults[I]->kappa,FitDataResults[I]->sigmaSusceptibility);      
      
      error[0] += sqr(FitDataResults[I]->sigmaSusceptibility * (para2[0]-para[0])/h);
      error[1] += sqr(FitDataResults[I]->sigmaSusceptibility * (para2[1]-para[1])/h);
      error[2] += sqr(FitDataResults[I]->sigmaSusceptibility * (para2[2]-para[2])/h);
      error[3] += sqr(FitDataResults[I]->sigmaSusceptibility * (para2[3]-para[3])/h);
    
      FitDataResults[I]->Susceptibility -= h;
    }
    error[0] = sqrt(error[0]);
    error[1] = sqrt(error[1]);
    error[2] = sqrt(error[2]);
    error[3] = sqrt(error[3]);
    delete[] para2;
  }
  
  return b;
}


bool findStaggeredPhiPhaseTransitionFromSusceptibility(double* para, double* error, double& redChiSqr) {
  int I;
  double high = -1E10;
  double M1,M2,M3,M0;

  for (I=0; I<FitDataCount; I++) {
    if (high<FitDataResults[I]->StaggeredSusceptibility) {
      high = FitDataResults[I]->StaggeredSusceptibility;
      para[0] = FitDataResults[I]->kappa+1E-6;
    }
  }
  para[2] = para[3] = 20;
  para[1] = high*exp(-(Critical_Gamma/Critical_Nu)*log(FitDataResults[0]->Lavg));
  error[0] = error[1] = error[2] = error[3] = 0;
  M0 = para[0];
  M1 = para[1];
  M2 = para[2];
  M3 = para[3];
  
  
  if (FitDataCount <= 3) return false;
  
  bool b = true;
  int count = 0;
  while (true) {
    b = true;
    para[0] = M0;
    para[1] = M1;
    para[2] = M2;
    para[3] = M3;
    
    if (count<4) b = b && GradientMinimization(&FitChiStaggeredPhi, 4, 1E-4, 1E-11, 1E-12, para, NULL, NULL,NULL,4,1000);
    if (count<3) b = b && GradientMinimization(&FitChiStaggeredPhi, 4, 1E-4, 1E-11, 1E-12, para, NULL, NULL,NULL,8,1000);
    if (count<2) b = b && GradientMinimization(&FitChiStaggeredPhi, 4, 1E-4, 1E-11, 1E-12, para, NULL, NULL,NULL,3,1000);
    if (count<1) b = b && GradientMinimization(&FitChiStaggeredPhi, 4, 1E-4, 1E-11, 1E-12, para, NULL, NULL,NULL,15,1000);
    if (b) break;
    count++;
  }
  
  if (FLAG_GnuplotFit) {
    double* x = new double[FitDataCount];
    double* y = new double[FitDataCount];
    double* yErr = new double[FitDataCount];
    for (I=0; I<FitDataCount; I++) {
      x[I] = FitDataResults[I]->kappa;
      y[I] = FitDataResults[I]->StaggeredSusceptibility;
      yErr[I] = FitDataResults[I]->sigmaStaggeredSusceptibility;
    }
    redChiSqr = 0;
    double NennerAdd = exp(-2*log(FitDataResults[0]->Lavg)/Critical_Nu); 

    char* functionBody = new char[1000];
    snprintf(functionBody,1000,"A2*exp(-0.5*%1.15f*log(0.5*(1-sgn(x-A1))*A3*(x-A1)**2 + 0.5*(1+sgn(x-A1))*A4*(x-A1)**2 + %1.15f))",Critical_Gamma, NennerAdd);
  
    b = b && performGnuplotFit(functionBody, x, y, yErr, FitDataCount, 4, para, error, redChiSqr);
printf("ChiSqr: %f\n",redChiSqr);    
    delete[] x;
    delete[] y;
    delete[] yErr;
    delete[] functionBody;
  } else {
    double* para2 = new double[4];
  
    for (I=0; I<FitDataCount; I++) {
      double h = 1E-2 * FitDataResults[I]->sigmaStaggeredSusceptibility;
      FitDataResults[I]->StaggeredSusceptibility += h;
      para2[0] = para[0];
      para2[1] = para[1];
      para2[2] = para[2];
      para2[3] = para[3];
    
      b = b && GradientMinimization(&FitChiStaggeredPhi, 4, 1E-4, 1E-11, 1E-12, para2, NULL, NULL,NULL,15,1000);
      
//printf("abl: %f at kappa=%f, err=%f\n",(para2[0]-para[0])/h,FitDataResults[I]->kappa,FitDataResults[I]->sigmaStaggeredSusceptibility);      
      
      error[0] += sqr(FitDataResults[I]->sigmaStaggeredSusceptibility * (para2[0]-para[0])/h);
      error[1] += sqr(FitDataResults[I]->sigmaStaggeredSusceptibility * (para2[1]-para[1])/h);
      error[2] += sqr(FitDataResults[I]->sigmaStaggeredSusceptibility * (para2[2]-para[2])/h);
      error[3] += sqr(FitDataResults[I]->sigmaStaggeredSusceptibility * (para2[3]-para[3])/h);
    
      FitDataResults[I]->StaggeredSusceptibility -= h;
    }
    error[0] = sqrt(error[0]);
    error[1] = sqrt(error[1]);
    error[2] = sqrt(error[2]);
    error[3] = sqrt(error[3]);
  
    delete[] para2;
  }
  return b;
}


void makeKappaTransitionPoints(double thres) {
  int I,I2;
  double* x;
  double* y;
  double* y1;
  double* y2;
  double* err;
  double* err1;
  double* err2;
  char* plotCommand;
  char* plotCommand2;
  char* plotCommand3;
  double smallestX;
  double largestX;
  char* title;
  char* title1;
  char* title2;
  char* caption;
  double redChiSqr;
  
  FILE* file1 = fopen("phaseTransKappaPhi.dat","w");;
  FILE* file2 = fopen("phaseTransKappaStaggeredPhi.dat","w");;
 
  for (I=0; I<MCresultCount; I++) {
    MCresults[I].readout = false;
  }
  printf("\nCalculating kappa transition points...\n");
  for (I=0; I<MCresultCount; I++) {
    if (!MCresults[I].readout) {
      FitDataCount = analyzeKappaDep(false, NULL, MCresults[I].lambda, MCresults[I].Y, MCresults[I].rho, MCresults[I].r, MCresults[I].L0,MCresults[I].L1,MCresults[I].L2,MCresults[I].L3, MCresults[I].Nf);
      FitDataResults = new ResultDescriptorType*[FitDataCount];
      FitDataCount = analyzeKappaDep(true, FitDataResults, MCresults[I].lambda, MCresults[I].Y, MCresults[I].rho, MCresults[I].r, MCresults[I].L0,MCresults[I].L1,MCresults[I].L2,MCresults[I].L3, MCresults[I].Nf);
	
      double* fitPara = new double[4];
      double* fitParaError = new double[4];
      
      char* section = new char[1000];
      snprintf(section,1000, "Plots for $y_N = %1.3f$, $\\lambda=%1.3f$, $L=%dx%dx%dx%d$, $N_f=%d$", FitDataResults[0]->Y, FitDataResults[0]->lambda, FitDataResults[0]->L0,FitDataResults[0]->L1,FitDataResults[0]->L2,FitDataResults[0]->L3, FitDataResults[0]->Nf);
      FitControlLogger.addSection(section);
      delete[] section;
      plotCommand = new char[1000];
      plotCommand2 = new char[1000];
      plotCommand3 = new char[1000];
      x = new double[FitDataCount];
      y = new double[FitDataCount];
      err = new double[FitDataCount];
      smallestX = 1E10;
      largestX = -1E10;	
      for (I2=0; I2<FitDataCount; I2++) {
        x[I2] = FitDataResults[I2]->kappa;
        y[I2] = FitDataResults[I2]->Susceptibility;
        err[I2] = FitDataResults[I2]->sigmaSusceptibility;
        if (x[I2]<smallestX) smallestX = x[I2];
        if (x[I2]>largestX) largestX = x[I2];
      }

      if (findPhiPhaseTransitionFromSusceptibility(fitPara, fitParaError, redChiSqr)) {
        if (FLAG_TransCut == false) {
          printf("Kap: %1.3f+-%1.3f, Lam: %1.3f, yN: %1.3f, Rho: %1.3f, R: %1.3f, L: %dx%dx%dx%d, Nf: %d\n", 
           fitPara[0],fitParaError[0], FitDataResults[0]->lambda, FitDataResults[0]->Y, FitDataResults[0]->rho, FitDataResults[0]->r, 
	   FitDataResults[0]->L0, FitDataResults[0]->L1,FitDataResults[0]->L2,FitDataResults[0]->L3, FitDataResults[0]->Nf);
          fprintf(file1,"%1.15f %1.15f %1.15f %1.15f %1.15f %1.15f %d %d %d %d %d\n", 
           fitPara[0],fitParaError[0], FitDataResults[0]->lambda, FitDataResults[0]->Y, FitDataResults[0]->rho, FitDataResults[0]->r, 
	   FitDataResults[0]->L0,FitDataResults[0]->L1,FitDataResults[0]->L2,FitDataResults[0]->L3, FitDataResults[0]->Nf);
	 }
	snprintf(plotCommand,1000,"set parametric\n set trange [0:1]\n p1(x)=%f+x*(%f-%f)\np2(x)=%f+x*(%f-%f)\n",
	 smallestX,fitPara[0],smallestX,largestX,fitPara[0],largestX);
	snprintf(plotCommand2,1000,"%s\nreplot p1(t),%f*exp(-0.5*%f*log(exp(-2*log(%f)/%f) + %f*(p1(t)-%f)*(p1(t)-%f))) with lines lw 4 lt 1 notitle \n",
 	 plotCommand,fitPara[1],Critical_Gamma,FitDataResults[0]->Lavg,Critical_Nu,
	 fitPara[2],fitPara[0],fitPara[0]);
	snprintf(plotCommand3,1000,"%s\nreplot p2(t),%f*exp(-0.5*%f*log(exp(-2*log(%f)/%f) + %f*(p2(t)-%f)*(p2(t)-%f))) with lines lw 4 lt 1 notitle \n",
	 plotCommand2,fitPara[1],Critical_Gamma,FitDataResults[0]->Lavg,Critical_Nu,
	 fitPara[3],fitPara[0],fitPara[0]);
      } else {
	snprintf(plotCommand3,1000,"\n");
      }

      title = new char[1000];
      snprintf(title,1000,"Susceptibility at $y_N=%1.3f$",FitDataResults[0]->Y);
       caption = new char[1000];
      snprintf(caption,1000,"Susceptibility $\\chi$ versus $\\kappa_N$ fitted by finite size analysis with reduced $\\chi^2=%1.2f$.",redChiSqr);
	
      FitControlLogger.addPlot(title,NULL,caption,"$\\kappa_N$","Susceptibility $\\chi$",x, y, err, NULL, NULL,FitDataCount, plotCommand3);

      delete[] plotCommand;
      delete[] plotCommand2;
      delete[] plotCommand3;
      delete[] caption;
      delete[] title;
      delete[] x;
      delete[] y;
      delete[] err;


      plotCommand = new char[1000];
      plotCommand2 = new char[1000];
      plotCommand3 = new char[1000];
      x = new double[FitDataCount];
      y = new double[FitDataCount];
      err = new double[FitDataCount];
      smallestX = 1E10;
      largestX = -1E10;	
      for (I2=0; I2<FitDataCount; I2++) {
        x[I2] = FitDataResults[I2]->kappa;
        y[I2] = FitDataResults[I2]->StaggeredSusceptibility;
        err[I2] = FitDataResults[I2]->sigmaStaggeredSusceptibility;
        if (x[I2]<smallestX) smallestX = x[I2];
        if (x[I2]>largestX) largestX = x[I2];
      }
      if (findStaggeredPhiPhaseTransitionFromSusceptibility(fitPara, fitParaError,redChiSqr)) {
        if (FLAG_TransCut == false) {
          printf("Kap: %1.3f+-%1.3f, Lam: %1.3f, yN: %1.3f, Rho: %1.3f, R: %1.3f, L: %dx%dx%dx%d, Nf: %d\n", 
           fitPara[0],fitParaError[0], FitDataResults[0]->lambda, FitDataResults[0]->Y, FitDataResults[0]->rho, FitDataResults[0]->r, 
	   FitDataResults[0]->L0,FitDataResults[0]->L1,FitDataResults[0]->L2,FitDataResults[0]->L3, FitDataResults[0]->Nf);
          fprintf(file2,"%1.15f %1.15f %1.15f %1.15f %1.15f %1.15f %d %d %d %d %d\n", 
           fitPara[0],fitParaError[0], FitDataResults[0]->lambda, FitDataResults[0]->Y, FitDataResults[0]->rho, FitDataResults[0]->r, 
	   FitDataResults[0]->L0,FitDataResults[0]->L1,FitDataResults[0]->L2,FitDataResults[0]->L3, FitDataResults[0]->Nf);
	}
  	snprintf(plotCommand,1000,"set parametric\n set trange [0:1]\n p1(x)=%f+x*(%f-%f)\np2(x)=%f+x*(%f-%f)\n",
	 smallestX,fitPara[0],smallestX,largestX,fitPara[0],largestX);
        snprintf(plotCommand2,1000,"%s\nreplot p1(t),%f*exp(-0.5*%f*log(exp(-2*log(%f)/%f) + %f*(p1(t)-%f)*(p1(t)-%f))) with lines lw 4 lt 1 notitle \n",
	 plotCommand,fitPara[1],Critical_Gamma,FitDataResults[0]->Lavg,Critical_Nu,
	 fitPara[2],fitPara[0],fitPara[0]);
	snprintf(plotCommand3,1000,"%s\nreplot p2(t),%f*exp(-0.5*%f*log(exp(-2*log(%f)/%f) + %f*(p2(t)-%f)*(p2(t)-%f))) with lines lw 4 lt 1 notitle \n",
	 plotCommand2,fitPara[1],Critical_Gamma,FitDataResults[0]->Lavg,Critical_Nu,
	 fitPara[3],fitPara[0],fitPara[0]);
      } else {
        snprintf(plotCommand3,1000,"\n");
      }

      title = new char[1000];
      snprintf(title,1000,"Staggered susceptibility at $y_N=%1.3f$",FitDataResults[0]->Y);
      caption = new char[1000];
      snprintf(caption,1000,"Staggered susceptibility $\\chi$ versus $\\kappa_N$ fitted by finite size analysis with reduced $\\chi^2=%1.2f$.",redChiSqr);
	
      FitControlLogger.addPlot(title,NULL,caption,"$\\kappa_N$","Staggered susceptibility $\\chi$",x, y, err, NULL, NULL, FitDataCount, plotCommand3);

      delete[] plotCommand;
      delete[] plotCommand2;
      delete[] plotCommand3;
      delete[] caption;
      delete[] title;
      delete[] x;
      delete[] y;
      delete[] err;
      
      delete[] fitPara;
      delete[] fitParaError;

      ResultDescriptorType transP;
      ResultDescriptorType sigmaTransP;

      if (findPhiPhaseTransition(FitDataResults, FitDataCount, thres, transP, sigmaTransP)) {
        if (FLAG_TransCut == true) {
          printf("Kap: %1.3f+-%1.3f, Lam: %1.3f+-%1.3f, yN: %1.3f+-%1.3f, Rho: %1.3f+-%1.3f, R: %1.3f+-%1.3f, L: %dx%dx%dx%d, Nf: %d, Stat: %d, DSL: %d\n", 
           transP.kappa, sigmaTransP.kappa, transP.lambda, sigmaTransP.lambda, transP.Y, sigmaTransP.Y, transP.rho, sigmaTransP.rho, 
           transP.r, sigmaTransP.r, transP.L0, transP.L1,transP.L2,transP.L3,transP.Nf, transP.count, (int)(transP.DSL/log(10)));
          fprintf(file1,"%1.15f %1.15f %1.15f %1.15f %1.15f %1.15f %1.15f %1.15f %1.15f %1.15f %d %d %d %d %d %d %d %d %d %d\n", 
           transP.kappa, sigmaTransP.kappa, transP.lambda, sigmaTransP.lambda, transP.Y, sigmaTransP.Y, transP.rho, sigmaTransP.rho, 
           transP.r, sigmaTransP.r, transP.L0, sigmaTransP.L0, transP.L1, sigmaTransP.L1, transP.L2, sigmaTransP.L2,
	   transP.L3, sigmaTransP.L3, transP.Nf, sigmaTransP.Nf);
	 }
      }


      if (findStaggeredPhiPhaseTransition(FitDataResults, FitDataCount, thres, transP, sigmaTransP)) {
        if (FLAG_TransCut == true) {
          printf("Kap: %1.3f+-%1.3f, Lam: %1.3f+-%1.3f, yN: %1.3f+-%1.3f, Rho: %1.3f+-%1.3f, R: %1.3f+-%1.3f, L: %dx%dx%dx%d, Nf: %d, Stat: %d, DSL: %d\n", 
           transP.kappa, sigmaTransP.kappa, transP.lambda, sigmaTransP.lambda, transP.Y, sigmaTransP.Y, transP.rho, sigmaTransP.rho, 
           transP.r, sigmaTransP.r, transP.L0, transP.L1, transP.L2, transP.L3, transP.Nf, transP.count, (int)(transP.DSL/log(10)));
          fprintf(file2,"%1.15f %1.15f %1.15f %1.15f %1.15f %1.15f %1.15f %1.15f %1.15f %1.15f %d %d %d %d %d %d %d %d %d %d\n", 
           transP.kappa, sigmaTransP.kappa, transP.lambda, sigmaTransP.lambda, transP.Y, sigmaTransP.Y, transP.rho, sigmaTransP.rho, 
           transP.r, sigmaTransP.r, transP.L0, sigmaTransP.L0, transP.L1, sigmaTransP.L1, transP.L2, sigmaTransP.L2, 
	   transP.L3, sigmaTransP.L3, transP.Nf, sigmaTransP.Nf);
	}
      }

      x = new double[FitDataCount];
      y1 = new double[FitDataCount];
      err1 = new double[FitDataCount];
      y2 = new double[FitDataCount];
      err2 = new double[FitDataCount];
      for (I2=0; I2<FitDataCount; I2++) {
        x[I2] = FitDataResults[I2]->kappa;
	y1[I2] = FitDataResults[I2]->averagePhiNorm;
	err1[I2] = FitDataResults[I2]->sigmaPhiNorm;
	y2[I2] = FitDataResults[I2]->averageStaggeredPhiNorm;
	err2[I2] = FitDataResults[I2]->sigmaStaggeredPhiNorm;
      }
      plotCommand = new char[1000];
      if (FLAG_TransCut == true) {
        snprintf(plotCommand,1000,"replot %f title 'Cut value = %f'\n", thres, thres);  
      } else {
        snprintf(plotCommand,1000,"\n");  
      }

      title1 = new char[1000];
      title2 = new char[1000];
      snprintf(title1,1000,"Average $|\\Phi|$ at $\\tilde y_N=%1.3f$",FitDataResults[0]->Y);
      snprintf(title2,1000,"Average staggered $|\\Phi|$ at $\\tilde y_N=%1.3f$",FitDataResults[0]->Y);
      caption = new char[1000];
      snprintf(caption,1000,"(Staggered) average $|\\Phi|$ versus $\\kappa_N$ for $y_N=%1.3f$, $\\lambda_N=%1.3f$, $L=%dx%dx%dx%d$, $N_f=%d$.",FitDataResults[0]->Y,FitDataResults[0]->lambda, FitDataResults[0]->L0,FitDataResults[0]->L1,FitDataResults[0]->L2,FitDataResults[0]->L3,FitDataResults[0]->Nf);
    
      if (FLAG_FiniteSizeEffect) {
        double* fse1 = new double[FitDataCount];
        double* fse2 = new double[FitDataCount];
      
	int I3;
	for (I3=0; I3<FitDataCount; I3++) {
          findFiniteVolumeEffectiveActionGroundState(FitDataResults[0]->Nf, FitDataResults[0]->L0, FitDataResults[0]->L1, FitDataResults[0]->L2, FitDataResults[0]->L3, FitDataResults[I3]->kappa, FitDataResults[0]->lambda, fse1[I3], fse2[I3]);
	}
	
        FitControlLogger.addPlotWithAdditionalLinePlots(title1, title2, caption,"$\\kappa_N$", "(Staggered) $|\\Phi|$", x, y1, err1, y2, err2, fse1, fse2, FitDataCount, plotCommand, "");
        FitControlLogger.clearPage();
        FitControlLogger.newPage();

        delete[] fse1;
        delete[] fse2;
      } else {
        FitControlLogger.addPlot(title1, title2, caption,"$\\kappa_N$", "(Staggered) $|\\Phi|$", x, y1, err1, y2, err2, FitDataCount, plotCommand);
        FitControlLogger.clearPage();
        FitControlLogger.newPage();
      }
      
      delete[] plotCommand;
      delete[] caption;
      delete[] title1;
      delete[] title2;
      delete[] x;
      delete[] y1;
      delete[] err1;
      delete[] y2;
      delete[] err2;
     
      delete[] FitDataResults;
      
    }
  }
  fclose(file1);
  fclose(file2);
}


void printPhaseTransitionPoints() {
  FILE* file = fopen("phaseTransKappaPhi.dat","r");;
  printf("\n\n   *** Phi - Phase-Transition-Points ***\n");
  while (true) {
    double d1,d2,d3,d4,d5,d6;
    int i10,i11,i12,i13,i2;
    if (fscanf(file,"%lf %lf %lf %lf %lf %lf %d %d %d %d %d",&d1,&d2,&d3,&d4,&d5,&d6,&i10,&i11,&i12,&i13,&i2) != 11) break;
    printf("Kap: %1.3f+-%1.3f, Lam: %1.3f, yN: %1.3f, Rho: %1.3f, R: %1.3f, L: %dx%dx%dx%d, Nf: %d\n",d1,d2,d3,d4,d5,d6, i10,i11,i12,i13,i2);
  }
  fclose(file);

  file = fopen("phaseTransKappaStaggeredPhi.dat","r");;
  printf("\n\n   *** Staggered-Phi - Phase-Transition-Points ***\n");
  while (true) {
    double d1,d2,d3,d4,d5,d6;
    int i10,i11,i12,i13,i2;
    if (fscanf(file,"%lf %lf %lf %lf %lf %lf %d  %d %d %d%d",&d1,&d2,&d3,&d4,&d5,&d6,&i10,&i11,&i12,&i13,&i2) != 11) break;
    printf("Kap: %1.3f+-%1.3f, Lam: %1.3f, yN: %1.3f, Rho: %1.3f, R: %1.3f, L: %dx%dx%dx%d, Nf: %d\n",d1,d2,d3,d4,d5,d6, i10,i11,i12,i13,i2);
  }

  fclose(file);
}


void makeFiniteVolumeKappaPlot(int Nf, int L0, int L1, int L2, int L3, double lambdaTilde) {
  double* k = NULL;
  double* m = NULL;
  double* s = NULL;
  double* err = new double[100];

  makeFiniteVolumeKappaFunction(Nf, L0, L1, L2, L3, lambdaTilde, -0.05, 0.05, 100, k,  m, s);
  int I;
  for (I=0; I<100; I++) {
    err[I] = 0;
  }

  ControlLogger finiteVol;
  finiteVol.setLogging(true);
  finiteVol.setLatexFileName("FiniteVolume");
  finiteVol.addPlot("", "", "", "", "", k, m, err, s, err, 100, "");

  finiteVol.generateLatex();

  delete[] k;
  delete[] m;
  delete[] s;
  delete[] err;
}


void executeCommand(int commandNr, char* cmd) {
  if (commandNr == 0) {
    FLAG_nozeros = true;
  }
  if (commandNr == 1) {
    FLAG_nochi = true;
  }
  if (commandNr == 2) {
    FLAG_wilson = true;
  }  
  if (commandNr == 3) {
    FLAG_bootStrap = true;
    bootStrapIter = 100;
    int v=0;
    if (sscanf(cmd,"boot=%d",&v)==1) {
      bootStrapIter = v;
    }
  }  
  if (commandNr == 4) {
    extrapolateToNf = 1;
    FLAG_ExtrapolateNf = true;
    int v=0;
    if (sscanf(cmd,"exNf=%d",&v)==1) {
      extrapolateToNf = v;
    }
  }  
  if (commandNr == 4) {
    extrapolateToNf = 1;
    FLAG_ExtrapolateNf = true;
    int v=0;
    if (sscanf(cmd,"exNf=%d",&v)==1) {
      extrapolateToNf = v;
    }
  }  
  if (commandNr == 5) {
    transitionCutValue = 0.5;
    double v=0;
    FLAG_TransCut = true;
    if (sscanf(cmd,"TransCut=%lf",&v)==1) {
      transitionCutValue = v;
    }
  }  
  if (commandNr == 6) {
    FLAG_hybrid = true;
  }  
  if (commandNr == 7) {
    FLAG_tildeValues = true;
  }  
  if (commandNr == 8) {
    FLAG_ImpPhiDyn = true;
  }  
  if (commandNr == 9) {
    FLAG_TexLog = true;
  }  
  if (commandNr == 10) {
    FLAG_GnuplotFit = true;
  }  
  if (commandNr == 11) {
    FLAG_MCTrajectoryLog = true;
  }  
  if (commandNr == 12) {
    FLAG_AutoCorrLog = true;
  }  
  if (commandNr == 13) {
    FLAG_FiniteSizeEffect = true;
  }  
  if (commandNr == 14) {
    removeThermal = 0;
    int v = 0;
    if (sscanf(cmd,"removeT=%d",&v)==1) {
      removeThermal = v;
    }
  }  
  if (commandNr == 15) {
    FLAG_Test = true;
  }  
  if (commandNr == 16) {
    FLAG_pHMC = true;
  }  
}


int main(int argc,char **argv) {
iniTools(5517);
  int I;

  printf("\n                                      *** Evaluating MC - Data *** \n\n");
  for (I=1; I<argc; I++) {
    Parameter_yN_Selector.readSpecification(argv[I]);
    Parameter_rho_Selector.readSpecification(argv[I]);
    Parameter_r_Selector.readSpecification(argv[I]);
    Parameter_kappa_Selector.readSpecification(argv[I]);
    Parameter_lambda_Selector.readSpecification(argv[I]);
    Parameter_L0_Selector.readSpecification(argv[I]);
    Parameter_L1_Selector.readSpecification(argv[I]);
    Parameter_L2_Selector.readSpecification(argv[I]);
    Parameter_L3_Selector.readSpecification(argv[I]);
    Parameter_Nf_Selector.readSpecification(argv[I]);
  }
  
  CommandCount = 17;
  commands = new char*[CommandCount];
  for (I=0; I<CommandCount; I++) commands[I] = new char[100];
  commands[0] = "nozeros";
  commands[1] = "nochi";
  commands[2] = "wilson";  
  commands[3] = "boot";  
  commands[4] = "exNf";  
  commands[5] = "TransCut";  
  commands[6] = "hybrid";  
  commands[7] = "tilde";  
  commands[8] = "ImpPhiDyn";  
  commands[9] = "texLog";  
  commands[10] = "gnuFit";  
  commands[11] = "traLog";  
  commands[12] = "autoLog";  
  commands[13] = "finiteSE";  
  commands[14] = "removeT";  
  commands[15] = "test";  
  commands[16] = "pHMC";  

  for (I=1; I<argc; I++) {
    int I2;
    for (I2=0; I2<CommandCount; I2++) {
      if (strncmp(commands[I2],argv[I],4)==0) {
        executeCommand(I2,argv[I]);
      }
    }
  }
  FitControlLogger.setLogging(FLAG_TexLog);
  FitControlLogger.setLatexFileName("FitControlLog");
  MCTrajectoryLogger.setLogging(FLAG_MCTrajectoryLog);
  MCTrajectoryLogger.setLatexFileName("MCTrajectoryLog");
  AutoCorrLogger.setLogging(FLAG_AutoCorrLog);
  AutoCorrLogger.setLatexFileName("AutoCorrelationLog");
  
//makeFiniteVolumeKappaPlot(10, 8, 0.1);  
//exit(0);

  loadFileList();
  printf("\n");
  readMCresults();
  
  makeKappaTransitionPoints(transitionCutValue);
  FitControlLogger.generateLatex();
  MCTrajectoryLogger.generateLatex();  
  AutoCorrLogger.generateLatex();  
  printMCresults();
  printPhaseTransitionPoints();
}
