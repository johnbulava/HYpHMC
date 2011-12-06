#include "Force.h"

Force::Force(FermionMatrixOperations* fOps, bool loc, int id) {
  if (LogLevel>2) printf("Initializing Force-Calculator with LocalMode = %d and ID = %d\n", loc, id);
  local = loc;
  ID = id;
  fermiOps = fOps;  

  if (fermiOps->getYN()<=0) {
    printf("Force-Calculator: ERROR: yN=0!!!\n");
    exit(0);
  }
  if (local) {
    omegaField = fermiOps->createFermionVector();
  } else {
    omegaField = NULL;
  }
  dSdPhi = (vector4D*) fermiOps->createFermionVector(2);
}



void Force::desini() {
  if (LogLevel>0) printf("Desinitializing HMC-Force-Calculator...");
  fermiOps->destroyFermionVector(omegaField);
  
  Complex* c1 = (Complex*) dSdPhi;
  fermiOps->destroyFermionVector(c1);
  dSdPhi = NULL;

  desiniAdditionalFields();
  
  if (LogLevel>0) printf("sucessfully.\n");      
}



vector4D* Force::getdSdPhi() {
  return dSdPhi;
}


Complex* Force::getOmega() {
  return omegaField;
}


bool Force::isLocal() {
  return local;
}


void Force::writeOmegaFieldToDisk(char *confFileName) {
  if (LogLevel>2) printf("Writing Omega-Field to file %s on Force %d on node %d:\n",confFileName,ID, ownNodeID);
  std::fstream confFile;
  confFile.open(confFileName, std::ios::out);
  int VL = fermiOps->getVectorLength();

  fermiOps->transformFromXtraSizeArray(omegaField, omegaField);  

  int I;
  double avgNorm = 0;
  for (I=0; I<VL; I++) {
    avgNorm += sqr(omegaField[I].x) + sqr(omegaField[I].y); 
  }
  avgNorm = sqrt(avgNorm/VL);
  if (LogLevel>0) printf("Writing Omega Field to disk (avgNorm=%f) on force %d on node %d to file %s...\n",avgNorm,ID,ownNodeID,confFileName);   
  
  confFile.write((char*)omegaField, 16*VL);
  confFile.flush();
  confFile.close();

  fermiOps->transformToXtraSizeArray(omegaField, omegaField);  
}


bool Force::readOmegaFieldFromDisk(char *confFileName) {
  std::fstream confFile;
  confFile.open(confFileName, std::ios::in);
  int VL = fermiOps->getVectorLength();

  confFile.read((char*)omegaField, 16*VL);
  int I;
  double avgNorm = 0;
  for (I=0; I<VL; I++) {
    avgNorm += sqr(omegaField[I].x) + sqr(omegaField[I].y); 
  }
  avgNorm = sqrt(avgNorm/VL);
    
  int dataRead = 0;
  dataRead = confFile.gcount();
  bool ok = true;
  if (dataRead != 16*VL) {
    if (LogLevel>0) printf("Error while reading Omega Field on force %d on node %d from file %s...\n",ID,ownNodeID,confFileName);  
    ok = false;
  } else {
    if (LogLevel>0) printf("Omega Field successfully read (avgNorm=%f) on force %d on node %d from file %s...\n",avgNorm,ID,ownNodeID,confFileName);   
  }
  confFile.close();
  
  fermiOps->transformToXtraSizeArray(omegaField, omegaField);    
  
  return ok;
}


void Force::setOmegaFieldToZero() {
  if (local) {
    if (LogLevel>3) printf("Setting Omega-Field to zero in Force-Calculator %d on node %d:\n",ID, ownNodeID);
    int VLxtr = fermiOps->getVectorLengthXtrSize();
    int I;
    for (I=0; I<VLxtr; I++) {
      omegaField[I].x = 0;
      omegaField[I].y = 0;
    }
  }
}
