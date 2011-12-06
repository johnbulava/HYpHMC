#ifndef Force_included
#define Force_included

#include <math.h>
#include <fstream>
#include "Complex.h"
#include "Quat.h"
#include "Global.h"
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include "Tools.h"
#include "FermionMatrixOperations.h"



class Force {
protected:
  FermionMatrixOperations* fermiOps;
  Complex* omegaField;
  vector4D* dSdPhi;
  bool local;
  
  virtual void iniAdditionalFields() {};
  virtual void desiniAdditionalFields() {};
  int ID;
  

public:
  Force(FermionMatrixOperations* fOps, bool loc, int id); 
  virtual ~Force() {};
  void desini();

  vector4D* getdSdPhi();
  Complex* getOmega();
  bool isLocal();
  void writeOmegaFieldToDisk(char *confFileName);
  bool readOmegaFieldFromDisk(char *confFileName);
  void setOmegaFieldToZero();
};


#endif
