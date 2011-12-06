#ifndef SimulationParameterSet_included
#define SimulationParameterSet_included

#include <stdlib.h>

#include "Global.h"
#include "StateDescriptorReader.h"


#define SimulationParameterSet_ContinuumNotation 1
#define SimulationParameterSet_NfNotation 2 
#define SimulationParameterSet_TildeNotation 3 



class SimulationParameterSet {
private:
  double kappaOrMassSquared;
  double lambda;
  double y;
  int Nf;
  int ParameterType;

  void ini(double kapOrMass, double lam, double py, int nf, int type);
    
public:
  SimulationParameterSet(double kapOrMass, double lam, double py, int nf, int type);
  SimulationParameterSet(StateDescriptorReader* SDreader, int type);
  ~SimulationParameterSet();

  double reparametrize_HiggsField(int targetType);
  double getY0();
  double getLambda0();
  double getM0Squared();
  double getYN();
  double getLambdaN();
  double getKappaN();

};

#endif
