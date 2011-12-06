#include "SimulationParameterSet.h"
#include "Tools.h"

void SimulationParameterSet::ini(double kapOrMass, double lam, double py, int nf, int type) {
  if (LogLevel>2) printf("SimulationParameterSet initialized with kappaOrMassSquared=%f, lambda=%f, y=%f, Nf=%d, and type=%d\n", kapOrMass, lam, py, nf, type);
  kappaOrMassSquared = kapOrMass;
  lambda = lam;
  y = py;
  Nf = nf;
  ParameterType = type;
  if ((ParameterType<1) || (ParameterType>3)) {
    printf("ERROR in SimulationParameterSet Constructor: Invalid Parameter-Type\n");
    exit(0);  
  }
}


SimulationParameterSet::SimulationParameterSet(double kapOrMass, double lam, double py, int nf, int type) {
  ini(kapOrMass, lam, py, nf, type);
}


SimulationParameterSet::SimulationParameterSet(StateDescriptorReader* SDreader, int type) {
  if (type==SimulationParameterSet_NfNotation) {
    ini(SDreader->getKappa(), SDreader->getLambda(), SDreader->getYN(), SDreader->getNf(), type);
  } else {
    printf("SimulationParameterSet cannot be constructed from StateDescriptorReader for Type=%d\n",type);
  }
}


SimulationParameterSet::~SimulationParameterSet() {
}


double SimulationParameterSet::reparametrize_HiggsField(int targetType) {
  if ((ParameterType==SimulationParameterSet_NfNotation) && (targetType==SimulationParameterSet_ContinuumNotation)) {
    return sqrt(2*kappaOrMassSquared);  
  }
  if ((ParameterType==SimulationParameterSet_ContinuumNotation) && (targetType==SimulationParameterSet_NfNotation)) {
    double d = 2*getKappaN();  
    return 1/sqrt(d);  
  }
  
  return NaN;
}


double SimulationParameterSet::getY0() {
  if (ParameterType==SimulationParameterSet_NfNotation) {
    return y / sqrt(2*kappaOrMassSquared);  
  }
  return NaN;
}


double SimulationParameterSet::getLambda0() {
  if (ParameterType==SimulationParameterSet_NfNotation) {
    return lambda / sqr(2*kappaOrMassSquared);  
  }
  return NaN;
}


double SimulationParameterSet::getM0Squared() {
  if (ParameterType==SimulationParameterSet_NfNotation) {
    return (1 - 2*Nf*lambda - 8*kappaOrMassSquared) / kappaOrMassSquared;
  }
  return NaN;
}


double SimulationParameterSet::getYN() {
  if (ParameterType==SimulationParameterSet_NfNotation) {
    return y;
  }
  if (ParameterType==SimulationParameterSet_ContinuumNotation) {
    return y*sqrt(2*getKappaN());
  }
  return NaN;
}


double SimulationParameterSet::getLambdaN() {
  if (ParameterType==SimulationParameterSet_NfNotation) {
    return lambda;
  }
  if (ParameterType==SimulationParameterSet_ContinuumNotation) {
    double d = 2*getKappaN();
    return lambda*d*d;
  }
  return NaN;
}


double SimulationParameterSet::getKappaN() {
  if (ParameterType==SimulationParameterSet_NfNotation) {
    return kappaOrMassSquared;
  }
  if (ParameterType==SimulationParameterSet_ContinuumNotation) {
    double d1 = 8+kappaOrMassSquared;
    double d2 = 16*Nf*lambda;
    if (lambda==0) {
      return 1.0/d1;
    }
    return (sqrt(d1*d1+2*d2) - d1) / d2;
  }
  
  return NaN;
}
