#include "MassCorrelationMatrixAnalyzer.h"

MassCorrelationMatrixAnalyzer::MassCorrelationMatrixAnalyzer(int dim, int timeEx, bool opCorMode, bool genEWprob, bool iterEValign, bool SubVac, bool FitVac, int fitRangeRed, int datMx, char* uniqueID) {
  if (LogLevel>2) printf("Initializing MassCorrelationMatrixAnalyzer with dim = %d, timeEx = %d, opCorMode = %d, genEWprob = %d, iterEValign = %d , SubVac = %d, FitVac = %d, fitRangeRed = %d, dataMax = %d and uniqueID = %s\n", dim, timeEx, opCorMode, genEWprob, iterEValign, SubVac, FitVac, fitRangeRed, datMx, uniqueID);

  dimension = dim;
  timeExtent = timeEx;
  fittingRangeReduction = fitRangeRed;
  
  uniqueAnalyzerFileNameExtension = new char[1000];
  snprintf(uniqueAnalyzerFileNameExtension,1000,"%s",uniqueID);
  operatorCorrelationMode = opCorMode;
  generalEWproblem = genEWprob;
  iterativeEigenVectorAlignment = iterEValign;
  subtractVacuumExpectation = SubVac;  
  fitVacuumExpectation = FitVac;    
  dataMax = datMx;
  dataCount = 0;
  dataID = new int[dataMax];
  dataWeights = new double[dataMax];
  MassCorrelationEigenvalues = new Complex*[timeExtent+1];
  MassCorrelationEigenvaluesError = new Complex*[timeExtent+1];
  for (int I=0; I<timeExtent+1; I++) {
    MassCorrelationEigenvalues[I] = new Complex[dimension];
    MassCorrelationEigenvaluesError[I] = new Complex[dimension];
  }
  effectiveMasses = new double*[timeExtent];
  effectiveMassesError = new double*[timeExtent];
  for (int I=0; I<timeExtent; I++) {
    effectiveMasses[I] = new double[dimension];
    effectiveMassesError[I] = new double[dimension];
  }
  asymptoticEffectiveMasses = new double[dimension];
  asymptoticEffectiveMassesError = new double[dimension];    
  fittedMasses = new double*[dimension];
  fittedMassesError = new double*[dimension];
  fittedCoefficients = new double*[dimension];
  fittedCoefficientsError = new double*[dimension];
  for (int I=0; I<dimension; I++) {
    fittedMasses[I] = NULL;
    fittedMassesError[I] = NULL;
    fittedCoefficients[I] = NULL;
    fittedCoefficientsError[I] = NULL;
  }
  fittedVacuumExpectationValue = new double[dimension];
  fittedVacuumExpectationValueError = new double[dimension];
  fitChiSquare = new double[dimension];
  fitChiSquareError = new double[dimension];
  
  OperatorData = NULL;
  if (!operatorCorrelationMode) {
    OperatorData = new Complex**[dataMax];
  }
  OperatorCorrelationData = NULL;
  if (operatorCorrelationMode) {
    OperatorData = new Complex**[dataMax]; 
    OperatorCorrelationData = new Complex***[dataMax];
  }
  MassCorrelationMatrices = new ComplexMatrix*[timeExtent+1];
  for (int I=0; I<timeExtent+1; I++) {
    MassCorrelationMatrices[I] = new ComplexMatrix(dimension);
  }
  GeneralEWNormalizationMatrix = NULL;
  clearData();
}


MassCorrelationMatrixAnalyzer::~MassCorrelationMatrixAnalyzer() {
  if (LogLevel>2) printf("De-Initializing MassCorrelationMatrixAnalyzer\n");

  clearData();
  delete[] uniqueAnalyzerFileNameExtension;
  delete[] dataID;
  delete[] dataWeights;
  if (!operatorCorrelationMode) {
    delete[] OperatorData;
  }
  if (operatorCorrelationMode) {
    delete[] OperatorData;
    delete[] OperatorCorrelationData;
  }
  OperatorData = NULL;
  OperatorCorrelationData = NULL;

  for (int I=0; I<timeExtent+1; I++) {
    delete[] MassCorrelationEigenvalues[I];
    delete[] MassCorrelationEigenvaluesError[I];
  }
  delete[] MassCorrelationEigenvalues;
  delete[] MassCorrelationEigenvaluesError;
  for (int I=0; I<timeExtent; I++) {
    delete[] effectiveMasses[I];
    delete[] effectiveMassesError[I];
  }
  delete[] effectiveMasses;
  delete[] effectiveMassesError;  
  delete[] asymptoticEffectiveMasses;
  delete[] asymptoticEffectiveMassesError;  

  for (int I=0; I<timeExtent+1; I++) {
    delete MassCorrelationMatrices[I];
  }
  delete[] MassCorrelationMatrices;
  delete GeneralEWNormalizationMatrix;
  for (int I=0; I<dimension; I++) {
    delete[] fittedMasses[I];
    delete[] fittedMassesError[I];
    delete[] fittedCoefficients[I];
    delete[] fittedCoefficientsError[I];
  }
  delete[] fittedMasses;
  delete[] fittedMassesError;
  delete[] fittedCoefficients;
  delete[] fittedCoefficientsError;
  delete[] fittedVacuumExpectationValue;
  delete[] fittedVacuumExpectationValueError;
  delete[] fitChiSquare;
  delete[] fitChiSquareError;
}


void MassCorrelationMatrixAnalyzer::clearData() {
  if (LogLevel>2) printf("Clearing Data in MassCorrelationMatrixAnalyzer\n");
  if (!operatorCorrelationMode) {
    for (int I=0; I<dataCount; I++) {
      for (int I2=0; I2<timeExtent; I2++) {
        delete[] OperatorData[I][I2];
      }
      delete[] OperatorData[I];
    }
  }
  if (operatorCorrelationMode) {
    for (int I=0; I<dataCount; I++) {
      for (int I2=0; I2<2*timeExtent; I2++) {
        delete[] OperatorData[I][I2];
      }
      delete[] OperatorData[I];
      for (int I2=0; I2<timeExtent+1; I2++) {
        for (int I3=0; I3<dimension; I3++) {
          delete[] OperatorCorrelationData[I][I2][I3];
        }
        delete[] OperatorCorrelationData[I][I2];
      }
      delete[] OperatorCorrelationData[I];
    }
  }
  dataCount = 0;
  for (int I=0; I<timeExtent+1; I++) {
    MassCorrelationMatrices[I]->setZero();
  }
  for (int I=0; I<timeExtent+1; I++) {
    for (int I2=0; I2<dimension; I2++) {
      MassCorrelationEigenvalues[I][I2].x = 0;
      MassCorrelationEigenvalues[I][I2].y = 0;
      MassCorrelationEigenvaluesError[I][I2].x = 0;
      MassCorrelationEigenvaluesError[I][I2].y = 0;
    }
  }
  for (int I=0; I<timeExtent; I++) {
    for (int I2=0; I2<dimension; I2++) {
      effectiveMasses[I][I2] = 0;
      effectiveMassesError[I][I2] = 0;
    }
  }
  for (int I2=0; I2<dimension; I2++) {
    asymptoticEffectiveMasses[I2] = 0;
    asymptoticEffectiveMassesError[I2] = 0;
  }

  for (int I=0; I<dimension; I++) {
    fittedVacuumExpectationValue[I] = 0;
    fittedVacuumExpectationValueError[I] = 0;
    delete[] fittedMasses[I];
    delete[] fittedMassesError[I];
    delete[] fittedCoefficients[I];
    delete[] fittedCoefficientsError[I];
    fittedMasses[I] = NULL;
    fittedMassesError[I] = NULL;
    fittedCoefficients[I] = NULL;
    fittedCoefficientsError[I] = NULL;
    fitChiSquare[I] = 0;
    fitChiSquareError[I] = 0;
  }
  fittedMassCount = 0;
  delete GeneralEWNormalizationMatrix;
  GeneralEWNormalizationMatrix = NULL;
}


void MassCorrelationMatrixAnalyzer::addOperatorData(int ID, double weight, Complex** opData) {
  if (dataCount>=dataMax) {
    printf("ERROR: DataMax = %d too small in MassCorrelationMatrixAnalyzer\n",dataMax);
    exit(0);
  }
  if (operatorCorrelationMode) {
    printf("ERROR: operatorCorrelationMode activated but addOperatorData called!!!\n");
    exit(0);
  }

  if (LogLevel>3) printf("Adding Data to MassCorrelationMatrixAnalyzer with ID = %d and weight = %f\n",ID,weight);
  dataID[dataCount] = ID;
  dataWeights[dataCount] = weight;
  OperatorData[dataCount] = new Complex*[timeExtent];
  for (int I=0; I<timeExtent; I++) {
    OperatorData[dataCount][I] = new Complex[dimension];
    for (int I2=0; I2<dimension; I2++) {
      OperatorData[dataCount][I][I2].x = opData[I][I2].x;
      OperatorData[dataCount][I][I2].y = opData[I][I2].y;
    }
  }
  dataCount++;
}


void MassCorrelationMatrixAnalyzer::addOperatorCorrelationData(int ID, double weight, Complex*** opCorData) {
  if (dataCount>=dataMax) {
    printf("ERROR: DataMax = %d too small in MassCorrelationMatrixAnalyzer\n",dataMax);
    exit(0);
  }
  if (!operatorCorrelationMode) {
    printf("ERROR: operatorCorrelationMode deactivated but addOperatorCorrelationData called!!!\n");
    exit(0);
  }

  if (LogLevel>3) printf("Adding Data to MassCorrelationMatrixAnalyzer with ID = %d and weight = %f\n",ID,weight);
  dataID[dataCount] = ID;
  dataWeights[dataCount] = weight;
  OperatorCorrelationData[dataCount] = new Complex**[timeExtent+1];
  for (int I=0; I<timeExtent+1; I++) {
    OperatorCorrelationData[dataCount][I] = new Complex*[dimension];
    for (int I2=0; I2<dimension; I2++) {
      OperatorCorrelationData[dataCount][I][I2] = new Complex[dimension];
      for (int I3=0; I3<dimension; I3++) {
        OperatorCorrelationData[dataCount][I][I2][I3].x = opCorData[I][I2][I3].x;
        OperatorCorrelationData[dataCount][I][I2][I3].y = opCorData[I][I2][I3].y;
      }
    }
  }
  OperatorData[dataCount] = new Complex*[2*timeExtent];
  for (int I=0; I<2*timeExtent; I++) {
    OperatorData[dataCount][I] = new Complex[dimension];
    for (int I2=0; I2<dimension; I2++) {
      OperatorData[dataCount][I][I2].x = 0;
      OperatorData[dataCount][I][I2].y = 0;
    }
  }
  
  dataCount++;
}


void MassCorrelationMatrixAnalyzer::addOperatorCorrelationData(int ID, double weight, Complex*** opCorData, Complex** opDagData, Complex** opData) {
  addOperatorCorrelationData(ID, weight, opCorData);
  
  for (int I=0; I<timeExtent; I++) {
    for (int I2=0; I2<dimension; I2++) {
      OperatorData[dataCount-1][I][I2] = opDagData[I][I2];
    }
  }
  for (int I=0; I<timeExtent; I++) {
    for (int I2=0; I2<dimension; I2++) {
      OperatorData[dataCount-1][I+timeExtent][I2] = opData[I][I2];
    }
  }
}


int MassCorrelationMatrixAnalyzer::getDataCount() {
  return dataCount;
}


void MassCorrelationMatrixAnalyzer::sortData() {
  if (LogLevel>3) printf("Sorting Data in MassCorrelationMatrixAnalyzer\n");
  
printf("...not implemented yet!!!\n");  
exit(0);
}


void MassCorrelationMatrixAnalyzer::calcMassCorrelationMatrices(int ignoreStart, int ignoreEnd, bool boot) {
  if (LogLevel>3) printf("Calculating MassCorrelation-Matrix in MassCorrelationMatrixAnalyzer with ignoreStart = %d, ignoreEnd = %d and boot = %d\n",ignoreStart, ignoreEnd, boot);
  for (int I=0; I<timeExtent+1; I++) {
    MassCorrelationMatrices[I]->setZero();
  }
  if (dataCount <= 0) return;
  if ((!boot) && (ignoreStart<=0) && (ignoreEnd>=dataCount-1)) return;

  //Effective DataCount and Index-Array
  int* indices = NULL;
  int effectiveDataCount = 0;
  if (boot) {
    effectiveDataCount = dataCount;
    indices = new int[effectiveDataCount];
    for (int I=0; I<dataCount; I++) {
      int ind = (int) (dataCount * AdvancedZufall(AdvancedSeed));
      if (ind<0) ind = 0;
      if (ind >= dataCount) ind = dataCount - 1;
      indices[I] = ind;
    }  
  } else {
    for (int I=0; I<dataCount; I++) {
      if ((I<ignoreStart) || (I>ignoreEnd)) {
        effectiveDataCount++;
      }
    }
    indices = new int[effectiveDataCount];
    effectiveDataCount = 0;
    for (int I=0; I<dataCount; I++) {
      if ((I<ignoreStart) || (I>ignoreEnd)) {
        indices[effectiveDataCount] = I;
        effectiveDataCount++;
      }
    }  
  }
  
  //Average Weight
  double avgWeight = 0;
  for (int I=0; I<effectiveDataCount; I++) {
    avgWeight += dataWeights[indices[I]];
  }
  avgWeight /= effectiveDataCount;

  if (operatorCorrelationMode) {
    //Operator-Correlation-Mode
    //Correlation Matrix
    for (int I=0; I<effectiveDataCount; I++) {
      for (int t=0; t<1+timeExtent; t++) {
        for (int d1=0; d1<dimension; d1++) {
          for (int d2=0; d2<dimension; d2++) {
            Complex c0 = dataWeights[indices[I]] * OperatorCorrelationData[indices[I]][t][d1][d2];
            Complex c1 = MassCorrelationMatrices[t]->matrix[d1][d2];
            MassCorrelationMatrices[t]->matrix[d1][d2] = c1 + c0;
          }
        }
      }
    }
    
    for (int t=0; t<timeExtent+1; t++) {
      for (int d1=0; d1<dimension; d1++) {
        for (int d2=0; d2<dimension; d2++) {
          MassCorrelationMatrices[t]->matrix[d1][d2].x /= (effectiveDataCount * avgWeight);
          MassCorrelationMatrices[t]->matrix[d1][d2].y /= (effectiveDataCount * avgWeight);
        }
      }
    }
    
    //Subtract VEV of separate operators
    if (subtractVacuumExpectation) {
      Complex** avgOpDag = new Complex*[timeExtent];
      Complex** avgOp = new Complex*[timeExtent];
      for (int t=0; t<timeExtent; t++) {
        avgOpDag[t] = new Complex[dimension];
        avgOp[t] = new Complex[dimension];
        for (int d=0; d<dimension; d++) {  
          avgOpDag[t][d].x = 0;
          avgOpDag[t][d].y = 0;
          avgOp[t][d].x = 0;
          avgOp[t][d].y = 0;
        }
      }
      for (int I=0; I<effectiveDataCount; I++) {
        for (int t=0; t<timeExtent; t++) {
          for (int d=0; d<dimension; d++) {
            avgOpDag[t][d] = avgOpDag[t][d] + (dataWeights[indices[I]] * OperatorData[indices[I]][t][d]);
            avgOp[t][d] = avgOp[t][d] + (dataWeights[indices[I]] * OperatorData[indices[I]][t+timeExtent][d]);  
    	  }
        }
      }
      for (int t=0; t<timeExtent; t++) {
        for (int d=0; d<dimension; d++) {
          avgOpDag[t][d] = (1.0 / (effectiveDataCount * avgWeight)) * avgOpDag[t][d];
          avgOp[t][d] = (1.0 / (effectiveDataCount * avgWeight)) * avgOp[t][d];		
        }
      }

      for (int t1=0; t1<timeExtent; t1++) {
        for (int t2=0; t2<timeExtent; t2++) {
          for (int d1=0; d1<dimension; d1++) {
            for (int d2=0; d2<dimension; d2++) {
              int deltaT = t1-t2;
              if (deltaT<0) deltaT = -deltaT;
              int deltaT2 = timeExtent-deltaT; 
              double normFac = 2*timeExtent;
              if (t1==t2) normFac = timeExtent;

	      MassCorrelationMatrices[deltaT]->matrix[d1][d2] = MassCorrelationMatrices[deltaT]->matrix[d1][d2] - ((1.0/normFac) * (avgOpDag[t1][d1] * avgOp[t2][d2]));
	      MassCorrelationMatrices[deltaT2]->matrix[d1][d2] = MassCorrelationMatrices[deltaT2]->matrix[d1][d2] - ((1.0/normFac) * (avgOpDag[t1][d1] * avgOp[t2][d2]));
	    }
 	  }
        }
      }
    
      for (int t=0; t<timeExtent; t++) {
        delete[] avgOpDag[t];
        delete[] avgOp[t];
      }
      delete[] avgOpDag;
      delete[] avgOp;
    }
    
  } else {
    //No Operator-Correlator-Mode
    Complex** avgOp = new Complex*[timeExtent];
    for (int I=0; I<timeExtent; I++) {
      avgOp[I] = new Complex[dimension];
      for (int I2=0; I2<dimension; I2++) {
        avgOp[I][I2].x = 0;
        avgOp[I][I2].y = 0;
      }
    }
    //Average Operators
    for (int I=0; I<effectiveDataCount; I++) {
      for (int I2=0; I2<timeExtent; I2++) {
        for (int I3=0; I3<dimension; I3++) {
          avgOp[I2][I3].x += dataWeights[indices[I]] * OperatorData[indices[I]][I2][I3].x;
          avgOp[I2][I3].y += dataWeights[indices[I]] * OperatorData[indices[I]][I2][I3].y;
        }
      }
    }
    for (int I2=0; I2<timeExtent; I2++) {
      for (int I3=0; I3<dimension; I3++) {
        avgOp[I2][I3].x /= (effectiveDataCount * avgWeight);
        avgOp[I2][I3].y /= (effectiveDataCount * avgWeight);
      }
    }
    //Correlation Matrix
    for (int t1=0; t1<timeExtent; t1++) {
      for (int t2=0; t2<timeExtent; t2++) {
        int deltaT = t2-t1;
        if (deltaT<0) deltaT = timeExtent + deltaT;

        for (int d1=0; d1<dimension; d1++) {
          for (int d2=0; d2<dimension; d2++) {
            for (int I=0; I<effectiveDataCount; I++) {
              Complex c0 = dataWeights[indices[I]] * (adj(OperatorData[indices[I]][t1][d1]) * (OperatorData[indices[I]][t2][d2]));
              Complex c1 = MassCorrelationMatrices[deltaT]->matrix[d1][d2];
              MassCorrelationMatrices[deltaT]->matrix[d1][d2] = c1 + c0;
              if (deltaT==0) {
                Complex c2 = MassCorrelationMatrices[timeExtent]->matrix[d1][d2];
                MassCorrelationMatrices[timeExtent]->matrix[d1][d2] = c2 + c0;
              }
            }
          }
        }
      }
    }
    
    for (int DeltaT=0; DeltaT<timeExtent+1; DeltaT++) {
      for (int d1=0; d1<dimension; d1++) {
        for (int d2=0; d2<dimension; d2++) {
          MassCorrelationMatrices[DeltaT]->matrix[d1][d2].x /= (timeExtent*effectiveDataCount * avgWeight);
          MassCorrelationMatrices[DeltaT]->matrix[d1][d2].y /= (timeExtent*effectiveDataCount * avgWeight);

          if (subtractVacuumExpectation) {
            for (int t1=0; t1<timeExtent; t1++) {
              int t2 = t1+DeltaT;
              if (t2>=timeExtent) t2 -= timeExtent;
              MassCorrelationMatrices[DeltaT]->matrix[d1][d2] = MassCorrelationMatrices[DeltaT]->matrix[d1][d2] - ((1.0/timeExtent)*(adj(avgOp[t1][d1]) * avgOp[t2][d2]));
            }
	  }
        }
      }
    }

    for (int I=0; I<timeExtent; I++) {
      delete[] avgOp[I];
    }
    delete[] avgOp;
  }
  delete[] indices;
  
  //Bring Correlation-Matrices to generalized form if applicable and possible
  if (generalEWproblem) {
    if ((GeneralEWNormalizationMatrix == NULL) && (effectiveDataCount==dataCount)) {
      delete GeneralEWNormalizationMatrix;      
      GeneralEWNormalizationMatrix  = new ComplexMatrix(1);
      bool b = potentiateHermiteanComplexMatrix(*(MassCorrelationMatrices[0]), (*GeneralEWNormalizationMatrix), -0.5);
      bool b2 = checkMat1IsAnInverseSquareRootOfMat2((*GeneralEWNormalizationMatrix), *(MassCorrelationMatrices[0]));
      if ((!b) || (!b2)) {
        delete GeneralEWNormalizationMatrix;      
	GeneralEWNormalizationMatrix = NULL;
      }
    }
  
    if (GeneralEWNormalizationMatrix != NULL) {
      for (int DeltaT=0; DeltaT<=timeExtent; DeltaT++) {
        (*(MassCorrelationMatrices[DeltaT])) = (*GeneralEWNormalizationMatrix) * (*(MassCorrelationMatrices[DeltaT])) * (*GeneralEWNormalizationMatrix);
      }
    }
  }
}


void MassCorrelationMatrixAnalyzer::calcOptimalEVassignment(ComplexVector** v1, ComplexVector** v2, int* assignment) {
  ComplexMatrix mat(dimension);

  for (int I=0; I<dimension; I++) {
    for (int I2=0; I2<dimension; I2++) {
      mat.matrix[I][I2].x = 0;
      mat.matrix[I][I2].y = 0;
      for (int I3=0; I3<dimension; I3++) {
        mat.matrix[I][I2] = mat.matrix[I][I2] + adj(v1[I]->vectorElements[I3]) * v2[I2]->vectorElements[I3];
      }
      mat.matrix[I][I2] = adj(mat.matrix[I][I2]) * mat.matrix[I][I2];
      mat.matrix[I][I2].y = 0;
    }
  }
  for (int ewCount=0; ewCount<dimension; ewCount++) {
    double bestVal = -1;
    int bestI = -1;
    int bestI2 = -1;
    for (int I=0; I<dimension; I++) {
      for (int I2=0; I2<dimension; I2++) {
        if ((mat.matrix[I][I2].y > -1) && (mat.matrix[I][I2].x > bestVal)) {
          bestI = I;
          bestI2 = I2;
          bestVal = mat.matrix[I][I2].x;
        }
      }
    }
    assignment[bestI2] = bestI;
    for (int I=0; I<dimension; I++) {
      mat.matrix[I][bestI2].x = -2;
      mat.matrix[I][bestI2].y = -2;
      mat.matrix[bestI][I].x = -2;
      mat.matrix[bestI][I].y = -2;
    }
  }
}



bool MassCorrelationMatrixAnalyzer::calcMassCorrelationEigenvalues(int ignoreStart, int ignoreEnd, bool boot, Complex** eigenvalues, ComplexVector** sortOrder) {
  calcMassCorrelationMatrices(ignoreStart, ignoreEnd, boot);

  for (int I=0; I<timeExtent+1; I++) {
    if (!MassCorrelationMatrices[I]->calcEigenvaluesAndEigenvectors()) return false;
  }
  int* optimalAssignment = new int[dimension];
  int* optimalPermutation = new int[dimension];
  int* dummyPermutation = new int[dimension];
  for (int I=0; I<dimension; I++) {
    optimalPermutation[I] = I;
  }

  ComplexVector** vv = sortOrder;
  for (int I=0; I<(timeExtent/2)+1; I++) {
    if ((!generalEWproblem) || (I>=1)) {
      calcOptimalEVassignment(vv, MassCorrelationMatrices[I]->rightEigenVectors, optimalAssignment);

      if (iterativeEigenVectorAlignment) {
        for (int I2=0; I2<dimension; I2++) {
          dummyPermutation[I2] = optimalPermutation[optimalAssignment[I2]];
        }
        for (int I2=0; I2<dimension; I2++) {
          optimalPermutation[I2] = dummyPermutation[I2];
        }
        vv = MassCorrelationMatrices[I]->rightEigenVectors;
      } else {
        for (int I2=0; I2<dimension; I2++) {
          optimalPermutation[I2] = optimalAssignment[I2];
        }
      }
    }

    for (int I2=0; I2<dimension; I2++) {
      eigenvalues[I][optimalPermutation[I2]] = MassCorrelationMatrices[I]->eigenvalues[I2];
    }    
  }

  vv = sortOrder;
  for (int I=0; I<dimension; I++) {
    optimalPermutation[I] = I;
  }
  for (int I=timeExtent; I>=(timeExtent/2)+1; I--) {
    if ((!generalEWproblem) || (I<=timeExtent-1)) {
      calcOptimalEVassignment(vv, MassCorrelationMatrices[I]->rightEigenVectors, optimalAssignment);

      if (iterativeEigenVectorAlignment) {
        for (int I2=0; I2<dimension; I2++) {
          dummyPermutation[I2] = optimalPermutation[optimalAssignment[I2]];
        }
        for (int I2=0; I2<dimension; I2++) {
          optimalPermutation[I2] = dummyPermutation[I2];
        }
        vv = MassCorrelationMatrices[I]->rightEigenVectors;
      } else {
        for (int I2=0; I2<dimension; I2++) {
          optimalPermutation[I2] = optimalAssignment[I2];
        }    
      }
    }
    
    for (int I2=0; I2<dimension; I2++) {
      eigenvalues[I][optimalPermutation[I2]] = MassCorrelationMatrices[I]->eigenvalues[I2];
    }
  }

  delete[] optimalPermutation;
  delete[] optimalAssignment;
  delete[] dummyPermutation;
  
  return true;
}


void MassCorrelationMatrixAnalyzer::calcEigenvaluesAndErrors(bool boot, int iterations)  {
  for (int I=0; I<timeExtent+1; I++) {
    for (int I2=0; I2<dimension; I2++) {
      MassCorrelationEigenvalues[I][I2].x = 0;
      MassCorrelationEigenvalues[I][I2].y = 0;
      MassCorrelationEigenvaluesError[I][I2].x = 0;
      MassCorrelationEigenvaluesError[I][I2].y = 0;
    }
  }

  if (dataCount<=0) return;
  if (iterations<=0) return;

  calcMassCorrelationMatrices(-1,-1, false);

  int tRef = 0;
  if (generalEWproblem) tRef = 1;    
  ComplexMatrix sortMat(*MassCorrelationMatrices[tRef]);
  sortMat.calcEigenvaluesAndEigenvectors();
  Complex** eigenvalues = new Complex*[timeExtent+1];
  for (int I=0; I<timeExtent+1; I++) {
    eigenvalues[I] = new Complex[dimension];
  }
  
  int successIter = 0;
  for (int iter=0; iter<iterations; iter++) {
    int ignoreStart = (iter*dataCount) / iterations;
    int ignoreEnd = ((iter+1)*dataCount) / iterations;

    bool b = calcMassCorrelationEigenvalues(ignoreStart, ignoreEnd, boot, eigenvalues, sortMat.rightEigenVectors);
    if (b) {
      //Average Eigenvalues
      for (int I=0; I<timeExtent+1; I++) {
        for (int I2=0; I2<dimension; I2++) {
          MassCorrelationEigenvalues[I][I2].x += eigenvalues[I][I2].x;
          MassCorrelationEigenvalues[I][I2].y += eigenvalues[I][I2].y;
          MassCorrelationEigenvaluesError[I][I2].x += sqr(eigenvalues[I][I2].x);
          MassCorrelationEigenvaluesError[I][I2].y += sqr(eigenvalues[I][I2].y);
        }
      }
      successIter++;
    }
  }

  for (int I=0; I<timeExtent+1; I++) {
    delete[] eigenvalues[I];
  }
  delete[] eigenvalues;

  if (successIter==0) {
    printf("ERROR in calcEigenvaluesAndErrors: Could not determine results!\n");
    exit(0);
  }

  double errorRescale = calcErrorRescaleFactor(boot, iterations);
  for (int I=0; I<timeExtent+1; I++) {
    for (int I2=0; I2<dimension; I2++) {
      MassCorrelationEigenvalues[I][I2].x /= successIter;
      MassCorrelationEigenvalues[I][I2].y /= successIter;
      MassCorrelationEigenvaluesError[I][I2].x = errorRescale * sqrt((MassCorrelationEigenvaluesError[I][I2].x/successIter) - sqr(MassCorrelationEigenvalues[I][I2].x));
      MassCorrelationEigenvaluesError[I][I2].y = errorRescale * sqrt((MassCorrelationEigenvaluesError[I][I2].y/successIter) - sqr(MassCorrelationEigenvalues[I][I2].y));
    }
  }

  bool b = calcMassCorrelationEigenvalues(-1,-1, false, MassCorrelationEigenvalues, sortMat.rightEigenVectors);
  if (!b) {
    printf("ERROR in calcEigenvaluesAndErrors: Could not determine final results!\n");
    exit(0);
  }
}


void MassCorrelationMatrixAnalyzer::calcMassesAndErrors(bool boot, int iterations, int massCount) {
  for (int I=0; I<timeExtent; I++) {
    for (int I2=0; I2<dimension; I2++) {
      effectiveMasses[I][I2] = 0;
      effectiveMassesError[I][I2] = 0;
    }
  }
  for (int I2=0; I2<dimension; I2++) {
    asymptoticEffectiveMasses[I2] = 0;
    asymptoticEffectiveMassesError[I2] = 0;
  }
  
  for (int I=0; I<dimension; I++) {
    fittedVacuumExpectationValue[I] = 0;
    fittedVacuumExpectationValueError[I] = 0;
    delete[] fittedMasses[I];
    delete[] fittedMassesError[I];
    delete[] fittedCoefficients[I];
    delete[] fittedCoefficientsError[I];
    fitChiSquare[I] = 0;
    fitChiSquareError[I] = 0;
  }
  fittedMassCount = 0;
  
  if (dataCount <= 0) return;
  if (iterations <= 0) return;
  if (massCount <= 0) return;

  fittedMassCount = massCount;
  for (int I=0; I<dimension; I++) {
    fittedMasses[I] = new double[fittedMassCount];
    fittedMassesError[I] = new double[fittedMassCount];
    fittedCoefficients[I] = new double[fittedMassCount];
    fittedCoefficientsError[I] = new double[fittedMassCount];  
    for (int I2=0; I2<fittedMassCount; I2++) {
      fittedMasses[I][I2] = 0;
      fittedMassesError[I][I2] = 0;      
      fittedCoefficients[I][I2] = 0;
      fittedCoefficientsError[I][I2] = 0;
    }
  }
  
  double* fitX = new double[timeExtent+1];
  double* fitY = new double[timeExtent+1];
  double* fitError = new double[timeExtent+1];
  
  char** fitFunctionBody = new char*[fittedMassCount+1];
  for (int I=0; I<fittedMassCount+1; I++) fitFunctionBody[I] = new char[1000];
  int fitVarCount = 2*fittedMassCount;
  double* fitVar = new double[2*fittedMassCount+1];
  double* fitVarErrors = new double[2*fittedMassCount+1];
  snprintf(fitFunctionBody[0],1000,"A1*cosh((x-%f)*A2) ",0.5*timeExtent);
  for (int I=1; I<fittedMassCount; I++) {
    snprintf(fitFunctionBody[I],1000,"%s + A%d*cosh((x-%f)*A%d) ",fitFunctionBody[I-1],2*I+1,0.5*timeExtent,2*I+2);
  }
  if (fitVacuumExpectation) {
    snprintf(fitFunctionBody[fittedMassCount],1000,"%s + A%d",fitFunctionBody[fittedMassCount-1],2*fittedMassCount+1);
    fitVarCount++;
  }

  calcEigenvaluesAndErrors(boot, iterations);

  int tRef = 0;
  if (generalEWproblem) tRef = 1;    
  ComplexMatrix sortMat(*MassCorrelationMatrices[tRef]);
  sortMat.calcEigenvaluesAndEigenvectors();

  Complex** eigenvalues = new Complex*[timeExtent+1];
  for (int I=0; I<timeExtent+1; I++) {
    eigenvalues[I] = new Complex[dimension];
  }
  
  int* successIter = new int[dimension];
  for (int I=0; I<dimension; I++) successIter[I] = 0;
  for (int iter=0; iter<iterations; iter++) {
    int ignoreStart = (iter*dataCount) / iterations;
    int ignoreEnd = ((iter+1)*dataCount) / iterations;

    bool b = calcMassCorrelationEigenvalues(ignoreStart, ignoreEnd, boot, eigenvalues, sortMat.rightEigenVectors);
    if (b) {
      //Global Fit to cosh - functions
      for (int I2=0; I2<dimension; I2++) {
        int fitDataCount = 0;
        for (int I=0; I<timeExtent+1; I++) {
	  if ((I>=fittingRangeReduction) && (I<=timeExtent-fittingRangeReduction)) {
            fitX[fitDataCount] = I;	  
            fitY[fitDataCount] = eigenvalues[I][I2].x;
            fitError[fitDataCount] = MassCorrelationEigenvaluesError[I][I2].x;	  
	    fitDataCount++;
	  }
        }  
        double fitChi = 0;
        for (int I=0; I<fitVarCount; I++) {
	  fitVar[I] = 1E-6;
	  if ((I%2)==1) {
	    fitVar[I] = 0.05;
	    if ((fitY[0]>0) && ((fitY[1]>0))) {
	      double dummy = log(fitY[0] / fitY[1]);
	      if ((!isNaN(dummy)) && (dummy>0) && (dummy<10)) {
      	        fitVar[I] = 2*fitY[0]/(exp(-0.5*dummy*timeExtent) + exp(0.5*dummy*timeExtent));
  	        fitVar[I] = dummy;
	      }
            }
	  }
        }

        for (int I=0; I<fittedMassCount; I++) {
          b = b & performGnuplotFit(fitFunctionBody[I], fitX, fitY, fitError, fitDataCount, 2*I+2, fitVar, fitVarErrors, fitChi);
	}
        if (fitVacuumExpectation) {
          b = b & performGnuplotFit(fitFunctionBody[fittedMassCount], fitX, fitY, fitError, fitDataCount, fitVarCount, fitVar, fitVarErrors, fitChi);
	}

        if (b) {
          fitChiSquare[I2] += fitChi;
          fitChiSquareError[I2] += sqr(fitChi);
          for (int I=0; I<fittedMassCount; I++) {
            fittedCoefficients[I2][I] += fitVar[2*I+0];
            fittedCoefficientsError[I2][I] += sqr(fitVar[2*I+0]);	
            fittedMasses[I2][I] += fitVar[2*I+1];
            fittedMassesError[I2][I] += sqr(fitVar[2*I+1]);  
  	  }
          if (fitVacuumExpectation) {
            fittedVacuumExpectationValue[I2] += fitVar[2*fittedMassCount];
            fittedVacuumExpectationValueError[I2] += sqr(fitVar[2*fittedMassCount]);
  	  }
	
          //Average Effective Masses
          for (int I=0; I<timeExtent/2; I++) {
            double sub = 0;
            if (fitVacuumExpectation) {
  	      sub = fitVar[2*fittedMassCount];
	    }
            double effMass = EffectiveMassSolver(I, I+1, eigenvalues[I][I2].x-sub, eigenvalues[I+1][I2].x-sub, timeExtent);
            effectiveMasses[I][I2] += effMass;
            effectiveMassesError[I][I2] += sqr(effMass);
          }
          for (int I=timeExtent-1; I>=timeExtent/2; I--) {
            double sub = 0;
            if (fitVacuumExpectation) {
  	      sub = fitVar[2*fittedMassCount];
	    }
            double effMass = EffectiveMassSolver(I+1, I, eigenvalues[I+1][I2].x-sub, eigenvalues[I][I2].x-sub, timeExtent);
            effectiveMasses[I][I2] += effMass;
            effectiveMassesError[I][I2] += sqr(effMass);
          }
          successIter[I2]++;
	}
      }
    }
  }

  for (int I2=0; I2<dimension; I2++) {
    if (successIter[I2]==0) {
      printf("ERROR in calcMassesAndErrors: Could not determine results!\n");
      exit(0);
    }
  }

  //Calculate Averages and Errors
  double errorRescale = calcErrorRescaleFactor(boot, iterations);
  
  for (int I=0; I<timeExtent; I++) {
    for (int I2=0; I2<dimension; I2++) {
      effectiveMasses[I][I2] /= successIter[I2];
      effectiveMassesError[I][I2] = errorRescale * sqrt((effectiveMassesError[I][I2]/successIter[I2]) - sqr(effectiveMasses[I][I2]));
    }
  }
  
  for (int I2=0; I2<dimension; I2++) {
    fittedVacuumExpectationValue[I2] /= successIter[I2];
    fittedVacuumExpectationValueError[I2] = errorRescale * sqrt((fittedVacuumExpectationValueError[I2]/successIter[I2]) - sqr(fittedVacuumExpectationValue[I2]));
    for (int I=0; I<fittedMassCount; I++) {
      fittedMasses[I2][I] /= successIter[I2];
      fittedMassesError[I2][I] = errorRescale * sqrt((fittedMassesError[I2][I]/successIter[I2]) - sqr(fittedMasses[I2][I]));
      fittedCoefficients[I2][I] /= successIter[I2];
      fittedCoefficientsError[I2][I] = errorRescale * sqrt((fittedCoefficientsError[I2][I]/successIter[I2]) - sqr(fittedCoefficients[I2][I]));
    }
    fitChiSquare[I2] /= successIter[I2];
    fitChiSquareError[I2] = errorRescale * sqrt((fitChiSquareError[I2]/successIter[I2]) - sqr(fitChiSquare[I2]));
  }

  //Finaler Cosh-Fit
  for (int I2=0; I2<dimension; I2++) {
    int fitDataCount = 0;
    for (int I=0; I<timeExtent+1; I++) {
      if ((I>=fittingRangeReduction) && (I<=timeExtent-fittingRangeReduction)) {
        fitX[fitDataCount] = I;	  
        fitY[fitDataCount] = MassCorrelationEigenvalues[I][I2].x;
        fitError[fitDataCount] = MassCorrelationEigenvaluesError[I][I2].x;	  
	fitDataCount++;
      }
    }  
    for (int I=0; I<fitVarCount; I++) {
      fitVar[I] = 1E-6;
      if ((I%2)==1) {
        fitVar[I] = 0.05;
        if ((fitY[0]>0) && ((fitY[1]>0))) {
          double dummy = log(fitY[0] / fitY[1]);
	  if ((!isNaN(dummy)) && (dummy>0) && (dummy<10)) {
  	    fitVar[I] = 2*fitY[0]/(exp(-0.5*dummy*timeExtent) + exp(0.5*dummy*timeExtent));
  	    fitVar[I] = dummy;
  	  }
	}
      }
    }
    
    bool b = true;
    for (int I=0; I<fittedMassCount; I++) {
      b = b & performGnuplotFit(fitFunctionBody[I], fitX, fitY, fitError, fitDataCount, 2*I+2, fitVar, fitVarErrors, fitChiSquare[I2]);
    }
    if (fitVacuumExpectation) {
      b = b & performGnuplotFit(fitFunctionBody[fittedMassCount], fitX, fitY, fitError, fitDataCount, fitVarCount, fitVar, fitVarErrors, fitChiSquare[I2]);
    }

    if (!b) {
      printf("ERROR in calcMassesAndErrors: Could not determine final results\n");
      exit(0);    
    }

    for (int I=0; I<fittedMassCount; I++) {
      fittedCoefficients[I2][I] = fitVar[2*I+0];
      fittedMasses[I2][I] = fitVar[2*I+1];
    }
    if (fitVacuumExpectation) {
      fittedVacuumExpectationValue[I2] = fitVar[2*fittedMassCount];
    }
  }
  
  //Determine asymptotic Effective Masses
  for (int I2=0; I2<dimension; I2++) {
    asymptoticEffectiveMasses[I2] = effectiveMasses[(timeExtent/2)-1][I2] ;
    asymptoticEffectiveMassesError[I2] = effectiveMassesError[(timeExtent/2)-1][I2];  
  }


  delete[] successIter;
  delete[] fitX;
  delete[] fitY;
  delete[] fitError;
  for (int I=0; I<fittedMassCount+1; I++) delete[] fitFunctionBody[I];
  delete[] fitFunctionBody;
  delete[] fitVar;
  delete[] fitVarErrors;

  for (int I=0; I<timeExtent+1; I++) {
    delete[] eigenvalues[I];
  }
  delete[] eigenvalues;
}


double MassCorrelationMatrixAnalyzer::calcErrorRescaleFactor(bool boot, int iterations) {
  double scale = 1.0;
  
  if (!boot) {
    scale = sqrt(iterations) * (1.0 - 1.0/iterations);
  }
  
  return scale;
}


int MassCorrelationMatrixAnalyzer::calcDecentJackKnifeIterations() {
  if (dataCount <= 0) return 0;
  if (dataCount <= 10) return dataCount;
  if (dataCount <= 100) return 10;
  if (dataCount >= 10000) return 100;  
  return ((int) sqrt(dataCount));
}


void MassCorrelationMatrixAnalyzer::calcEigenvaluesWithBootStrapAnalysis(int iterations) {
  calcEigenvaluesAndErrors(true, iterations);
}


void MassCorrelationMatrixAnalyzer::calcEigenvaluesWithJackKnifeAnalysis(int iterations) {
  calcEigenvaluesAndErrors(false, iterations);
}


void MassCorrelationMatrixAnalyzer::calcEigenvaluesWithJackKnifeAnalysis() {
  int iterations = calcDecentJackKnifeIterations();
  calcEigenvaluesWithJackKnifeAnalysis(iterations);
}


void MassCorrelationMatrixAnalyzer::calcEigenvaluesAndMassesWithBootStrapAnalysis(int iterations, int massCount) {
  calcMassesAndErrors(true, iterations, massCount);
}


void MassCorrelationMatrixAnalyzer::calcEigenvaluesAndMassesWithJackKnifeAnalysis(int iterations, int massCount) {
  calcMassesAndErrors(false, iterations, massCount);
}


void MassCorrelationMatrixAnalyzer::calcEigenvaluesAndMassesWithJackKnifeAnalysis(int massCount) {
  int iterations = calcDecentJackKnifeIterations();
  calcEigenvaluesAndMassesWithJackKnifeAnalysis(iterations, massCount);  
}


void MassCorrelationMatrixAnalyzer::plotEigenvalues(bool logScale) {
  char* dataFileName = new char[1000];
  snprintf(dataFileName,1000,"data/MassCorrMatrixEW_%s.dat",uniqueAnalyzerFileNameExtension);
  FILE* file = fopen(dataFileName,"w");  
  for (int I=0; I<timeExtent+1; I++) {
    fprintf(file,"%d ",I);
    for (int I2=0; I2<dimension; I2++) {
      fprintf(file,"%1.15f %1.15f ",MassCorrelationEigenvalues[I][I2].x,MassCorrelationEigenvaluesError[I][I2].x);
    }
    fprintf(file,"\n");
  }
  fclose(file);
 
  
  char* gnuFileName = new char[1000];
  snprintf(gnuFileName,1000,"pics/MassCorrMatrixEW_%s.gnu",uniqueAnalyzerFileNameExtension);
  file = fopen(gnuFileName,"w");  
  if (logScale) {
    fprintf(file,"set logscale y\n");
  }
  fprintf(file,"set xlabel '$\\Delta t = |t_2 - t_1|$'\n");
  fprintf(file,"plot [0:%d] '%s' using 1:%d:%d with errorbars notitle\n",timeExtent,dataFileName, 2,3);
  for (int I2=1; I2<dimension; I2++) {
    fprintf(file,"replot '../pics/%s' using 1:%d:%d with errorbars notitle\n",dataFileName, 2*I2, 2*I2+1);
  }
  fprintf(file,"\n");
  for (int I2=0; I2<dimension; I2++) {
    fprintf(file,"f%d(x) = 0 ",I2);
    for (int I=0; I<fittedMassCount; I++) {
      fprintf(file,"+ %1.15f*cosh((x-%1.15f)*%1.15f) ",fittedCoefficients[I2][I],0.5*timeExtent,fittedMasses[I2][I]);
    }
    if (fitVacuumExpectation) {
      fprintf(file,"+ %1.15f ",fittedVacuumExpectationValue[I2]);
    }
    fprintf(file,"\n");
    fprintf(file,"replot f%d(x) title '$\\chi^2/dof = %1.3f$'\n",I2, fitChiSquare[I2]);
    fprintf(file,"\n");
  }
  
  
  fclose(file);
  
  delete[] gnuFileName;  
  delete[] dataFileName;  
}


void MassCorrelationMatrixAnalyzer::plotEffectiveMasses() {
  char* dataFileName = new char[1000];
  snprintf(dataFileName,1000,"data/MassCorrMatrixEffectiveMasses_%s.dat",uniqueAnalyzerFileNameExtension);
  FILE* file = fopen(dataFileName,"w");  
  for (int I=0; I<timeExtent; I++) {
    fprintf(file,"%f ",0.5 + I);
    for (int I2=0; I2<dimension; I2++) {
      fprintf(file,"%1.15f %1.15f ",effectiveMasses[I][I2],effectiveMassesError[I][I2]);
    }
    fprintf(file,"\n");
  }
  fclose(file);
 
  
  char* gnuFileName = new char[1000];
  snprintf(gnuFileName,1000,"pics/MassCorrMatrixEffectiveMasses_%s.gnu",uniqueAnalyzerFileNameExtension);
  file = fopen(gnuFileName,"w");  
  fprintf(file,"set xlabel '$\\Delta t = |t_2 - t_1|$'\n");
  fprintf(file,"set ylabel 'Effective mass'\n");  
  fprintf(file,"plot [0:%d][0:] '%s' using 1:%d:%d with errorbars notitle\n",timeExtent,dataFileName, 2,3);
  for (int I2=1; I2<dimension; I2++) {
    fprintf(file,"replot '../pics/%s' using 1:%d:%d with errorbars notitle\n",dataFileName, 2*I2, 2*I2+1);
  }
  fprintf(file,"\n");
  for (int I2=0; I2<dimension; I2++) {
    for (int I=0; I<fittedMassCount; I++) {
      fprintf(file,"replot %1.15f notitle \n",fittedMasses[I2][I]);
      fprintf(file,"\n");
    }
  }
  
  
  fclose(file);
  
  delete[] gnuFileName;  
  delete[] dataFileName; 
}


void MassCorrelationMatrixAnalyzer::plotEigenvalues(ControlLogger* logger, bool logScale) {
  double* x = new double[timeExtent+1];
  double* y = new double[timeExtent+1];
  double* err = new double[timeExtent+1];
  char* title = new char[1000];
  char* fitCommand = new char[1000];
  char* dummyStr = new char[1000];

  for (int I2=0; I2<dimension; I2++) {
    for (int I=0; I<timeExtent+1; I++) {
      x[I] = I;
      y[I] = MassCorrelationEigenvalues[I][I2].x;
      err[I] = MassCorrelationEigenvaluesError[I][I2].x;
    }

    snprintf(fitCommand,1000, "replot 0");
    for (int I=0; I<fittedMassCount; I++) {
      snprintf(dummyStr,1000,"%s + %1.15f*cosh((x-%1.15f)*%1.15f) ",fitCommand,fittedCoefficients[I2][I],0.5*timeExtent,fittedMasses[I2][I]);
      snprintf(fitCommand,1000,"%s",dummyStr);
    }    
    if (fitVacuumExpectation) {
      snprintf(dummyStr,1000,"%s + %1.15f ",fitCommand,fittedVacuumExpectationValue[I2]);
      snprintf(fitCommand,1000,"%s",dummyStr);
    }
    snprintf(dummyStr,1000,"%s title '$\\chi^2/dof = %1.2f$'",fitCommand,fitChiSquare[I2]);
    snprintf(fitCommand,1000,"%s",dummyStr);

    logger->addPlot("", NULL, "Eigenvalues of mass-correlation matrix", "", "", x, y, err, NULL, NULL, timeExtent+1, fitCommand);
  }
  
  delete[] x;
  delete[] y;
  delete[] err;
  delete[] title;
  delete[] fitCommand;
  delete[] dummyStr;
}


void MassCorrelationMatrixAnalyzer::plotEffectiveMasses(ControlLogger* logger) {
  double* x = new double[timeExtent];
  double* y = new double[timeExtent];
  double* err = new double[timeExtent];
  char* title = new char[1000];
  char* fitCommand = new char[1000];
  char* dummyStr = new char[1000];

  for (int I2=0; I2<dimension; I2++) {
    for (int I=0; I<timeExtent; I++) {
      x[I] = I;
      y[I] = effectiveMasses[I][I2];
      err[I] = effectiveMassesError[I][I2];
    }

    snprintf(fitCommand,1000,"replot %1.15f notitle \n",fittedMasses[I2][0]);
    for (int I=1; I<fittedMassCount; I++) {
      snprintf(dummyStr,1000,"%s replot %1.15f notitle \n",fitCommand,fittedMasses[I2][I]);
      snprintf(fitCommand,1000,"%s",dummyStr);
    }

    logger->addPlot("", NULL, "Effective masses from mass-correlation matrix", "", "", x, y, err, NULL, NULL, timeExtent, fitCommand);
  }
  
  delete[] x;
  delete[] y;
  delete[] err;
  delete[] title;
  delete[] fitCommand;
  delete[] dummyStr;
}


double MassCorrelationMatrixAnalyzer::getFittedMass(int operatorIndex, int massIndex) {
  if (operatorIndex<0) return NaN;
  if (operatorIndex>=dimension) return NaN;
  if (massIndex<0) return NaN;
  if (massIndex>fittedMassCount) return NaN;
  return fittedMasses[operatorIndex][massIndex];
}


double MassCorrelationMatrixAnalyzer::getFittedMassError(int operatorIndex, int massIndex) {
  if (operatorIndex<0) return NaN;
  if (operatorIndex>=dimension) return NaN;
  if (massIndex<0) return NaN;
  if (massIndex>fittedMassCount) return NaN;
  return fittedMassesError[operatorIndex][massIndex];
}


int MassCorrelationMatrixAnalyzer::getFittedMassCount() {
  return fittedMassCount;
}


int MassCorrelationMatrixAnalyzer::getDimension() {
  return dimension;
}


int MassCorrelationMatrixAnalyzer::getTimeExtent() {
  return timeExtent;
}


double MassCorrelationMatrixAnalyzer::getFittedChiSquare(int operatorIndex) {
  if (operatorIndex<0) return NaN;
  if (operatorIndex>=dimension) return NaN;
  return fitChiSquare[operatorIndex];
}


double MassCorrelationMatrixAnalyzer::getFittedChiSquareError(int operatorIndex) {
  if (operatorIndex<0) return NaN;
  if (operatorIndex>=dimension) return NaN;
  return fitChiSquareError[operatorIndex];
}


double MassCorrelationMatrixAnalyzer::getFittedCoefficient(int operatorIndex, int massIndex) {
  if (operatorIndex<0) return NaN;
  if (operatorIndex>=dimension) return NaN;
  if (massIndex<0) return NaN;
  if (massIndex>=fittedMassCount) return NaN;
  return fittedCoefficients[operatorIndex][massIndex];
}


double MassCorrelationMatrixAnalyzer::getFittedCoefficientError(int operatorIndex, int massIndex) {
  if (operatorIndex<0) return NaN;
  if (operatorIndex>=dimension) return NaN;
  if (massIndex<0) return NaN;
  if (massIndex>=fittedMassCount) return NaN;
  return fittedCoefficientsError[operatorIndex][massIndex];
}


double MassCorrelationMatrixAnalyzer::getFittedVacuumExpectationValue(int operatorIndex) {
  if (operatorIndex<0) return NaN;
  if (operatorIndex>=dimension) return NaN;
  return fittedVacuumExpectationValue[operatorIndex];
}


double MassCorrelationMatrixAnalyzer::getFittedVacuumExpectationValueError(int operatorIndex) {
  if (operatorIndex<0) return NaN;
  if (operatorIndex>=dimension) return NaN;
  return fittedVacuumExpectationValueError[operatorIndex];
}


double MassCorrelationMatrixAnalyzer::getFittedEffectiveMasses(int timeIndex, int operatorIndex) {
  if (operatorIndex<0) return NaN;
  if (operatorIndex>=dimension) return NaN;
  if (timeIndex<0) return NaN;
  if (timeIndex>=timeExtent) return NaN;  
  return effectiveMasses[timeIndex][operatorIndex];
}


double MassCorrelationMatrixAnalyzer::getFittedEffectiveMassesError(int timeIndex, int operatorIndex) {
  if (operatorIndex<0) return NaN;
  if (operatorIndex>=dimension) return NaN;
  if (timeIndex<0) return NaN;
  if (timeIndex>=timeExtent) return NaN;  
  return effectiveMassesError[timeIndex][operatorIndex];
}


double MassCorrelationMatrixAnalyzer::getFittedAsymptoticEffectiveMasses(int operatorIndex) {
  if (operatorIndex<0) return NaN;
  if (operatorIndex>=dimension) return NaN;
  return asymptoticEffectiveMasses[operatorIndex];
}


double MassCorrelationMatrixAnalyzer::getFittedAsymptoticEffectiveMassesError(int operatorIndex) {
  if (operatorIndex<0) return NaN;
  if (operatorIndex>=dimension) return NaN;
  return asymptoticEffectiveMassesError[operatorIndex];
}


double MassCorrelationMatrixAnalyzer::getMassCorrelationEigenvalue(int timeIndex, int operatorIndex) {
  if (operatorIndex<0) return NaN;
  if (operatorIndex>=dimension) return NaN;
  if (timeIndex<0) return NaN;
  if (timeIndex>timeExtent) return NaN;  
  return MassCorrelationEigenvalues[timeIndex][operatorIndex].x;
}


double MassCorrelationMatrixAnalyzer::getMassCorrelationEigenvalueError(int timeIndex, int operatorIndex) {
  if (operatorIndex<0) return NaN;
  if (operatorIndex>=dimension) return NaN;
  if (timeIndex<0) return NaN;
  if (timeIndex>timeExtent) return NaN;  
  return MassCorrelationEigenvaluesError[timeIndex][operatorIndex].x;
}


void MassCorrelationMatrixAnalyzer::plotEigenvalues() {
calcMassesAndErrors(false, 10, 1);

  char* fileName = new char[1000];
  snprintf(fileName,1000,"EWcorr%s.dat",uniqueAnalyzerFileNameExtension);
  FILE* file = fopen(fileName,"w");  
  for (int I=0; I<timeExtent+1; I++) {
    printf("DeltaT = %d\n",I);
    fprintf(file,"%d ",I);
    for (int I2=0; I2<dimension; I2++) {
      fprintf(file,"%1.15f %1.15f ",MassCorrelationEigenvalues[I][I2].x,MassCorrelationEigenvaluesError[I][I2].x);
      printf("%e+-%e %e+-%e\n",MassCorrelationEigenvalues[I][I2].x,MassCorrelationEigenvaluesError[I][I2].x,MassCorrelationEigenvalues[I][I2].y,MassCorrelationEigenvaluesError[I][I2].y);
    }
    fprintf(file,"\n");

  }
  fclose(file);

  snprintf(fileName,1000,"EffectiveMasses%s.dat",uniqueAnalyzerFileNameExtension);
  file = fopen(fileName,"w");  
  for (int I=0; I<timeExtent; I++) {
    printf("DeltaT = %1.1f\n",I+0.5);
    fprintf(file,"%f ",I+0.5);
    for (int I2=0; I2<dimension; I2++) {
      if ((!isNaN(effectiveMasses[I][I2])) && (!isNaN(effectiveMassesError[I][I2]))) {
        fprintf(file,"%1.15f %1.15f ",effectiveMasses[I][I2],effectiveMassesError[I][I2]);
      } else {
        fprintf(file,"-1 -1 ");
      }
      printf("%e+-%e\n",effectiveMasses[I][I2],effectiveMassesError[I][I2]);
    }
    fprintf(file,"\n");
  }
  fclose(file);

  snprintf(fileName,1000,"FitParameter%s.dat",uniqueAnalyzerFileNameExtension);
  file = fopen(fileName,"w");  
  for (int I2=0; I2<dimension; I2++) {
    for (int I=0; I<fittedMassCount; I++) {
      fprintf(file,"%1.15f %1.15f %1.15f %1.15f ", fittedMasses[I2][I], fittedMassesError[I2][I], fittedCoefficients[I2][I], fittedCoefficientsError[I2][I]);
    }
    if (fitVacuumExpectation) {
      fprintf(file,"%1.15f %1.15f ",fittedVacuumExpectationValue[I2],fittedVacuumExpectationValueError[I2]);
    }
    fprintf(file,"%1.15f %1.15f",fitChiSquare[I2],fitChiSquareError[I2]);
    
    fprintf(file,"\n");
  }
  fclose(file);
  
  
  delete[] fileName;
}
