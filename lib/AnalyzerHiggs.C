#include "AnalyzerHiggs.h"
#include "Tools.h"

void AnalyzerHiggs::ini(int l0, int l1, int l2, int l3, ControlLogger* log, bool FLAG_GFit) {
  totalN = 0;
  L0 = l0;
  L1 = l1;
  L2 = l2;
  L3 = l3;  
  FLAG_GnuplotFit = FLAG_GFit;
  LatticeResult_VEV = NaN;
  LatticeResult_VEVsigma = NaN;
  LatticeResult_HiggsPropagatorMass = NaN;
  LatticeResult_HiggsPropagatorMassSigma = NaN;
  LatticeResult_GoldstoneZFactor = NaN;
  LatticeResult_GoldstoneZFactorSigma = NaN;
  LatticeResult_PhysicalHiggsMass = NaN;
  LatticeResult_PhysicalHiggsMassSigma = NaN;

  timeDirection = 0;
  LargestL = L0;
  spaceDirection1 = 1;
  spaceDirection2 = 2;
  spaceDirection3 = 3;
  if (L1>LargestL) {
    LargestL = L1;
    timeDirection = 1;
    spaceDirection1 = 0;
    spaceDirection2 = 2;
    spaceDirection3 = 3;
  }
  if (L2>LargestL) {
    LargestL = L2;
    timeDirection = 2;
    spaceDirection1 = 0;
    spaceDirection2 = 1;
    spaceDirection3 = 3;
  }
  if (L3>LargestL) {
    LargestL = L3;
    timeDirection = 3;
    spaceDirection1 = 0;
    spaceDirection2 = 1;
    spaceDirection3 = 2;
  }
  
  MassControlLog = log;
  weightData = new double[AnalyzerHiggsDataMAX];
  phiFieldBuffer = new Complex[4*L0*L1*L2*L3];
  phiMomentumBuffer = new Complex[4*L0*L1*L2*L3];
  sinPSqr = new double[L0*L1*L2*L3];
  HiggsPropagator = new double[L0*L1*L2*L3];
  GoldstonePropagator = new double[L0*L1*L2*L3];
  HiggsTimeSliceData = new double*[AnalyzerHiggsDataMAX];
  GoldstoneTimeSliceData = new double*[AnalyzerHiggsDataMAX];
  HiggsTimeSliceCorrelatorData = new double*[AnalyzerHiggsDataMAX];
  GoldstoneTimeSliceCorrelatorData = new double*[AnalyzerHiggsDataMAX];
  HiggsTimeSliceCorrelator = new AutoCorrelation*[1+LargestL];
  GoldstoneTimeSliceCorrelator = new AutoCorrelation*[1+LargestL];
  HiggsVEVdata = new double[AnalyzerHiggsDataMAX];
  HiggsVEV = new AutoCorrelation(2,100);

  MomentumSqrSlotLocations = new double[AnalyzerHiggsDataMAX];
  MomentumSqrSlotCount = 0;
  HiggsPropagatorSlottedData = new double*[AnalyzerHiggsDataMAX];
  GoldstonePropagatorSlottedData = new double*[AnalyzerHiggsDataMAX];

  GoldstoneTwoParticleMassAnalyzer = new MassCorrelationMatrixAnalyzer(8, LargestL, false, false, true, true, false,  0, 100000,"Goldstone2Particle");
  HiggsGoldstoneMassAnalyzer = new MassCorrelationMatrixAnalyzer(4, LargestL, false, false, true, true, true, 0, 100000,"HiggsGoldstone");
  HiggsMassAnalyzer = new MassCorrelationMatrixAnalyzer(1, LargestL, false, false, true, true, true, 0, 100000,"Higgs");

  
  int I;
  for (I=0; I<1+LargestL; I++) {
    HiggsTimeSliceCorrelator[I] = new AutoCorrelation(2,100);
    GoldstoneTimeSliceCorrelator[I] = new AutoCorrelation(2,100);
  }

  LatticeResult_PhysicalEffectiveHiggsMasses = new double[LargestL/2-1];
  LatticeResult_PhysicalEffectiveHiggsMassesSigmas = new double[LargestL/2-1];
  for (I=0; I<LargestL/2-1; I++) {
    LatticeResult_PhysicalEffectiveHiggsMasses[I] = NaN;
    LatticeResult_PhysicalEffectiveHiggsMassesSigmas[I] = NaN;    
  }
  
  
  for (I=0; I<L0*L1*L2*L3; I++) {
    HiggsPropagator[I] = 0;
    GoldstonePropagator[I] = 0;   
  }
  
  int I0,I1,I2,I3;
  vector4D p;
  int count = 0;
  for (I0=0; I0<L0; I0++) {
    p[0] = 2*I0*pi/L0;
    for (I1=0; I1<L1; I1++) {
      p[1] = 2*I1*pi/L1;
      for (I2=0; I2<L2; I2++) {
        p[2] = 2*I2*pi/L2;
        for (I3=0; I3<L3; I3++) {
          p[3] = 2*I3*pi/L3;
	  
	  sinPSqr[count] = 4.0 * (sqr(sin(0.5*p[0])) + sqr(sin(0.5*p[1])) + sqr(sin(0.5*p[2])) + sqr(sin(0.5*p[3])));
	  
	  int slotNr = findMomentumSqrSlot(sinPSqr[count]);
	  if (slotNr<0) {
	    MomentumSqrSlotLocations[MomentumSqrSlotCount] = sinPSqr[count];
	    MomentumSqrSlotCount++;
	  }
	  count++;
	}
      }
    }
  }
  
  HiggsPropagatorSlottedDataCorrelation = new AutoCorrelation*[MomentumSqrSlotCount];
  GoldstonePropagatorSlottedDataCorrelation = new AutoCorrelation*[MomentumSqrSlotCount];  
  for (I=0; I<MomentumSqrSlotCount; I++) {
    HiggsPropagatorSlottedDataCorrelation[I] = new AutoCorrelation(2,100);
    GoldstonePropagatorSlottedDataCorrelation[I] = new AutoCorrelation(2,100);  
  }
  
    
  int* n = new int[4];
  n[0] = L0;
  n[1] = L1;
  n[2] = L2;
  n[3] = L3;
  
  int rank = 4;
  int howmany = 4;
  int* inembed = new int[4];
  inembed[0] = L0;
  inembed[1] = L1;
  inembed[2] = L2;
  inembed[3] = L3;  
  int istride = 4;
  int idist = 1;
  
  int* onembed = new int[4];
  onembed[0] = L0;
  onembed[1] = L1;
  onembed[2] = L2;
  onembed[3] = L3;    
  int ostride = 4;
  int odist = 1;

  phiFieldFourierForwardPlan4Components = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*) phiFieldBuffer, inembed,
                                                istride, idist,
	   			                (fftw_complex*)phiMomentumBuffer, onembed, ostride, odist,
	  				        FFTW_FORWARD, FFTW_MEASURE);
						
  howmany = 1;
  
  phiFieldFourierForwardPlan1Components = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*) phiFieldBuffer, inembed,
                                                istride, idist,
	   			                (fftw_complex*)phiMomentumBuffer, onembed, ostride, odist,
	  				        FFTW_FORWARD, FFTW_MEASURE);
    						
  delete[] n;
  delete[] inembed;
  delete[] onembed;
}


void AnalyzerHiggs::desini() {
  delete[] phiFieldBuffer;
  delete[] phiMomentumBuffer;
  delete[] sinPSqr;
  delete[] HiggsPropagator;
  delete[] GoldstonePropagator;
  fftw_destroy_plan(phiFieldFourierForwardPlan1Components);
  fftw_destroy_plan(phiFieldFourierForwardPlan4Components);
  int I;
  for (I=0; I<1+LargestL/2; I++) {
    delete HiggsTimeSliceCorrelator[I];
    delete GoldstoneTimeSliceCorrelator[I];
  }
  for (I=0; I<totalN; I++) {
    delete[] HiggsTimeSliceCorrelatorData[I];
    delete[] GoldstoneTimeSliceCorrelatorData[I];
    delete[] HiggsTimeSliceData[I];
    delete[] GoldstoneTimeSliceData[I];
    delete[] HiggsPropagatorSlottedData[I];
    delete[] GoldstonePropagatorSlottedData[I];
  }  
  delete[] LatticeResult_PhysicalEffectiveHiggsMasses;
  delete[] LatticeResult_PhysicalEffectiveHiggsMassesSigmas;
  delete[] HiggsPropagatorSlottedData;
  delete[] GoldstonePropagatorSlottedData;
  delete[] MomentumSqrSlotLocations;
  for (I=0; I<MomentumSqrSlotCount; I++) {
    delete HiggsPropagatorSlottedDataCorrelation[I];
    delete GoldstonePropagatorSlottedDataCorrelation[I];
  }
  delete[] HiggsVEVdata;
  delete[] HiggsTimeSliceData;
  delete[] GoldstoneTimeSliceData;
  delete[] HiggsTimeSliceCorrelator;
  delete[] GoldstoneTimeSliceCorrelator;
  delete[] HiggsTimeSliceCorrelatorData;
  delete[] GoldstoneTimeSliceCorrelatorData;
  delete[] weightData;
  delete HiggsVEV;
  delete GoldstoneTwoParticleMassAnalyzer;
  delete HiggsGoldstoneMassAnalyzer;
  delete HiggsMassAnalyzer;
}


AnalyzerHiggs::AnalyzerHiggs(int l0, int l1, int l2, int l3, ControlLogger* log, bool FLAG_GFit) {
  ini(l0,l1,l2,l3, log, FLAG_GFit);
}


AnalyzerHiggs::~AnalyzerHiggs() {
  desini();
}


int AnalyzerHiggs::getTotalN() {
  return totalN;
}


int AnalyzerHiggs::findMomentumSqrSlot(double momSqr) {
  int I;
  for(I=0; I<MomentumSqrSlotCount; I++) {
    if (fabs(MomentumSqrSlotLocations[I]-momSqr) < AnalyzerHiggsMomentumSqrSlotSize) {
      return I;    
    }  
  }
  return -1;
}


void AnalyzerHiggs::measureMags(vector4D* phiField) {
  int I1,I2,I3,I4;
  int count = 0;
  vector4D measurePhi, measureStaggeredPhi;
  double avgNorm, sigmaNorm;
  double dummy;
  
  measurePhi[0] = 0;
  measurePhi[1] = 0;
  measurePhi[2] = 0;
  measurePhi[3] = 0;
  measureStaggeredPhi[0] = 0;
  measureStaggeredPhi[1] = 0;
  measureStaggeredPhi[2] = 0;
  measureStaggeredPhi[3] = 0;
  avgNorm = 0;
  sigmaNorm = 0;
  
  count = 0;
  for (I1=0; I1<L0; I1++) {
    for (I2=0; I2<L1; I2++) {
      for (I3=0; I3<L2; I3++) {
        for (I4=0; I4<L3; I4++) {
          measurePhi[0] += phiField[count][0];
          measurePhi[1] += phiField[count][1];
          measurePhi[2] += phiField[count][2];
          measurePhi[3] += phiField[count][3];
	  
	  avgNorm += dummy = sqrt( phiField[count][0]*phiField[count][0] 
	                          +phiField[count][1]*phiField[count][1]
	                          +phiField[count][2]*phiField[count][2]
	                          +phiField[count][3]*phiField[count][3]);
          sigmaNorm += dummy*dummy;
	  
	  double stagFac = 1;
	  if (((I1+I2+I3+I4) % 2) == 1) stagFac = -1;
          measureStaggeredPhi[0] += stagFac*phiField[count][0];
          measureStaggeredPhi[1] += stagFac*phiField[count][1];
          measureStaggeredPhi[2] += stagFac*phiField[count][2];
          measureStaggeredPhi[3] += stagFac*phiField[count][3];
	  count++;
	}
      }
    }
  }
  double norm = L0*L1*L2*L3;
  measurePhi[0] = measurePhi[0]/norm;
  measurePhi[1] = measurePhi[1]/norm;
  measurePhi[2] = measurePhi[2]/norm;
  measurePhi[3] = measurePhi[3]/norm;
  measureStaggeredPhi[0] = measureStaggeredPhi[0]/norm;
  measureStaggeredPhi[1] = measureStaggeredPhi[1]/norm;
  measureStaggeredPhi[2] = measureStaggeredPhi[2]/norm;
  measureStaggeredPhi[3] = measureStaggeredPhi[3]/norm;
  
  avgNorm /= norm;
  sigmaNorm = sqrt(sigmaNorm/norm - avgNorm*avgNorm);
  
  measurePhiNorm = measurePhi[0]*measurePhi[0]
                 + measurePhi[1]*measurePhi[1]
                 + measurePhi[2]*measurePhi[2]
                 + measurePhi[3]*measurePhi[3];
  measureStaggeredPhiNorm = measureStaggeredPhi[0]*measureStaggeredPhi[0]
                          + measureStaggeredPhi[1]*measureStaggeredPhi[1]
                          + measureStaggeredPhi[2]*measureStaggeredPhi[2]
                          + measureStaggeredPhi[3]*measureStaggeredPhi[3];
		 
  measurePhiNorm = sqrt(measurePhiNorm);
  measureStaggeredPhiNorm = sqrt(measureStaggeredPhiNorm);
  
  if (LogLevel>4) {
    printf("m = %1.15f\n", measurePhiNorm);
    printf("s = %1.15f\n", measureStaggeredPhiNorm);
  }
  if (LogLevel>4) printf("\n");
}


void AnalyzerHiggs::analyzeHiggsField(vector4D* phiField, double weight) {
  if (LogLevel>2) printf("Analyzing Higgs Field with weight=%f and time direction = %d...\n",weight,timeDirection);

  weightData[totalN] = weight;
  measureMags(phiField);
  HiggsVEVdata[totalN] = measurePhiNorm;

  vector4D averageVector;
  averageVector[0] = 0;
  averageVector[1] = 0;
  averageVector[2] = 0;
  averageVector[3] = 0;
  
  int I;
  for (I=0; I<L0*L1*L2*L3; I++) {
    averageVector[0] += phiField[I][0];
    averageVector[1] += phiField[I][1];
    averageVector[2] += phiField[I][2];
    averageVector[3] += phiField[I][3];
  }  
  double norm = sqr(averageVector[0]);
  norm += sqr(averageVector[1]);
  norm += sqr(averageVector[2]);
  norm += sqr(averageVector[3]);
  norm = sqrt(norm);
  averageVector[0] /= norm;
  averageVector[1] /= norm;
  averageVector[2] /= norm;
  averageVector[3] /= norm;
  
  
  //Higgs-Goldstone-Analysator
  int count = 0;
  Complex** higgsGoldstoneData = new Complex*[LargestL];
  for (I=0; I<LargestL; I++) {
    higgsGoldstoneData[I] = new Complex[4];
    for (int I2=0; I2<4; I2++) {
      higgsGoldstoneData[I][I2].x = 0;
      higgsGoldstoneData[I][I2].y = 0;
    }
  }
  int indy[4];
  for (indy[0]=0; indy[0]<L0; indy[0]++) {
    for (indy[1]=0; indy[1]<L1; indy[1]++) {
      for (indy[2]=0; indy[2]<L2; indy[2]++) {
        for (indy[3]=0; indy[3]<L3; indy[3]++) {
	  higgsGoldstoneData[indy[timeDirection]][0].x += phiField[count][0];
	  higgsGoldstoneData[indy[timeDirection]][1].x += phiField[count][1];
	  higgsGoldstoneData[indy[timeDirection]][2].x += phiField[count][2];
	  higgsGoldstoneData[indy[timeDirection]][3].x += phiField[count][3];	
	  count++;
	}
      }
    }
  }
  for (I=0; I<LargestL; I++) {
    for (int I2=0; I2<4; I2++) {
      higgsGoldstoneData[I][I2].x /= L0*L1*L2*L3;
      higgsGoldstoneData[I][I2].y /= L0*L1*L2*L3;
      higgsGoldstoneData[I][I2].x *= LargestL;
      higgsGoldstoneData[I][I2].y *= LargestL;
    }
  }
  
  HiggsGoldstoneMassAnalyzer->addOperatorData(0, weight, higgsGoldstoneData);
  for (I=0; I<LargestL; I++) {
    delete[] higgsGoldstoneData[I];
  }
  delete[] higgsGoldstoneData;
  
  
  //Higgs-Mode: Propagator
  count = 0;
  for (I=0; I<L0*L1*L2*L3; I++) {
    double scalar = phiField[I][0]*averageVector[0];
    scalar += phiField[I][1]*averageVector[1];
    scalar += phiField[I][2]*averageVector[2];
    scalar += phiField[I][3]*averageVector[3];
    phiFieldBuffer[count].x = scalar;
    phiFieldBuffer[count].y = 0;
    phiFieldBuffer[count+1].x = NaN;
    phiFieldBuffer[count+1].y = NaN;
    phiFieldBuffer[count+2].x = NaN;
    phiFieldBuffer[count+2].y = NaN;
    phiFieldBuffer[count+3].x = NaN;
    phiFieldBuffer[count+3].y = NaN;
    count += 4;
  }
    
  fftw_execute(phiFieldFourierForwardPlan1Components);
 
  count = 0;
  double normFac = 1.0 / ((L0*L1*L2*L3));
  double* tempDataHiggsProp = new double[MomentumSqrSlotCount];
  int* tempDataHiggsPropCount = new int[MomentumSqrSlotCount];
  for (I=0; I<MomentumSqrSlotCount; I++) {
    tempDataHiggsProp[I] = 0;
    tempDataHiggsPropCount[I] = 0;
  }
  
  for (I=0; I<L0*L1*L2*L3; I++) {
    phiMomentumBuffer[count].x = normFac * (sqr(phiMomentumBuffer[count].x) + sqr(phiMomentumBuffer[count].y));
    HiggsPropagator[I] += weightData[totalN] * phiMomentumBuffer[count].x;
    int slotNr = findMomentumSqrSlot(sinPSqr[I]);
    tempDataHiggsProp[slotNr] += phiMomentumBuffer[count].x;
    tempDataHiggsPropCount[slotNr]++;    
    phiMomentumBuffer[count].y = 0;
    count += 4;
  }
  
  HiggsPropagatorSlottedData[totalN] = new double[MomentumSqrSlotCount];
  for (I=0; I<MomentumSqrSlotCount; I++) {
    HiggsPropagatorSlottedData[totalN][I] = tempDataHiggsProp[I] / tempDataHiggsPropCount[I];
  }
  delete[] tempDataHiggsProp;
  delete[] tempDataHiggsPropCount;
  
  //Higgs-Mode: Time-Slice Data
  HiggsTimeSliceData[totalN] = new double[LargestL];
  for (I=0; I<LargestL; I++) HiggsTimeSliceData[totalN][I] = 0;
  int ind[4];
  count = 0;
  for (ind[0]=0; ind[0]<L0; ind[0]++) {
    for (ind[1]=0; ind[1]<L1; ind[1]++) {
      for (ind[2]=0; ind[2]<L2; ind[2]++) {
        for (ind[3]=0; ind[3]<L3; ind[3]++) {
          HiggsTimeSliceData[totalN][ind[timeDirection]] += phiFieldBuffer[count].x;	
          count += 4;
	}
      }
    }
  }
  for (I=0; I<LargestL; I++) {
    HiggsTimeSliceData[totalN][I] /= L0*L1*L2*L3;
    HiggsTimeSliceData[totalN][I] *= LargestL;
  }
  
  
  //Higgs-Mode: Time-Slice Data - Unitaere Eichung
  vector4D* avgVec = new vector4D[LargestL];
  for (I=0; I<LargestL; I++) {
    avgVec[I][0] = 0;
    avgVec[I][1] = 0;
    avgVec[I][2] = 0;
    avgVec[I][3] = 0;    
  }
  count = 0;
  for (ind[0]=0; ind[0]<L0; ind[0]++) {
    for (ind[1]=0; ind[1]<L1; ind[1]++) {
      for (ind[2]=0; ind[2]<L2; ind[2]++) {
        for (ind[3]=0; ind[3]<L3; ind[3]++) {
          avgVec[ind[timeDirection]][0] += phiField[count][0];
          avgVec[ind[timeDirection]][1] += phiField[count][1];
          avgVec[ind[timeDirection]][2] += phiField[count][2];
          avgVec[ind[timeDirection]][3] += phiField[count][3];
          count++;
	}
      }
    }
  }
  
  for (I=0; I<LargestL; I++) {
    double norm = sqrt(sqr(avgVec[I][0]) + sqr(avgVec[I][1]) + sqr(avgVec[I][2]) + sqr(avgVec[I][3]));
    avgVec[I][0] /= norm;
    avgVec[I][1] /= norm;
    avgVec[I][2] /= norm;
    avgVec[I][3] /= norm;
  }
  
  
  Complex** higgsData = new Complex*[LargestL];
  for (int t=0; t<LargestL; t++) {
    higgsData[t] = new Complex[1];
    higgsData[t][0].x = 0;
    higgsData[t][0].y = 0;
  }
  count = 0;
  for (ind[0]=0; ind[0]<L0; ind[0]++) {
    for (ind[1]=0; ind[1]<L1; ind[1]++) {
      for (ind[2]=0; ind[2]<L2; ind[2]++) {
        for (ind[3]=0; ind[3]<L3; ind[3]++) {
//          higgsData[ind[timeDirection]][0].x += sqrt(sqr(phiField[count][0]) + sqr(phiField[count][1]) + sqr(phiField[count][2]) + sqr(phiField[count][3]));
//          count++;
          higgsData[ind[timeDirection]][0].x += phiFieldBuffer[count].x;	
          count+=4;
//          higgsData[ind[timeDirection]][0].x += avgVec[ind[timeDirection]][0]*phiField[count][0] + avgVec[ind[timeDirection]][1]*phiField[count][1] + avgVec[ind[timeDirection]][2]*phiField[count][2] + avgVec[ind[timeDirection]][3]*phiField[count][3];
//          count++;
	}
      }
    }
  }
  for (I=0; I<LargestL; I++) {
    higgsData[I][0].x /= L0*L1*L2*L3;
    higgsData[I][0].x *= LargestL;
  }
  
  HiggsMassAnalyzer->addOperatorData(0, weight, higgsData);
  
  for (int t=0; t<LargestL; t++) {
    delete[] higgsData[t];
  }
  delete[] higgsData;
  
  
  //Goldstone-Modes
  count = 0;
  for (ind[0]=0; ind[0]<L0; ind[0]++) {
    for (ind[1]=0; ind[1]<L1; ind[1]++) {
      for (ind[2]=0; ind[2]<L2; ind[2]++) {
        for (ind[3]=0; ind[3]<L3; ind[3]++) {
          double scalar = 0;
          scalar += phiField[count][0]*averageVector[0];
          scalar += phiField[count][1]*averageVector[1];
          scalar += phiField[count][2]*averageVector[2];
          scalar += phiField[count][3]*averageVector[3];
    
    
/*          scalar = 0;
          scalar += phiField[count][0]*avgVec[ind[timeDirection]][0];
          scalar += phiField[count][1]*avgVec[ind[timeDirection]][1];
          scalar += phiField[count][2]*avgVec[ind[timeDirection]][2];
          scalar += phiField[count][3]*avgVec[ind[timeDirection]][3];*/
    
    
          phiFieldBuffer[4*count+0].x = phiField[count][0] - scalar*averageVector[0];
          phiFieldBuffer[4*count+0].y = 0;
          phiFieldBuffer[4*count+1].x = phiField[count][1] - scalar*averageVector[1];
          phiFieldBuffer[4*count+1].y = 0;
          phiFieldBuffer[4*count+2].x = phiField[count][2] - scalar*averageVector[2];
          phiFieldBuffer[4*count+2].y = 0;
          phiFieldBuffer[4*count+3].x = phiField[count][3] - scalar*averageVector[3];
          phiFieldBuffer[4*count+3].y = 0;

/*          phiFieldBuffer[4*count+0].x = phiField[count][0] - scalar*avgVec[ind[timeDirection]][0];
          phiFieldBuffer[4*count+0].y = 0;
          phiFieldBuffer[4*count+1].x = phiField[count][1] - scalar*avgVec[ind[timeDirection]][1];
          phiFieldBuffer[4*count+1].y = 0;
          phiFieldBuffer[4*count+2].x = phiField[count][2] - scalar*avgVec[ind[timeDirection]][2];
          phiFieldBuffer[4*count+2].y = 0;
          phiFieldBuffer[4*count+3].x = phiField[count][3] - scalar*avgVec[ind[timeDirection]][3];
          phiFieldBuffer[4*count+3].y = 0;*/
          count++;
	}
      }
    }
  }
  
  fftw_execute(phiFieldFourierForwardPlan4Components);
  
  
  
  //Goldstone-2-Particle-Time-Slice-Correlator from Fourier-Transformation
  int relImp[4];
  int relImpMax[4];
  relImpMax[0] = 2;
  relImpMax[1] = 2;
  relImpMax[2] = 2;
  relImpMax[3] = 2;
  relImpMax[timeDirection] = 1;

  Complex** opData = new Complex*[LargestL];
  for (int t=0; t<LargestL; t++) {
    opData[t] = new Complex[relImpMax[0]*relImpMax[1]*relImpMax[2]*relImpMax[3]];
  }
  int relImpCount = 0;
  for (relImp[0]=0; relImp[0]<relImpMax[0]; relImp[0]++) {
    for (relImp[1]=0; relImp[1]<relImpMax[1]; relImp[1]++) {
      for (relImp[2]=0; relImp[2]<relImpMax[2]; relImp[2]++) {
        for (relImp[3]=0; relImp[3]<relImpMax[3]; relImp[3]++) {
          for (int t=0; t<LargestL; t++) {
            opData[t][relImpCount].x = 0;
            opData[t][relImpCount].y = 0;
	  
	  
  	    for (int pt1=0; pt1<LargestL; pt1++) {
  	      for (int pt2=0; pt2<LargestL; pt2++) {
                relImp[timeDirection] = pt1;
	        int pos1 = (L3+relImp[3]) % L3;
	        pos1 += L3*relImp[2];
	        pos1 += L2*L3*relImp[1];
	        pos1 += L1*L2*L3*relImp[0];
		pos1 *= 4;

                relImp[timeDirection] = pt2;
	        int pos2 = (L3+relImp[3]) % L3;
	        pos2 += L3*relImp[2];
	        pos2 += L2*L3*relImp[1];
	        pos2 += L1*L2*L3*relImp[0];
		pos2 *= 4;
	      
                double expArg = (2*pi*t*(pt2-pt1))/LargestL;
                Complex fac = exp(expArg*ComplexI);
	      
                Complex scalar = adj(phiMomentumBuffer[pos1+0])*phiMomentumBuffer[pos2+0]
                               + adj(phiMomentumBuffer[pos1+1])*phiMomentumBuffer[pos2+1] 
                               + adj(phiMomentumBuffer[pos1+2])*phiMomentumBuffer[pos2+2] 
                               + adj(phiMomentumBuffer[pos1+3])*phiMomentumBuffer[pos2+3];

                opData[t][relImpCount] = opData[t][relImpCount] + (scalar * fac);
                relImp[timeDirection] = 0;
	      }
	    }
	  }

          relImpCount++;
        }
      }
    }
  }
  for (int I=0; I<LargestL; I++) {
    for (int I2=0; I2<relImpCount; I2++) {
      opData[I][I2].x /= L0*L1*L2*L3;
      opData[I][I2].x /= L0*L1*L2*L3;
      opData[I][I2].y /= L0*L1*L2*L3;
      opData[I][I2].y /= L0*L1*L2*L3;
    }
  }

  GoldstoneTwoParticleMassAnalyzer->addOperatorData(0, weight, opData);
  printf("\n");

  for (int I=0; I<LargestL; I++) {
    delete[] opData[I];
  }
  delete[] opData;
  delete[] avgVec;
  
  
  
  //Goldstone - Propagator  
  count = 0;
  normFac = 1.0 / (3.0*(L0*L1*L2*L3));
  double* tempDataGoldstoneProp = new double[MomentumSqrSlotCount];
  int* tempDataGoldstonePropCount = new int[MomentumSqrSlotCount];
  for (I=0; I<MomentumSqrSlotCount; I++) {
    tempDataGoldstoneProp[I] = 0;
    tempDataGoldstonePropCount[I] = 0;
  }
  
  for (I=0; I<L0*L1*L2*L3; I++) {
    phiMomentumBuffer[count].x  = sqr(phiMomentumBuffer[count+0].x) + sqr(phiMomentumBuffer[count+0].y);
    phiMomentumBuffer[count].x += sqr(phiMomentumBuffer[count+1].x) + sqr(phiMomentumBuffer[count+1].y);
    phiMomentumBuffer[count].x += sqr(phiMomentumBuffer[count+2].x) + sqr(phiMomentumBuffer[count+2].y);
    phiMomentumBuffer[count].x += sqr(phiMomentumBuffer[count+3].x) + sqr(phiMomentumBuffer[count+3].y);
    phiMomentumBuffer[count].x *= normFac,
    GoldstonePropagator[I] += weightData[totalN] * phiMomentumBuffer[count].x;
    int slotNr = findMomentumSqrSlot(sinPSqr[I]);
    tempDataGoldstoneProp[slotNr] += phiMomentumBuffer[count].x;
    tempDataGoldstonePropCount[slotNr]++;    
    
    phiMomentumBuffer[count].y = 0;
    count += 4;
  }
  
  GoldstonePropagatorSlottedData[totalN] = new double[MomentumSqrSlotCount];
  for (I=0; I<MomentumSqrSlotCount; I++) {
    GoldstonePropagatorSlottedData[totalN][I] = tempDataGoldstoneProp[I] / tempDataGoldstonePropCount[I];
  }
  delete[] tempDataGoldstoneProp;
  delete[] tempDataGoldstonePropCount;      
  
  
  //Goldstone-Mode: Time-Slice Data
  GoldstoneTimeSliceData[totalN] = new double[4*LargestL];
  for (I=0; I<4*LargestL; I++) GoldstoneTimeSliceData[totalN][I] = 0;
  count = 0;
  for (ind[0]=0; ind[0]<L0; ind[0]++) {
    for (ind[1]=0; ind[1]<L1; ind[1]++) {
      for (ind[2]=0; ind[2]<L2; ind[2]++) {
        for (ind[3]=0; ind[3]<L3; ind[3]++) {
          GoldstoneTimeSliceData[totalN][4*ind[timeDirection]+0] += phiFieldBuffer[count+0].x;	
          GoldstoneTimeSliceData[totalN][4*ind[timeDirection]+1] += phiFieldBuffer[count+1].x;	
          GoldstoneTimeSliceData[totalN][4*ind[timeDirection]+2] += phiFieldBuffer[count+2].x;	
          GoldstoneTimeSliceData[totalN][4*ind[timeDirection]+3] += phiFieldBuffer[count+3].x;	
	  
          count += 4;
	}
      }
    }
  }
  for (I=0; I<4*LargestL; I++) {
    GoldstoneTimeSliceData[totalN][I] /= L0*L1*L2*L3;
    GoldstoneTimeSliceData[totalN][I] *= LargestL;
  }




  totalN++;  
}


void AnalyzerHiggs::calcHiggsVEV() {
  int* RunLengths = new int[1];
  RunLengths[0] = totalN;
  HiggsVEV->loadData(1, RunLengths, weightData, HiggsVEVdata);
  LatticeResult_VEV = HiggsVEV->getAverage(1);

  ComplexVector derivatives(5);
  derivatives.setZero();
  derivatives.vectorElements[1].x = 1;
  LatticeResult_VEVsigma = HiggsVEV->estimateCombinedError(derivatives);
  
  delete[] RunLengths;
}
  
  
void AnalyzerHiggs::calcHiggsTimeSliceCorrelator() {
  int I,I2,I3;
  double vev = HiggsVEV->getAverage(1);
  printf("vev = %f\n",vev);
  for (I3=0; I3<totalN; I3++) {
    HiggsTimeSliceCorrelatorData[I3] = new double[1+LargestL];
    for (I=0; I<1+LargestL; I++) HiggsTimeSliceCorrelatorData[I3][I] = 0;
    for (I=0; I<LargestL; I++) {
      for (I2=0; I2<LargestL; I2++) {
        int deltaT = I2-I;
        if (deltaT<0) deltaT = -deltaT;
        int deltaT2 = LargestL-deltaT; 
	
        HiggsTimeSliceCorrelatorData[I3][deltaT] += (HiggsTimeSliceData[I3][I] - vev) * (HiggsTimeSliceData[I3][I2] - vev);
        HiggsTimeSliceCorrelatorData[I3][deltaT2] += (HiggsTimeSliceData[I3][I] - vev) * (HiggsTimeSliceData[I3][I2] - vev);
      }   
    }
    for (I=0; I<1+LargestL; I++) {
      int normFac = 2*LargestL;
      if ((I==0) || (I==LargestL)) normFac = LargestL;
      HiggsTimeSliceCorrelatorData[I3][I] /= normFac;
    }
  }
}


void AnalyzerHiggs::calcGoldstoneTimeSliceCorrelator() {
  int I,I2,I3;

  for (I3=0; I3<totalN; I3++) {
    GoldstoneTimeSliceCorrelatorData[I3] = new double[1+LargestL];
    for (I=0; I<1+LargestL; I++) GoldstoneTimeSliceCorrelatorData[I3][I] = 0;
    for (I=0; I<LargestL; I++) {
      for (I2=0; I2<LargestL; I2++) {
        int deltaT = I2-I;
        if (deltaT<0) deltaT = -deltaT;
        int deltaT2 = LargestL-deltaT; 
		
	
        GoldstoneTimeSliceCorrelatorData[I3][deltaT] += (GoldstoneTimeSliceData[I3][4*I+0] - 0) * (GoldstoneTimeSliceData[I3][4*I2+0] - 0);
        GoldstoneTimeSliceCorrelatorData[I3][deltaT] += (GoldstoneTimeSliceData[I3][4*I+1] - 0) * (GoldstoneTimeSliceData[I3][4*I2+1] - 0);
        GoldstoneTimeSliceCorrelatorData[I3][deltaT] += (GoldstoneTimeSliceData[I3][4*I+2] - 0) * (GoldstoneTimeSliceData[I3][4*I2+2] - 0);
        GoldstoneTimeSliceCorrelatorData[I3][deltaT] += (GoldstoneTimeSliceData[I3][4*I+3] - 0) * (GoldstoneTimeSliceData[I3][4*I2+3] - 0);
        GoldstoneTimeSliceCorrelatorData[I3][deltaT2] += (GoldstoneTimeSliceData[I3][4*I+0] - 0) * (GoldstoneTimeSliceData[I3][4*I2+0] - 0);
        GoldstoneTimeSliceCorrelatorData[I3][deltaT2] += (GoldstoneTimeSliceData[I3][4*I+1] - 0) * (GoldstoneTimeSliceData[I3][4*I2+1] - 0);
        GoldstoneTimeSliceCorrelatorData[I3][deltaT2] += (GoldstoneTimeSliceData[I3][4*I+2] - 0) * (GoldstoneTimeSliceData[I3][4*I2+2] - 0);
        GoldstoneTimeSliceCorrelatorData[I3][deltaT2] += (GoldstoneTimeSliceData[I3][4*I+3] - 0) * (GoldstoneTimeSliceData[I3][4*I2+3] - 0);
      }   
    }
    for (I=0; I<1+LargestL; I++) {
      int normFac = 2*LargestL;
      if ((I==0) || (I==LargestL)) normFac = LargestL;
      GoldstoneTimeSliceCorrelatorData[I3][I] /= normFac;
    }
  }
}


double AnalyzerHiggs::getLastMagnetization() {
  return measurePhiNorm;
}


double AnalyzerHiggs::getLastStaggeredMagnetization() {
  return measureStaggeredPhiNorm;
}


void AnalyzerHiggs::plotHiggsPropagator() {
  if (LogLevel>2) printf("Plotting Higgs-Propagator...\n");
  FILE* file = fopen("data/InvPropHiggs.dat","w");
  double* x = new double[L0*L1*L2*L3];
  double* y = new double[L0*L1*L2*L3];
  double* err = new double[L0*L1*L2*L3];
  int I;
  double avgWeight = 0;
  for (I=0; I<totalN; I++) avgWeight += weightData[I];
  avgWeight /= totalN;
  for (I=0; I<L0*L1*L2*L3; I++) {
    double v = HiggsPropagator[I] / (totalN * avgWeight);
    v = 1.0/v;
    x[I] = sinPSqr[I];
    y[I] = v;
    err[I] = 0.01;
    
    fprintf(file,"%1.15f %1.15f\n",sinPSqr[I],v);
  }
  fclose(file);
  
  char* fitCommand = new char[1000];
  if (FLAG_GnuplotFit) {
    char* functionBody = new char[1000];
    double* fitRes = new double[2];
    double* fitErr = new double[2];
    double redChiSqr = NaN;
  
    snprintf(functionBody,1000,"A1+A2*x");
    fitRes[0] = 1.0;
    fitRes[1] = 1.0;
    performGnuplotFit(functionBody, &(x[1]), &(y[1]), &(err[1]), L0*L1*L2*L3-1, 2, fitRes, fitErr, redChiSqr);
    snprintf(fitCommand,1000,"replot %1.15f+%1.15f*x notitle",fitRes[0],fitRes[1]);
    LatticeResult_HiggsPropagatorMass = sqrt(fitRes[0]);
    printf("Determined Higgs Propagator Mass: %f +- %f\n",sqrt(fitRes[0]),fitErr[0]);
    
    delete[] functionBody;
    delete[] fitRes;
    delete[] fitErr;
  } else {
    snprintf(fitCommand,1000,"\n");
    LatticeResult_HiggsPropagatorMass = NaN;  
  }
  MassControlLog->addSection("Higgs - Propagator");
  MassControlLog->addPlot("", NULL, "Higgs Propagator", "$\\hat p^2$", "$G_{\\phi}^{-1}(\\hat p^2)$", &(x[1]), &(y[1]), &(err[1]), NULL, NULL, L0*L1*L2*L3-1, fitCommand);
  
  delete[] fitCommand;
  delete[] x;
  delete[] y;
  delete[] err;
}


void AnalyzerHiggs::plotSlottedHiggsPropagator() {
  int reps = (int) (0.5*sqrt(totalN));
  if (reps>totalN) reps = totalN;
  int blockSize = totalN / reps;
  if (blockSize<1) blockSize = 1;
  int I;

  double avg = 0;
  double sigma = 0;
 
  for (I=0; I<reps; I++) {
    int igStart = I * blockSize;
    int igEnd = (I+1) * blockSize;
    if (igEnd>=totalN) igEnd = totalN-1;
    plotSlottedHiggsPropagator(igStart, igEnd);

    avg += LatticeResult_HiggsPropagatorMass;
    sigma += sqr(LatticeResult_HiggsPropagatorMass);
  }

  plotSlottedHiggsPropagator(-1, -1);
  LatticeResult_HiggsPropagatorMassSigma = sqrt(sigma/reps - sqr(avg/reps));
}


void AnalyzerHiggs::plotSlottedHiggsPropagator(int ignoreStart, int ignoreEnd) {
  if (LogLevel>2) printf("Plotting slotted Higgs-Propagator...\n");
  double* tempData = new double[totalN];
  double* tempWeights = new double[totalN];
  int I,I2;
  int* RunLengths = new int[1];
  for (I=0; I<MomentumSqrSlotCount; I++) {
    RunLengths[0] = 0;
    for (I2=0; I2<totalN; I2++) {
      if ((I2<ignoreStart) || (I2>ignoreEnd)) {
        tempData[RunLengths[0]] = HiggsPropagatorSlottedData[I2][I];
        tempWeights[RunLengths[0]] = weightData[I2];
        RunLengths[0]++;
      }
    }
    HiggsPropagatorSlottedDataCorrelation[I]->loadData(1, RunLengths, tempWeights, tempData);    
  }
  delete[] tempData;  
  delete[] tempWeights;
  delete[] RunLengths;
  
  
  FILE* file = fopen("data/InvPropHiggsSlotted.dat","w");
  double* x = new double[MomentumSqrSlotCount];
  double* y = new double[MomentumSqrSlotCount];
  double* err = new double[MomentumSqrSlotCount];
  double* autoCorr = new double[MomentumSqrSlotCount];  
  ComplexVector derivatives(5);
  derivatives.setZero();
  derivatives.vectorElements[1].x = 1;
  for (I=0; I<MomentumSqrSlotCount; I++) {
    double v = HiggsPropagatorSlottedDataCorrelation[I]->getAverage(1);
    x[I] = MomentumSqrSlotLocations[I];
    y[I] = 1/v;
    err[I] = HiggsPropagatorSlottedDataCorrelation[I]->estimateCombinedError(derivatives)/(v*v);
    autoCorr[I] = HiggsPropagatorSlottedDataCorrelation[I]->estimateAutoCorrelationTime();
    
    fprintf(file,"%1.15f %1.15f %1.15f %1.15f\n",x[I],y[I],err[I],autoCorr[I]);
  }
  fclose(file);
  
  char* fitCommand = new char[1000];
  if (FLAG_GnuplotFit) {
    char* functionBody = new char[1000];
    double* fitRes = new double[2];
    double* fitErr = new double[2];
    double redChiSqr = NaN;
  
    snprintf(functionBody,1000,"A1+A2*x");
    fitRes[0] = 1.0;
    fitRes[1] = 1.0;
    performGnuplotFit(functionBody, &(x[1]), &(y[1]), &(err[1]), MomentumSqrSlotCount-1, 2, fitRes, fitErr, redChiSqr);
    snprintf(fitCommand,1000,"replot %1.15f+%1.15f*x notitle",fitRes[0],fitRes[1]);
    LatticeResult_HiggsPropagatorMass = sqrt(fitRes[0]);
    printf("Determined (slotted) Higgs Propagator Mass: %f +- %f\n",sqrt(fitRes[0]),fitErr[0]);
    
    delete[] functionBody;
    delete[] fitRes;
    delete[] fitErr;
  } else {
    snprintf(fitCommand,1000,"\n");
    LatticeResult_HiggsPropagatorMass = NaN;  
  }
  if ((ignoreStart<0) && (ignoreEnd<0)) {
    MassControlLog->addSection("Higgs - Propagator (slotted)");
    MassControlLog->addPlot("", NULL, "Higgs Propagator (slotted)", "$\\hat p^2$", "$G_{\\phi}^{-1}(\\hat p^2)$", &(x[1]), &(y[1]), &(err[1]), NULL, NULL, MomentumSqrSlotCount-1, fitCommand);
  }

  delete[] fitCommand;
  delete[] x;
  delete[] y;
  delete[] err;
  delete[] autoCorr;  
}


void AnalyzerHiggs::plotGoldstonePropagator() {
  if (LogLevel>2) printf("Plotting Goldstone-Propagator...\n");
  FILE* file = fopen("data/InvPropGoldstone.dat","w");
  double* x = new double[L0*L1*L2*L3];
  double* y = new double[L0*L1*L2*L3];
  double* err = new double[L0*L1*L2*L3];
  int I;
  double avgWeight = 0;
  for (I=0; I<totalN; I++) avgWeight += weightData[I];
  avgWeight /= totalN;
  for (I=0; I<L0*L1*L2*L3; I++) {
    double v = GoldstonePropagator[I] / (totalN * avgWeight);
    v = 1.0/v;
    x[I] = sinPSqr[I];
    y[I] = v;
    err[I] = 0.01;
    
    fprintf(file,"%1.15f %1.15f\n",sinPSqr[I],v);
  }
  fclose(file);
  
  char* fitCommand = new char[1000];
  if (FLAG_GnuplotFit) {
    char* functionBody = new char[1000];
    double* fitRes = new double[2];
    double* fitErr = new double[2];
    double redChiSqr = NaN;
  
    snprintf(functionBody,1000,"A1+A2*x");
    fitRes[0] = 1.0;
    fitRes[1] = 1.0;
    performGnuplotFit(functionBody, &(x[1]), &(y[1]), &(err[1]), L0*L1*L2*L3-1, 2, fitRes, fitErr, redChiSqr);
    LatticeResult_GoldstoneZFactor = 1.0 / fitRes[1];
    printf("Determined Goldstone Z-Faktor: %f +- %f\n",1.0 / fitRes[1],fitErr[1]);
  
    snprintf(fitCommand,1000,"replot %1.15f+%1.15f*x notitle",fitRes[0],fitRes[1]);
    
    delete[] functionBody;
    delete[] fitRes;
    delete[] fitErr;    
  } else {
     LatticeResult_GoldstoneZFactor = NaN,
     snprintf(fitCommand,1000,"\n"); 
  }

  MassControlLog->addSection("Goldstone - Propagator");
  MassControlLog->addPlot("", NULL, "Goldstone Propagator", "$\\hat p^2$", "$G_{\\pi}^{-1}(\\hat p^2)$", &(x[1]), &(y[1]), &(err[1]), NULL, NULL, L0*L1*L2*L3-1, fitCommand);
   
  delete[] fitCommand;
  delete[] x;
  delete[] y;
  delete[] err;  
}


void AnalyzerHiggs::plotSlottedGoldstonePropagator() {
  int reps = (int) (0.5*sqrt(totalN));
  if (reps>totalN) reps = totalN;
  int blockSize = totalN / reps;
  if (blockSize<1) blockSize = 1;
  int I;

  double avg = 0;
  double sigma = 0;

  for (I=0; I<reps; I++) {
    int igStart = I * blockSize;
    int igEnd = (I+1) * blockSize;
    if (igEnd>=totalN) igEnd = totalN-1;
    plotSlottedGoldstonePropagator(igStart, igEnd);

    avg += LatticeResult_GoldstoneZFactor;
    sigma += sqr(LatticeResult_GoldstoneZFactor);
  }

  plotSlottedGoldstonePropagator(-1, -1);
  LatticeResult_GoldstoneZFactorSigma = sqrt(sigma/reps - sqr(avg/reps));
}


void AnalyzerHiggs::plotSlottedGoldstonePropagator(int ignoreStart, int ignoreEnd) {
  if (LogLevel>2) printf("Plotting slotted Goldstone-Propagator...\n");
  double* tempData = new double[totalN];
  double* tempWeights = new double[totalN];
  int I,I2;
  int* RunLengths = new int[1];
  for (I=0; I<MomentumSqrSlotCount; I++) {
    RunLengths[0] = 0;
    for (I2=0; I2<totalN; I2++) {
      if ((I2<ignoreStart) || (I2>ignoreEnd)) {
        tempData[RunLengths[0]] = GoldstonePropagatorSlottedData[I2][I];
        tempWeights[RunLengths[0]] = weightData[I2];
        RunLengths[0]++;
      }
    }
    GoldstonePropagatorSlottedDataCorrelation[I]->loadData(1, RunLengths, tempWeights, tempData);    
  }
  delete[] tempData;  
  delete[] tempWeights;
  delete[] RunLengths;
  
  
  FILE* file = fopen("data/InvPropGoldstoneSlotted.dat","w");
  double* x = new double[MomentumSqrSlotCount];
  double* y = new double[MomentumSqrSlotCount];
  double* err = new double[MomentumSqrSlotCount];
  double* autoCorr = new double[MomentumSqrSlotCount];  
  ComplexVector derivatives(5);
  derivatives.setZero();
  derivatives.vectorElements[1].x = 1;
  for (I=0; I<MomentumSqrSlotCount; I++) {
    double v = GoldstonePropagatorSlottedDataCorrelation[I]->getAverage(1);
    x[I] = MomentumSqrSlotLocations[I];
    y[I] = 1/v;
    err[I] = GoldstonePropagatorSlottedDataCorrelation[I]->estimateCombinedError(derivatives)/(v*v);
    autoCorr[I] = GoldstonePropagatorSlottedDataCorrelation[I]->estimateAutoCorrelationTime();
    
    fprintf(file,"%1.15f %1.15f %1.15f %1.15f\n",x[I],y[I],err[I],autoCorr[I]);
  }
  fclose(file);
  
  char* fitCommand = new char[1000];
  if (FLAG_GnuplotFit) {
    char* functionBody = new char[1000];
    double* fitRes = new double[2];
    double* fitErr = new double[2];
    double redChiSqr = NaN;
  
    snprintf(functionBody,1000,"A1+A2*x");
    fitRes[0] = 1.0;
    fitRes[1] = 1.0;
    performGnuplotFit(functionBody, &(x[1]), &(y[1]), &(err[1]), MomentumSqrSlotCount-1, 2, fitRes, fitErr, redChiSqr);
    LatticeResult_GoldstoneZFactor = 1.0 / fitRes[1];
    printf("Determined Goldstone Z-Faktor (slotted): %f +- %f\n",1.0 / fitRes[1],fitErr[1]);
  
    snprintf(fitCommand,1000,"replot %1.15f+%1.15f*x notitle",fitRes[0],fitRes[1]);
    
    delete[] functionBody;
    delete[] fitRes;
    delete[] fitErr;    
  } else {
     LatticeResult_GoldstoneZFactor = NaN,
     snprintf(fitCommand,1000,"\n"); 
  }

  if ((ignoreStart<0) && (ignoreEnd<0)) {
    MassControlLog->addSection("Goldstone - Propagator (slotted)");
    MassControlLog->addPlot("", NULL, "Goldstone Propagator (slotted)", "$\\hat p^2$", "$G_{\\pi}^{-1}(\\hat p^2)$", &(x[1]), &(y[1]), &(err[1]), NULL, NULL, MomentumSqrSlotCount-1, fitCommand);
  }

  delete[] fitCommand;
  delete[] x;
  delete[] y;
  delete[] err;
  delete[] autoCorr;  
}


void AnalyzerHiggs::plotHiggsTimeSliceCorrelator() {
  int reps = (int) (0.5*sqrt(totalN));
  if (reps>totalN) reps = totalN;
  int blockSize = totalN / reps;
  if (blockSize<1) blockSize = 1;
  int I,I2;

  double avg = 0;
  double sigma = 0;
  double* effectiveMassesAvg = new double[LargestL/2+1];
  double* effectiveMassesSigma = new double[LargestL/2+1];
  int* repCount = new int[LargestL/2+1];
  
  for (I2=0; I2<LargestL/2-1; I2++) {
    effectiveMassesAvg[I2] = 0;
    effectiveMassesSigma[I2] = 0; 
    repCount[I2] = 0;
  }

  for (I=0; I<reps; I++) {
    int igStart = I * blockSize;
    int igEnd = (I+1) * blockSize;
    if (igEnd>=totalN) igEnd = totalN-1;
    plotHiggsTimeSliceCorrelator(igStart, igEnd);

    avg += LatticeResult_PhysicalHiggsMass;
    sigma += sqr(LatticeResult_PhysicalHiggsMass);
    for (I2=0; I2<LargestL/2-1; I2++) {
      if (LatticeResult_PhysicalEffectiveHiggsMasses[I2] == LatticeResult_PhysicalEffectiveHiggsMasses[I2]) {
        //i.e. not NaN
        effectiveMassesAvg[I2] += LatticeResult_PhysicalEffectiveHiggsMasses[I2];
        effectiveMassesSigma[I2] += sqr(LatticeResult_PhysicalEffectiveHiggsMasses[I2]);
        repCount[I2]++;
      }
    }
  }

  plotHiggsTimeSliceCorrelator(-1, -1);
  LatticeResult_PhysicalHiggsMassSigma = sqrt(sigma/reps - sqr(avg/reps));
  double* x = new double[LargestL/2+1];
  double* y = new double[LargestL/2+1];
  double* err = new double[LargestL/2+1];
  int effMassCount = 0;
  for (I2=0; I2<LargestL/2-1; I2++) {
    if (repCount[I2]>0) {
      LatticeResult_PhysicalEffectiveHiggsMasses[I2] = effectiveMassesAvg[I2] / repCount[I2];
      LatticeResult_PhysicalEffectiveHiggsMassesSigmas[I2] = sqrt(effectiveMassesSigma[I2] / repCount[I2] - sqr(effectiveMassesAvg[I2] / repCount[I2]));
      x[effMassCount] = I2;
      y[effMassCount] = LatticeResult_PhysicalEffectiveHiggsMasses[I2];
      err[effMassCount] = LatticeResult_PhysicalEffectiveHiggsMassesSigmas[I2];      
      effMassCount++;
    } else {
      LatticeResult_PhysicalEffectiveHiggsMasses[I2] = NaN;
      LatticeResult_PhysicalEffectiveHiggsMassesSigmas[I2] = NaN;    
    }  
  }
  x[effMassCount] = LargestL/2;
  y[effMassCount] = 0;
  err[effMassCount] = 0;
  effMassCount++; 
  
  MassControlLog->addSection("Physical Effective Higgs Masses");
  MassControlLog->addPlot("", NULL, "Effective Masses", "$\\Delta t=|t_2-t_1|$", "$m_{eff}(\\Delta t)$", x, y, err, NULL, NULL, effMassCount, "");
    
  delete[] x,
  delete[] y;
  delete[] err;
  delete[] effectiveMassesAvg;
  delete[] effectiveMassesSigma;
  delete[] repCount;
}


void AnalyzerHiggs::plotHiggsTimeSliceCorrelator(int ignoreStart, int ignoreEnd) {
  if (LogLevel>2) printf("Plotting Higgs-Time-Slice Correlator...\n");
  FILE* file = fopen("data/TimeSliceCorrHiggs.dat","w");
  int I,I2;
  double* x = new double[1+LargestL];
  double* y = new double[1+LargestL];
  double* err = new double[1+LargestL];
  int* RunLengths = new int[1];
  double* dummyData = new double[totalN];
  double* tempWeights = new double[totalN];
  ComplexVector derivatives(5);
  derivatives.setZero();
  derivatives.vectorElements[1].x = 1;
  for (I=0; I<1+LargestL; I++) {
    RunLengths[0] = 0;
    for (I2=0; I2<totalN; I2++) {
      if ((I2<ignoreStart) || (I2>ignoreEnd)) {
        dummyData[RunLengths[0]] = HiggsTimeSliceCorrelatorData[I2][I];
        tempWeights[RunLengths[0]] = weightData[I2];
        RunLengths[0]++;
      }
    }
    HiggsTimeSliceCorrelator[I]->loadData(1, RunLengths, tempWeights, dummyData);
  
    x[I] = I;
    y[I] = HiggsTimeSliceCorrelator[I]->getAverage(1);
    err[I] = HiggsTimeSliceCorrelator[I]->estimateCombinedError(derivatives);
    double autocorr = HiggsTimeSliceCorrelator[I]->estimateAutoCorrelationTime();
    printf("Auto Correlation Time: %f\n",autocorr);
    fprintf(file,"%1.15f %1.15f %1.15f\n",x[I],y[I], err[I]);
  }
  fclose(file);
  delete[] RunLengths;
  delete[] tempWeights;
  delete[] dummyData;  
  
  double zoomFac = 1.0 / y[1];
  for (I=0; I<1+LargestL; I++) {
    y[I] *= zoomFac;
    err[I] *= zoomFac;
  }  
  
  char* fitCommand = new char[1000];
  if (FLAG_GnuplotFit) {
    char* functionBody = new char[1000];
    double* fitRes = new double[2];
    double* fitErr = new double[2];
    double redChiSqr = NaN;
  

    snprintf(functionBody,1000,"A1*cosh(A2*(x-%d))",LargestL/2);
    double v1 = y[1];
    double v2 = y[2];
    
    fitRes[1] = log(v1/v2);
    fitRes[0] = v1 / cosh(fitRes[1]*(LargestL/2 - 1));
    performGnuplotFit(functionBody, &(x[1]), &(y[1]), &(err[1]), LargestL-2, 2, fitRes, fitErr, redChiSqr);

    LatticeResult_PhysicalHiggsMass = fitRes[1];
    printf("Determined Physical Higgs-Mass: %f +- %f\n",fitRes[1],fitErr[1]);
    snprintf(fitCommand,1000,"replot %1.15f*cosh(%1.15f*(x-%d)) notitle",fitRes[0]/zoomFac,fitRes[1],LargestL/2);

    double save0 = fitRes[0] ;       
    double save1 = fitRes[1] ; 
    for (I=0; I<LargestL/2-1; I++) {
      if ((y[I]>0) && (y[I+1]>0)) {
        performGnuplotFit(functionBody, &(x[I]), &(y[I]), &(err[I]), 3, 2, fitRes, fitErr, redChiSqr);
        LatticeResult_PhysicalEffectiveHiggsMasses[I] = fitRes[1];
        LatticeResult_PhysicalEffectiveHiggsMassesSigmas[I] = fitErr[1];	
      } else {
        LatticeResult_PhysicalEffectiveHiggsMasses[I] = NaN;
        LatticeResult_PhysicalEffectiveHiggsMassesSigmas[I] = NaN;          
      }
      fitRes[0] = save0;
      fitRes[1] = save1;      
    }
      
    delete[] functionBody;
    delete[] fitRes;
    delete[] fitErr;    
  } else {
    LatticeResult_PhysicalHiggsMass = NaN;
    snprintf(fitCommand,1000,"\n"); 
    for (I=0; I<LargestL/2-1; I++) {
      LatticeResult_PhysicalEffectiveHiggsMasses[I] = NaN;
      LatticeResult_PhysicalEffectiveHiggsMassesSigmas[I] = NaN;    
    }
  }

  for (I=0; I<1+LargestL; I++) {
    y[I] /= zoomFac;
    err[I] /= zoomFac;
  }  

  if ((ignoreStart<0) && (ignoreEnd<0)) {
    MassControlLog->addSection("Physical Higgs Mass");
    MassControlLog->addPlot("", NULL, "Time slice correlator", "$\\Delta t=|t_2-t_1|$", "$\\langle\\Phi_{t_1}\\Phi_{t_2}\\rangle$", x, y, err, NULL, NULL, 1+LargestL, fitCommand);
  }
   
  delete[] fitCommand;
  delete[] x;
  delete[] y;
  delete[] err;  
}


void AnalyzerHiggs::plotGoldstoneTimeSliceCorrelator() {
  if (LogLevel>2) printf("Plotting Goldstone-Time-Slice Correlator...\n");
  FILE* file = fopen("data/TimeSliceCorrGoldstone.dat","w");
  int I,I2;
  int* RunLengths = new int[1];
  RunLengths[0] = totalN;
  double* dummyData = new double[totalN];
  ComplexVector derivatives(5);
  derivatives.setZero();
  derivatives.vectorElements[1].x = 1;
  for (I=0; I<1+LargestL; I++) {
    for (I2=0; I2<totalN; I2++) dummyData[I2] = GoldstoneTimeSliceCorrelatorData[I2][I];
    GoldstoneTimeSliceCorrelator[I]->loadData(1, RunLengths, weightData, dummyData);
  
    double v = GoldstoneTimeSliceCorrelator[I]->getAverage(1);
    double sigma = GoldstoneTimeSliceCorrelator[I]->estimateCombinedError(derivatives);
    double autocorr = GoldstoneTimeSliceCorrelator[I]->estimateAutoCorrelationTime();
    printf("Auto Correlation Time: %f\n",autocorr);
    fprintf(file,"%d %1.15f %1.15f\n",I,v, sigma);
  }
  fclose(file);
  delete[] RunLengths;
  delete[] dummyData;  
}


void AnalyzerHiggs::plotGoldstone2ParticleMasses() {
  GoldstoneTwoParticleMassAnalyzer->plotEigenvalues();

}


void AnalyzerHiggs::plotHiggsGoldstoneMasses() {
  HiggsGoldstoneMassAnalyzer->plotEigenvalues();

}


void AnalyzerHiggs::plotHiggsMasses() {
  HiggsMassAnalyzer->plotEigenvalues();
  HiggsMassAnalyzer->plotEigenvalues(true);
  HiggsMassAnalyzer->plotEffectiveMasses();
  HiggsMassAnalyzer->plotEigenvalues(MassControlLog, true);
  HiggsMassAnalyzer->plotEffectiveMasses(MassControlLog);
  
  LatticeResult_PhysicalHiggsMass = HiggsMassAnalyzer->getFittedMass(0,0);
  LatticeResult_PhysicalHiggsMassSigma = HiggsMassAnalyzer->getFittedMassError(0,0);
  
}
