#include "AnalyzerObservableScalarCondensate.h"

AnalyzerObservableScalarCondensate::AnalyzerObservableScalarCondensate(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "ScalarCondensate", "scond") { 
  ini(getAnalyzerResultsCount());
}


AnalyzerObservableScalarCondensate::~AnalyzerObservableScalarCondensate() {
}


bool AnalyzerObservableScalarCondensate::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  double* my_phi_field = phiFieldConf->getPhiFieldCopy();
  int my_L0 = SDReader->getL0();
  int my_L1 = SDReader->getL1();
  int my_L2 = SDReader->getL2();
  int my_L3 = SDReader->getL3();

  int var_ctr = 0;
  int x_neighbour_index = 0;
  int y_neighbour_index = 0;
  int z_neighbour_index = 0;
  int t_neighbour_index = 0;

  analyzerResults[0] = 0;
  for (int L0_ctr=0; L0_ctr<my_L0; L0_ctr++)
    {
      for (int L1_ctr=0; L1_ctr<my_L1; L1_ctr++)
	{
	  for (int L2_ctr=0; L2_ctr<my_L2; L2_ctr++)
	    {
	      for (int L3_ctr=0; L3_ctr<my_L3; L3_ctr++)
		{
		  for (int cmp_ctr=0; cmp_ctr<4; cmp_ctr++)
		    {

		      t_neighbour_index = cmp_ctr + ((L3_ctr+1) % my_L3)*4 + L2_ctr*4*my_L3 + L1_ctr*4*my_L3*my_L2 + L0_ctr*4*my_L3*my_L2*my_L1;
		      z_neighbour_index = cmp_ctr + L3_ctr*4 + ((L2_ctr+1) % my_L2)*4*my_L3 + L1_ctr*4*my_L3*my_L2 + L0_ctr*4*my_L3*my_L2*my_L1;
		      y_neighbour_index = cmp_ctr + L3_ctr*4 + L2_ctr*4*my_L3 + ((L1_ctr+1) % my_L1)*4*my_L3*my_L2 + L0_ctr*4*my_L3*my_L2*my_L1;
		      x_neighbour_index = cmp_ctr + L3_ctr*4 + L2_ctr*4*my_L3 + L1_ctr*4*my_L3*my_L2 + ((L0_ctr+1) % my_L0)*4*my_L3*my_L2*my_L1;

                      analyzerResults[0] += my_phi_field[var_ctr]*my_phi_field[x_neighbour_index]
                                          + my_phi_field[var_ctr]*my_phi_field[y_neighbour_index]
                                          + my_phi_field[var_ctr]*my_phi_field[z_neighbour_index]
                                          + my_phi_field[var_ctr]*my_phi_field[t_neighbour_index];
		      var_ctr = var_ctr + 1;
		    }
		}
	    }
	}
    }
  analyzerResults[0] /= (my_L0*my_L1*my_L2*my_L3);

  return true;
}


int AnalyzerObservableScalarCondensate::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservableScalarCondensate::getAnalyzerResultsCount() {
  return 1;
}
