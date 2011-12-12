
AM_CXXFLAGS = -I@top_srcdir@/lib 
AM_LDFLAGS = -L@top_builddir@/lib -L@top_builddir@/asmCode64
LDADD = -lhyphmc @LIBS@

bin_PROGRAMS = pHMC analyzer evaluateConfs PerturbativePredictions \
							 HiggsMassFromLargeNf PhaseDiagramAtLargeLambda

pHMC_SOURCES= pHMC.C

analyzer_SOURCES = analyzer.C 

evaluateConfs_SOURCES = evaluateConfs.C

PerturbativePredictions_SOURCES = PerturbativePredictions.C

HiggsMassFromLargeNf_SOURCES = HiggsMassFromLargeNf.C

PhaseDiagramAtLargeLambda_SOURCES = PhaseDiagramAtLargeLambda.C
