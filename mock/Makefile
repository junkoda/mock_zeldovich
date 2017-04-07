#
# CXX and OPT defined in ../Makefile
#
EXEC      = mock_zeldovich

default:
	cd .. && $(MAKE) mock

mock: $(EXEC)


CXXFLAGS  += -I../lib $(OPT)

#
# Linking libraries
#
LIBS += -L../lib -lfs -lgsl -lgslcblas -lboost_program_options

ifeq (,$(findstring -DDOUBLEPRECISION, $(OPT)))
  # Single precision FFTW
  FFTWSUF=f
endif
LIBS += -lfftw3$(FFTWSUF) -lfftw3$(FFTWSUF)_mpi

ifdef OPENMP
  LIBS += -lfftw3$(FFTWSUF)_omp
  #LIBS += -lfftw3$(FFTWSUF)_threads # for thread parallelization instead of omp
endif


OBJS     = mock_zeldovich.o


mock_zeldovich: $(OBJS)
	$(CXX) $(OBJS) $(LIBS) -o $@

.PHONY: clean
clean :
	rm -f $(EXEC) $(OBJS1)
