#
# Isotropic Distortion of Gaussian field
#

EXEC      = theory_zeldovich
DIRS      = $(GSL_DIR)

all: $(EXEC)

INCLDIRS  = $(foreach dir, $(DIRS), -I$(dir)/include)
LIBDIRS   = $(foreach dir, $(DIRS), -L$(dir)/lib) 
CXXFLAGS := $(INCLDIRS)

OBJS1 := theory_zeldovich.o power.o growth.o corr.o

LIBS1 := $(LIBDIRS) -lm -lboost_program_options
LIBS1 += -lgsl -lgslcblas

theory_zeldovich: $(OBJS1)
	$(CXX) $(OBJS1) $(LIBS1) -o $@


corr.o: corr.cpp corr.h power.h
growth.o: growth.cpp
isotropic_gaussian.o: isotropic_gaussian.cpp power.h growth.h corr.h
power.o: power.cpp power.h

.PHONY: clean
clean :
	rm -f $(EXEC) $(OBJS1) $(OBJS2) $(OBJS3)

