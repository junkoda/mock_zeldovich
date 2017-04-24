#
# Linear limit of zeldovich mock with discritised grid
#

EXEC      = theory_zeldovich_linear
DIRS      = $(GSL_DIR)

all: $(EXEC)

INCLDIRS  = $(foreach dir, $(DIRS), -I$(dir)/include)
LIBDIRS   = $(foreach dir, $(DIRS), -L$(dir)/lib) 
CXXFLAGS := $(INCLDIRS) -O2

OBJS := theory_zeldovich_linear.o power.o

LIBS := $(LIBDIRS) -lm -lboost_program_options
LIBS += -lgsl -lgslcblas

theory_zeldovich_linear: $(OBJS)
	$(CXX) $(OBJS) $(LIBS) -o $@


corr.o: corr.cpp corr.h power.h
growth.o: growth.cpp
power.o: power.cpp power.h

.PHONY: clean
clean :
	rm -f $(EXEC) $(OBJS)

