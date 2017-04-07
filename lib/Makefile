#
# fs: fast cosmological simulation system
#   


CXX := mpic++ -std=c++11
OPT := -DDOUBLEPRECISION

LIB = libfs.a #libfs.dylib
lib: $(LIB)


#
# Compile options
#

#OPT += -DCHECK  # slow assersions (turn on if you encounter a problem)
#OPT += -DNDEBUG # drop assertion (not recomended)

#
# Compiler configurations
#

# CXX is defined in ../Makefile

WOPT     ?= -Wall
CPPFLAGS += -O3 $(INCLS) $(WOPT) $(OPENMP) $(OPT)
#LIBS     += -lm #-lhdf5

OBJS := comm.o msg.o config.o fft.o mem.o particle.o util.o
OBJS += power.o cosmology.o lpt.o 
#OBJS += gadget_file.o hdf5_write.o





#
# Library libfs
#
libfs.a: $(OBJS)
	ar r $@ $(OBJS)

libfs.dylib: $(OBJS)
	echo $(OPT)
	$(CXX) -dynamiclib $(OBJS) -o $@


#
# Library
#

# Dependences
$(OBJS): Makefile

# g++ -MM -MG *.cpp
comm.o: comm.cpp msg.h comm.h error.h
config.o: config.cpp config.h msg.h
cosmology.o: cosmology.cpp msg.h const.h cosmology.h error.h
fft.o: fft.cpp config.h comm.h msg.h error.h mem.h util.h particle.h \
  fft.h
gadget_file.o: gadget_file.cpp msg.h comm.h error.h cosmology.h \
  gadget_file.h particle.h config.h
hdf5_write.o: hdf5_write.cpp config.h msg.h comm.h error.h cosmology.h \
  hdf5_io.h particle.h
lpt.o: lpt.cpp comm.h msg.h error.h mem.h config.h cosmology.h power.h \
  particle.h fft.h lpt.h
mem.o: mem.cpp config.h msg.h error.h util.h particle.h mem.h fft.h
msg.o: msg.cpp comm.h msg.h error.h
particle.o: particle.cpp msg.h comm.h error.h util.h particle.h config.h \
  fft.h mem.h
power.o: power.cpp comm.h msg.h error.h power.h
util.o: util.cpp util.h particle.h config.h


.PHONY: clean run dependence
clean:
	rm -f $(LIB) $(OBJS)

dependence:
	g++ -MM -MG *.cpp

