SRCS = ConfigReader.cxx FlowAnalyzer.cxx
OBJS = $(SRCS:.cxx=.o)
DEPS = FlowUtils.h
TARGET = FlowAnalyzer

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
ROOTLIBS     += -lEG

#INCFLAGS = -I$(ROOTSYS)/include -I/afs/rhic.bnl.gov/star/packages/SL19b/StRoot -I/afs/rhic.bnl.gov/star/packages/SL19b/StRoot/StPicoEvent -I/afs/rhic.bnl.gov/star/packages/SL19b/StRoot/StEpdUtil
INCFLAGS = -I$(ROOTSYS)/include -I./ -I./StRoot -I./StRoot/StEvent -I./StRoot/StBichsel
LIBFLAGS = -L$(ROOTSYS)/lib -L./ -L./libs -Wl,-R./libs
SOLIBS = -lStPicoDst -lStEpdUtil -lStBichsel
#for each non-standard dynamic library location -L a corresponding -Wl,-R should be specified
#the -L's to StPicoEvent, StEpdUtil, and StBichsel are for locating libStPicoEvent.so, libStEpdUtil.so, and libStBichsel.so
#the -l's are to link these .so files

#-L./StRoot/StPicoEvent -Wl,-R./StRoot/StPicoEvent -L./StRoot/StEpdUtil -Wl,-R./StRoot/StEpdUtil

CC = g++ #-m32
FLAGS = -Wall -g -fPIC

COMPILE = $(CC) $(FLAGS)


all: $(TARGET)

$(TARGET): $(OBJS)
	$(COMPILE) -o $@ $^ $(SOLIBS) $(ROOTCFLAGS) $(ROOTLIBS) $(INCFLAGS) $(LIBFLAGS)

FlowAnalyzer.o: FlowUtils.h ConfigReader.h FlowAnalyzer.cxx
	$(COMPILE) -o $@ -c FlowAnalyzer.cxx $(SOLIBS) $(ROOTCFLAGS) $(ROOTLIBS) $(INCFLAGS) $(LIBFLAGS)

ConfigReader.o: ConfigReader.h ConfigReader.cxx
	$(COMPILE) -o $@ -c ConfigReader.cxx $(ROOTCFLAGS) $(ROOTLIBS)

.PHONY: clean

clean:
	rm -f *.so $(TARGET) $(OBJS)



# target : dependencies
#	action

# $@ = current target
# $^ = current dependencies
# $< = name of the related file that caused the action
