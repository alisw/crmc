SYSTEM=$(shell uname)
SYSTEM2=$(shell uname -m)
m32=$(shell if [ "$(SYSTEM2)" = "i686" ]; then echo "-m32"; fi)

SOURCE_DIR=src/

VER=199

LIBDIR=obj/

CFLAGS     = $(m32) -ggdb3 -pipe -Wall
CXXFLAGS   = $(m32) -ggdb3 -pipe -Wall
CXXFLAGS  += -Wpacked -malign-double -mpreferred-stack-boundary=8

ifndef CXX
  CXX           = g++
endif
ifndef CC
  CC            = gcc
endif

ROOTCFLAGS   = $(shell root-config --cflags)
ROOTLIBS     = $(shell root-config --libs) -lEG

HEPCFLAGS    = -I$(HEP_ROOT)/include
HEPLIBS      = -L$(HEP_ROOT)/lib -lHepMC -lHepMCfio
BOOSTLIBS    = -lboost_iostreams -lboost_system

CFLAGS       = $(ROOTCFLAGS) $(HEPCFLAGS)
LIBS         = $(ROOTLIBS) $(HEPLIBS) $(BOOSTLIBS)

CXXFILES=analysis.cc

OBJS=$(FILES:%.f=$(LIBDIR)%.o)
CXXOBJS=$(CXXFILES:%.cc=$(LIBDIR)%.o)


all: check dirs bin_dir bin/analyse

.PHONY : check
check:
	@if [ -z "$(HEP_ROOT)" ]; then echo "Please set HEP_ROOT to the root directory of HepMC2"; exit 1; fi
	@if [ -z "$(ROOTSYS)" ]; then echo "Please set ROOTSYS to the root directory of root"; exit 1; fi

bin/analyse: $(OBJS) $(CXXOBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(CXXOBJS) -o $@ $(LIBS)
	@(echo "")
	@(echo "==> compilation successful!")
	@(echo "==> type bin/analyse")
	@(echo "")


$(OBJS) : $(LIBDIR)%.o : $(SOURCE_DIR)%.f
	$(FC) -o $@ -c $<
$(CXXOBJS) : $(LIBDIR)%.o : $(SOURCE_DIR)%.cc src/analysis.h
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ -c $<

dirs:
	@if [ ! -d $(LIBDIR)  ]  ;then \
	set -x; mkdir -p $(LIBDIR); set +x;\
	fi

bin_dir:
	@if [ ! -d bin  ]  ;then \
	set -x; mkdir -p bin; set +x; \
	fi

clean:
	rm -rf $(LIBDIR)
	rm -rf bin
	rm -f *.copy *.check fort*
	rm -f IO_GenEvent.dat
	rm -f *~ \#*\# */*~

#;;; Emacs Setup
#;;; Local Variables:
#;;; mode: Makefile
#;;; End:
