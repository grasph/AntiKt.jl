.PHONY: all

include ../Makefile.def

CXXFLAGS+=-MMD
CXXFLAGS+=-O2
#CXXFLAGS+=-g -Og

all: genevts

clean:
	$(RM) genevts *.o *.d

SRC = genevts.cc

DEPS:=$(addsuffix .d,$(basename $(SRC)))

genevts: genevts.cc
	$(LINK.cc) -o $@ -I $(PYTHIA_DIR)/include  -L $(PYTHIA_DIR)/lib -Wl,--rpath $(PYTHIA_DIR)/lib $<  -lpythia8 -lHepMC3 


-include $(DEPS)
