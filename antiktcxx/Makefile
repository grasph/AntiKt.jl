include ../Makefile.def

.PHONY: all clean

CXXFLAGS+=-MMD
CXXFLAGS+=-O2 -ggdb
CXXFLAGS+=-ggdb #For perf profiler
#CXXFLAGS+=-ggdb -Og

all: antikt-fastjet antikt

clean:
	$(RM) antikt-fastjet genevts antikt *.o *.d

SRC = antikt.cc  antikt-fastjet.cc  genevts.cc  PseudoJet.cc

DEPS:=$(addsuffix .d,$(basename $(SRC)))

antikt-fastjet: antikt-fastjet.cc
	$(LINK.cc) -o $@ -Wno-deprecated-declarations  -I $(FASTJET_DIR)/include -L $(FASTJET_DIR)/lib -Wl,--rpath $(FASTJET_DIR)/lib \
-I $(HEPMC3_DIR)/include -L $(HEPMC3_DIR)/lib -Wl,--rpath $(HEPMC3_DIR) $< -l fastjet -lHepMC3 -lm

antikt: antikt.o PseudoJet.o
	$(LINK.cc) --std=c++17  -o $@ -L $(HEPMC3_DIR)/lib -Wl,--rpath $(HEPMC3_DIR)/lib $^ -lHepMC3 -lm

antikt.o: antikt.cc
	$(COMPILE.cc) --std=c++17 -o $@ -I $(HEPMC3_DIR)/include/ $<

-include $(DEPS)

