include Makefile.def

.PHONY: show_setup

all:
	make -C genevts all
	make -C antiktcxx all

clean:
	make -C genevts clean
	make -C antiktcxx clean

show_setup:
	@echo HepMC3 software directory: $(HEPMC3_DIR)
	@echo Pythia8 software directory: $(PYTHIA_DIR)
	@echo FasetJet software directory: $(FASTJET_DIR)
