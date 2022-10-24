## Directory Contents
   This direcrory contains example code to run antikt algorithm on a HepMC3 file. Clustering parameters are fixed to distance parameter $R=0.4$ and minumum transverse momentum $p_{\textrm{T}}^{\textrm{min}} = 5\,\textrm{GeV}$.
   
   * `antikt` standalone version
   * `antikt-fastjet` version using the Fastjet shared library
   
Note: both codes are using Fastjet code. In the standalone version, the relevant code was copied in the single source file.
   
## Code dependencies

   * [`Fastjet`](https://www.fastjet.fr) for `antikt-fastjet`;
   * [`HepMC3`](https://gitlab.cern.ch/ hepmc/HepMC3).
   
The location of the external library is set in `Makefile.def`.
 