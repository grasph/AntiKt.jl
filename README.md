# Test implemenation of antikt fastjet algorithm in Julia

## Introduction
This code is an excercise to evaluate Julia for HEP application. The antikt hadronic jet clustering is implemented using the [fastjet](http://fastjet.fr/) tile N<sup>2</sup> algorithm&nbsp;\[[1](#Ref1)\].

We started with translating the C++ code extracted from the fastjet package. Then the code was profiled to improve the performance. We had little experience with the Julia programming language when translating the code.

It should be noted that the fastjet code has been written to run very fast, the code in addition to the algorithm itself has been highly optimized. The code use heavily pointers.

## Results

Timing has been peformed using 100 multijet events from proton-proton collisions at sqt{s} = 13TeV generated with [Pythia8](https://pythia.org)&nbsp;\[[2](#Ref2)\].

Hardware: laptop with an Intel(R) Core(TM) i7-8550U CPU, 8 GB of RAM.

Reference timing from C++ code: 263 ± 9 μs

**Time with the latest julia implementation: 315 ± 7 μs, 1.2x C++ implementation.**

JIT compiling is excluded from the measurement: function called once before the timing to trigger the compilation.

Optimisation history:

| Code version                                                              | Duration / event | Ratio of duration with C++ code |
|---------------------------------------------------------------------------|------------------|------------------------------------|
| First debugged implementation                                             |  1431 ± 50 μs    | 5.4 |
| Added a type annotations, that was missing in a struct                    |  1150 ± 32 μs    | 4.4 |
| Fixed type unstability of an iterator over a Union{Nothing, X} collection |   696 ± 14 μs    | 2.7 |
| Removed the Union{Nothing, X} (two of them)                               |   491 ± 14 μs    | 1.8 |
| Replaced a Vector by a static array                                       |   390 ±  9 μs    | 1.5 |
| Optimized tile neighbour iteration                                        |   331 ±  5 μs    | 1.3 | 
| Removed several boundary checks                                           |   315 ±  7 μs    | 1.2 | 
## Running the benchmark

### Prerequeries

The [HepMC3](https://gitlab.cern.ch/hepmc/HepMC3) is used by the c++ version of the code to read the data file. Install it in a directory you will call `externals` at the root directory of this project or set the directory where it is installed in the `Makefile.def` file. If you have access to LCG `cvmfs`, it can be found in `/cvmfs/sft.cern.ch/lcg/releases/hepmc3/<release>/<architecture>`

### One-line execution

A script to run the benchmark can be found in the directory `benchmark`.

   * `./run_benchmarks` will create the `events.hepmc3` input file and runs the benchmark on the C++ and Julia code.

### Running the individual code

Instead of using the `run_benchmarks` script the test codes can be executed directly. For this, you need first to:

1. Create a `events.hepmc3` or uncompress the one from `data/events.hepmc3`;
2. Either install the `AntiKt` and `HepMC3` packages located in the directory `antikjl` and `hepmc3jl` using the `dev` command of the Julia package manager or add these two directories in the `JULIA_LOAD_PATH` environment variable (as done in the `run_benchmark` script).

The Julia benchmark code is `benchmark/run_antikt_jl` and the c++ one is `antiktcxx/antikt`. The script `benchmark/run_antikt_cxx` can be used to ensure the c++ code is built and run it. Example of execution:

* `./run_antikit_jl` -m 100 `events.hepmc3`
* `./run_antikit_jl` -m 100 `events.hepmc3`

### Generating a data file

A compressed data file is located under the path [`data/events.hepmc3.gz`](data). The fata file can also be generated using the code from the [`genevts`](genevts)
 directory.
## References

<a id="Ref1"></a>\[1\] Pileup subtraction using jet areas, M.&nbsp;Cacciari and G.&nbsp;P.&nbsp;Salam, [Phys. Lett. B 641 (2006) 57](https://inspirehep.net/literature/755413)</br>
<a id="Ref2"></a>\[2\] A comprehensive guide to the physics and usage of PYTHIA 8.3, C.&nbsp;Bierlich et al., [arxXiv: 2203.11601 ](https://inspirehep.net/literature/755413)