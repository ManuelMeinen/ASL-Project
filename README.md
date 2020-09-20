# ASL Project 2020
### *ELLIPTIC CURVE PUBLIC KEY GENERATION: OPTIMIZING BIG NUMBERARITHMETIC IN RADIX REPRESENTATIONS FOR USE IN CURVE25519*

## Description
This repository contains multiple implementations of a scalar multiplication of a point on an elliptic curve (Curve25519). The different implementations are located in different folders. Each implementation uses either a different radix representation or different optimizations. The implementations called `opt_all` contain the code where all suitable optimizations were combined in order to achieve the highest possible performance. The implementations called `base` contain a naive implementation in a given radix representation.

Each implementation has a folder containing the 256-bit arithmetic functions and a folder containing the curve computations.

The folder `reference_implementations` contains several implementations, which we used as a reference for our project.

There are some additional folders containing things like plots, conversion tools etc.
The implementation is not tested to be safe against side-channel attacks.

## How to run it
From the folder `benchmark_framework` run the following shell script:

`$ ./run_benchmarking.sh`

In the shell script `run_benchmarking.sh` you can set the following parameters in order to configure the benchmarking framework.

* `compiler_flag_ids` the flags with which gcc compiles the code (no-opt, o2,  o3-no-vect, o3)
* `test_set_numbers` specifies the test sets that we use for the benchmarking. (1-5)
* `func_ids` the name of the implementations, which we want to benchmark. (base_256_17, opt_scalar_replacement, opt_scalar_repl_and_vectorize, opt_curve, opt_all, base_256_25, opt_all_25, base_256_51, opt_all_51)
* `run_mode` do you want to validate (v), benchmark (b) or do both i.e. full benchmarking (f)
* `be_verbose` print additional information if 1 and only the most essential information if 0.

## Report
A written report can be found in this repository. ([Report](31_report.pdf))

## Contributors
* Robin Burkhard
* Manuel Meinen
* Sabina Fischlin
* Supraja Sridhara
