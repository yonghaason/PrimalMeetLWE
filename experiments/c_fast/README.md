# Experiments for Primal-Meet-LWE

Tested on Ubuntu 22.04.1 with x86_64 architecture

## Compile Command
Requires cmake (>= 3.19)
```bash
mkdir build && cd build
cmake ..
make
```
This will make executables `adm_test`, `cons_test`, `nc_test` in the same directory.

### Experiment 1. Admissible Probability

A simple test for admissible probability.
```bash
./adm_test
```
will give how to insert arguments, for example
```bash
./adm_test -l 10 -e 20 -u 1
```
**Caution**: The box length $\ell$ is slightly different from the paper definition (where $2\ell$ is the box length).

### Experiment 2. Number of pairs in constraints

From Section 5.1, we claim that $R \cdot p_{rep}$ representation pairs lie in $\ell$-hypercube.  
This experiment validates this claim (and hence also indirectly validates underlying Heuristic 1).
```bash
./cons_test
```
will give how to insert arguments, for example
```bash
./cons_test -d 30 -h 10 -w 6 -u 0 -e 3 -l 10 -t 10
```

### Experiment 3. LSH-based Near-collision Finding

From Section 6, we choose the repetition number $R = \lceil 3/p_{good} \rceil$ in Algorithm 3 (LSH-based NCF)
(which comes from the assumption that error is freshly sampled for each repetition). 
This experiment shows such choice of $R$ really finds a sufficient portion of near-collisions.
```bash
./nc_test
```
will give how to insert arguments, for example
```bash
./nc_test -r 25 -e 3 -b 6 -d 15 -u 1
```