# Experiments for Primal-Meet-LWE

### Compile c++ experiments

```bash
mkdir build && cd build
cmake ..
make
```

This will make executables `expectation_exchange`, `ncf_uniform_test`, `split_uniform_test`, `meet_lwe_matrix` in the same directory.

## C. On Continuous vs. Discrete Gaussian

This experiment checks the computation of $p_{np}$ and $p_{adm}$ based on the herustic approximation 
$$\zeta_{Q}(\mathcal{D}_{\Z, \sqrt{2\pi}\sigma}) \approx \mathcal{G}_\sigma$$.

To run this experiment, run *gaussian_approx.sage* using [sagemath](https://www.sagemath.org/)

## E. Experimental Validation for Independence Heuristic
This experiment checks the lower bound 
$$ \Pr[E_{sp}^{(i)}] \ge \prod_{k\in [2^i]}\Pr[E_{sp, k}^{(i)}]$$
To reproduce the experiment in the paper, in **python** shell,
```python
from dependency_test import *
dep_test_lv1(l1=6, l2=6, r=50, error_param=1)
dep_test_lv2(l1=6, l2=6, l3=6, r=50, error_param=1)
```

## F. Toy Implementation of Meet-LWE (for matrix modulus)
This is a (proof-of-concept) implementation of Meet-LWE for matrix modulus (Algorithm 3 in the paper).
The experiment in the paper can be reproduced as follows:
First, run the following command
```bash
./meet_lwe_matrix -m 25 -d 20 -h 12 -q 30 -t 2 -C 2 --repeat 100
```
and then insert the split weights
```bash
> Insert w[1]: 8
> Insert w[2]: 4
```
