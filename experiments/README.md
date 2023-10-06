# Experiments for Primal-Meet-LWE

### Compile c++ experiments

```bash
mkdir build && cd build
cmake ..
make
```

This will make executables `expectation_exchange`, `ncf_uniform_test`, `split_uniform_test`, `meet_lwe_matrix` in the same directory.
  
## Section 9.1 Uniformity Modeling
  
### Uniformity Test 1
 
This experiment checks the uniformity of the set 
$$P_1(s_k^{(0)}) = \{[Ms_r]_{B, r^{(1)}} ~|~ (s_r, s_r'): w^{(1)}\text{-rep pair of~} s_k^{(0)}\} .$$
The experiment in the paper can be reproduced by
```bash
./split_uniform_test -m 46 -d 8488 -h 20 -w 14 --lastGS 331.8 --rhf 1.0048 # from q = 2^699 and beta = 300
```  

### Uniformity Test 2
This experiment checks the uniformity of the set 
$$L^{(1)} = \{[Mx]_{B, r^{(0)}} ~|~ x \in \mathcal{T}^d(w^{(1)}) \wedge \|[Ms]_{B, r^{(1)}}\| \le \ell^{(1)}\} .$$
The experiment in the paper can be reproduced by
```bash
./ncf_uniform_test -m 46 -d 8488 -w 14 -r 19 -b 18 --lastGS 331.8 --rhf 1.0048 # from q = 2^699 and beta = 300
```

## Section 9.2 Exchange of Expectation
This experiment measures the average $E_{x\leftarrow \mathcal D}\left[ 1 - (1 - p_{r}(x))^{R} \right]$ for $R \approx C/E_{x\leftarrow \mathcal D}[p_r(x)]$.
The experiments in the paper can be reproduced by
```bash
./expectation_exchange -r 50 -l 6 -b 1 -u 0 -C 10 # Gaussian distribution
./expectation_exchange -r 50 -l 6 -b 6 -u 1 -C 10 # Uniform distribution
```
while varying the dimension $r$.

## Section 9.3 Lower Bound by Neglicting Relation
This experiment checks the lower bound 
$$ \Pr[E_{sp}^{(i)}] \ge \prod_{k\in [2^i]}\Pr[E_{sp, k}^{(i)}]$$
To reproduce the experiment in the paper, in **python** shell,
```python
from dependency_test import *
dep_test_lv1(l1=6, l2=6, r=50, error_param=1)
dep_test_lv2(l1=6, l2=6, l3=6, r=50, error_param=1)
```

## Section 9.4 Meet-LWE for matrix modulus
This is a (proof-of-concept) implementation of Meet-LWE for matrix modulus (Algorithm 5 in the paper).
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