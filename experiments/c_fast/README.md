
# Experiments for Primal-Meet-LWE

  

Tested on Ubuntu 22.04.1 with x86_64 architecture

  

## Compile Command

Requires cmake (>= 3.19)

```bash

mkdir  build && cd  build

cmake  ..

make

```

This will make executables `unif_test`, `lsh_simulation`, `cons_simulation` in the same directory.

  

### Experiment 1. Uniformity Test

This experiment checks whether $\{[Ms]_B : HW(s) = w\}$ distributes uniformly at random in $\mathcal P(B^*)$, by comparing the ratio of volumes and the number of points.

```bash
./unif_test
```

will give how to insert arguments, for example

```bash
./unif_test -d 15 -h 6 -l 2 -q 8
```

>  **Caution**: The box length $\ell$ determines a box $[-\ell, \ell)$.

  

### Experiment 2. The probability of correct list definition
 
In Section 5, we show that the number of representation pairs $N_{pairs}(e) \sim 2\cdot B(R/2, p_{rep}(e))$ where with the ambiguity $R$ (under the uniform and independent heuristic). To estimate the overall success probability, we need to compute the expectation of the probability $\Pr[N_{pairs} \neq 0] = 1 - (1 - p_{rep}(e))^{R/2}$ where $e \leftarrow \mathcal D$ for some noise distribution $\mathcal D$, say $E_{e\leftarrow \mathcal D}\left[ 1 - (1 - p_{rep}(e))^{R/2} \right].$  However, the formula is extremely complicated and we instead approximate it by $$E_{e\leftarrow \mathcal D}\left[ 1 - (1 - p_{rep}(e))^{R/2} \right] \approx 1 - \left(1 - E_{e\leftarrow \mathcal D}[p_{rep}(e)]\right)^{R/2},$$ and use the RHS for our attack cost estimation. 
This experiment empirically computes the LHS and compare it with the RHS: Based on the results, we choose the constraint parameters that makes this approximation sufficiently close, which makes our attack cost estimation based on the RHS formula reliable to some extent.

```bash
./cons_simulation
```

will give how to insert arguments, for example

```bash
./cons_simulation -q 1024 -r 50 -l 20 -u 0 -b 2 --multiple 3
```

  

### Experiment 3. The success probability of near-collision finding  

In Section 6, our (torus-LSH based) near-collision finding algorithm with $R$ repetitions finds a near-collision pair (of distance e) with probabiltiy $1 - (1-p_{good}(e))^R$. Again, we approximate $$E_{e\leftarrow \mathcal D}\left[ 1 - (1 - p_{good}(e))^{R} \right] \approx 1 - \left(1 - E_{e\leftarrow \mathcal D}[p_{good}(e)]\right)^{R},$$ and choose $R = 3/E[p_{good}(e)]$ which makes the RHS $\approx 1 - e^{-3} \approx  0.95$
This experiment empirically computes the LHS and compare it with the RHS. Based on the results, we choose the LSH parameters that makes this approximation sufficiently close, in order to safely claim that we finds a near-collision with high probability, say $0.95$.

```bash
./lsh_simulation
```

will give how to insert arguments, for example

```bash
./lsh_simulation -q 1024 -r 50 -l 20 -u 1 -b 3 --multiple 3
```