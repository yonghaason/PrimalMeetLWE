# Primal Hybrid Estimator

Provide attack cost estimation for some hybrid attacks of (primal) lattice reduction and combinatorial approach, where the combinatorial approach is represented by "level" as folows:

- Lv 0: Primal + Exhaustive search (folklore)
- Lv 1: Primal + MitM [Howgrave-Graham, C07] with our optimization (See Sec 6.3)
- Lv $\ge$ 2: Primal + Meet-LWE [May, C21] = Ours

 
> This DON'T cover the non-hybrid (pure primal) attack [AGPW17], which is based on another clever idea that is incompatible with hybrid strategy.

## How to use?

In 'python (checked with 3.8)' shell (NOT 'sage' shell!):

```python
>>> from estimator import *
>>> # param = [(n), (q), (error_type), (error_param), (weight_param), (m)]
>>> # primal_may(param, (lv))
>>> param = [2**15, 2**768, 'gaussian', 3.19, 192, 2**15]
>>> primal_may(param, lv=2)
```

> It takes LONG LONG time (over several hours!) for high levels, or large LWE parameters. So if the reviewer is curious about the detailed parameters, we recommend to see the text files in `logs' directory.

 ### About the relation and consistency with `lattice-estimator'

Our level-0 attack strategy is exactly the same with `primal_hybrid' with options `mitm=False, babai=True' in lattice-estimator. We observe that both scripts output similar numbers for this strategy.

Meanwhile, our level-1 attack strategy is almost similar to the primal-MitM hytbrid attack [HG07], but has a small optimization on Babai NP call compared to it. 