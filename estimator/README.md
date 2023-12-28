# Primal Hybrid Estimator

Provide attack cost estimation for some hybrid attacks of (primal) lattice reduction and combinatorial approach, where the combinatorial approach is represented by "level" as folows:

- Lv 1: Primal + MitM [Howgrave-Graham, C07]
- Lv $\ge$ 2: Ours = Primal + Meet-LWE [May, C21]
 
> This DON'T cover the non-hybrid (pure primal) attack [[AGVW17](https://eprint.iacr.org/2017/815.pdf)], which is based on another clever idea that is incompatible with hybrid strategy. For this attack estimation, we used [lattice-estimator](https://github.com/malb/lattice-estimator/) that covers this attack as `usvp`.

## How to use?

In 'python (checked with 3.8)' shell (NOT 'sage' shell!):

```python
>>> from estimator import *
>>> # param = [(n), (q), (error_type), (error_param), (weight_param), (m)]
>>> # primal_may(param, (lv))
>>> param = [2**15, 2**768, 'gaussian', 3.19, 192, 2**15]
>>> _ = primal_may(param, t=2)
```

> It may takes LONG time (can be several hours!) especially when the weight is large. So if the reviewer is just curious about the detailed parameters, we recommend to see the log files in `logs' directory.