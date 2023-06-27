import numpy as np
from math import log

# from dependence_test import *
# dep_test(l = 5, r = 20, c = 40, b = 1, Rscale = 3)

def log2(a):
  return log(a, 2)

def p_adm_unif(b, ell):
  if b <= ell:
    return 1.0 - b/(2*ell)
  else:
    return ell/(2*b)


def dep_test(l, r, c, b, Rscale = 1, repeat=10000): 
  E_p_rep = 1.0
  for i in range(r):
    E_p_rep *= p_adm_unif(l, 2*c)
  
  scaler = 2**-10

  E_p_rep *= scaler
  
  print("******************************")
  print("p_sp(v): Prob of at least one rep pairs of v exists in [-%d, %d]^%d (in volume %d space)" % (c, c, r, 1/scaler))
  print("- v: sampled from U([-%d, %d]^%d)" % (l, l, r))  
  print("- e: sampled from Gaussian of std.dev", b)
  print()
  
  print("* Side Info")
  print("Log E[p_rep]: %.3f" % log2(E_p_rep))
  R = int(Rscale/E_p_rep)
  print("R =", Rscale, "* 1/p_rep =", R)
  print()
  
  # Experiment 1. Probability of v, v + e
  E_p_sp1 = 0.0
  E_p_sp2 = 0.0
  E_p_L = 0.0

  error = np.zeros(r)

  for i in range(r):
    error[i] = np.random.normal(scale=b)

  for _ in range(repeat):
    v1 = np.zeros(r)
    v2 = np.zeros(r)
    for i in range(r):
      while 1:
        v1[i] = np.random.uniform(-l, l)
        v2[i] = v1[i] + error[i]
        if abs(v2[i]) <= l:
          break

    p_rep1 = scaler
    p_rep2 = scaler
    for i in range(r):
      p_rep1 *= 1.0 - abs(v1[i]) / (2*c)
      p_rep2 *= 1.0 - abs(v2[i]) / (2*c)

    p_sp1 = 1-(1 - p_rep1)**R
    p_sp2 = 1-(1 - p_rep2)**R

    E_p_sp1 += p_sp1
    E_p_sp2 += p_sp2
    E_p_L += p_sp1*p_sp2

  print("Experiment 1: p_sp(v) and p_sp(v+e)")  
  print("- e in [%.3f, %.3f]" % (error.min(), error.max())) 
  print("- E[p_sp]: %.4f / %.4f" % (E_p_sp1/repeat, E_p_sp2/repeat))
  print("- E[p_L]: %.4f" % (E_p_L/repeat))
  print("- Log E[p_sp]: %.4f / %.4f" % (log2(E_p_sp1) - log2(repeat), log2(E_p_sp2) - log2(repeat)))
  print("- Log E[p_L]: %.4f" % (log2(E_p_L) - log2(repeat)))
  print()

  # Experiment 2. Probability of v1, v2 (indep)
  E_p_sp1 = 0.0
  E_p_sp2 = 0.0
  E_p_L = 0.0

  for _ in range(repeat):
    v1 = np.zeros(r)
    v2 = np.zeros(r)
    for i in range(r):
      v1[i] = np.random.uniform(-l, l)
      v2[i] = np.random.uniform(-l, l)

    p_rep1 = scaler
    p_rep2 = scaler
    for i in range(r):
      p_rep1 *= 1.0 - abs(v1[i]) / (2*c)
      p_rep2 *= 1.0 - abs(v2[i]) / (2*c)

    p_sp1 = 1-(1 - p_rep1)**R
    p_sp2 = 1-(1 - p_rep2)**R

    E_p_sp1 += p_sp1
    E_p_sp2 += p_sp2
    E_p_L += p_sp1*p_sp2
  
  print("Experiment 2: p_sp(v) and p_sp(v') for indep v and v'")
  print("- E[p_sp]: %.4f / %.4f" % (E_p_sp1/repeat, E_p_sp2/repeat))
  print("- E[p_L]: %.4f" % (E_p_L/repeat))
  print("- Log E[p_sp]: %.4f / %.4f" % (log2(E_p_sp1) - log2(repeat), log2(E_p_sp2) - log2(repeat)))
  print("- Log E[p_L]: %.4f" % (log2(E_p_L) - log2(repeat)))