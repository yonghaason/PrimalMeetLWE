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
  print("- E[p_L]=E[p_sp1*p_sp2]: %.4f" % (E_p_L/repeat))
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
  print("- E[p_L]=E[p_sp1*p_sp2]: %.4f" % (E_p_L/repeat))
  print("- Log E[p_sp]: %.4f / %.4f" % (log2(E_p_sp1) - log2(repeat), log2(E_p_sp2) - log2(repeat)))
  print("- Log E[p_L]: %.4f" % (log2(E_p_L) - log2(repeat)))

# double_dep_test(l1 = 5, l2 = 10, r = 20, c = 40, b = 1, Rscale = 2)

def double_dep_test(l1, l2, r, c, b, Rscale = 1, repeat=10000):
  assert l1 <= l2, "l1 should be smaller than l2"

  E_p_rep = 1.0
  for i in range(r):
    E_p_rep *= p_adm_unif(l2, 2*c)
  
  scaler = 2**-10

  E_p_rep *= scaler

  print("******* Goal of this experiment **********")
  print("p_sp(v) := At least one rep pairs of v in [-%d, %d]^%d (in volume %d space)" % (c, c, r, 1/scaler))
  print("Compare E[p_L] = E[p_sp(e_21)*p_sp(e_22)*p_sp(e_23)*p_sp(e_24)]")
  print("     where e_11 = e_21 - e_22 / e_12 = e_23 - e_24")
  print("           e_1*: sampled from U([-%d, %d]^%d)" % (l1, l1, r))  
  print("           e_2*: sampled from U([-%d, %d]^%d)" % (l2, l2, r))  
  print("******************************")
  print()
  
  print("* Side Info")
  print("Log E[p_rep]: %.3f" % log2(E_p_rep))
  R = int(Rscale/E_p_rep)
  print("R =", Rscale, "* 1/p_rep =", R)
  print()  
  
  # Experiment 1: e_11 & e_12 are dependent (e_0 = e_11 - e_12)
  E_p_sp1 = 0.0; E_p_sp2 = 0.0; E_p_sp3 = 0.0; E_p_sp4 = 0.0
  E_p_L = 0.0

  e_0 = np.zeros(r)

  for i in range(r):
    e_0[i] = np.random.normal(scale=b)

  for _ in range(repeat):
    e_11 = np.zeros(r); e_12 = np.zeros(r)
    for i in range(r):
      while 1:
        e_11[i] = np.random.uniform(-l1, l1)
        e_12[i] = e_11[i] + e_0[i]
        if abs(e_12[i]) <= l1:
          break

    e_21 = np.zeros(r); e_22 = np.zeros(r); e_23 = np.zeros(r); e_24 = np.zeros(r)
    for i in range(r):
      while 1:
        e_21[i] = np.random.uniform(-l2, l2)
        e_22[i] = e_21[i] + e_11[i]
        if abs(e_22[i]) <= l2:
          break
      while 1:
        e_23[i] = np.random.uniform(-l2, l2)
        e_24[i] = e_23[i] + e_12[i]
        if abs(e_24[i]) <= l2:
          break

    p_rep1 = scaler; p_rep2 = scaler; p_rep3 = scaler; p_rep4 = scaler
    for i in range(r):
      p_rep1 *= 1.0 - abs(e_21[i]) / (2*c)
      p_rep2 *= 1.0 - abs(e_22[i]) / (2*c)
      p_rep3 *= 1.0 - abs(e_23[i]) / (2*c)
      p_rep4 *= 1.0 - abs(e_24[i]) / (2*c)

    p_sp1 = 1-(1 - p_rep1)**R; p_sp2 = 1-(1 - p_rep2)**R; p_sp3 = 1-(1 - p_rep3)**R; p_sp4 = 1-(1 - p_rep4)**R

    E_p_sp1 += p_sp1; E_p_sp2 += p_sp2; E_p_sp3 += p_sp3; E_p_sp4 += p_sp4     
    E_p_L += p_sp1*p_sp2*p_sp3*p_sp4

  print("* Experiment 1: e_11 - e_12 = e_0 (Gaussian of std.dev %f)" % b)  
  print("- e_0 in [%.3f, %.3f]" % (e_0.min(), e_0.max())) 
  print("- E[p_sp]: %.4f / %.4f / %.4f / %.4f" % (E_p_sp1/repeat, E_p_sp2/repeat, E_p_sp3/repeat, E_p_sp4/repeat))
  print("- E[p_L]: %.4f" % (E_p_L/repeat))
  print("- Log E[p_sp]: %.4f / %.4f / %.4f / %.4f" % (log2(E_p_sp1) - log2(repeat), log2(E_p_sp2) - log2(repeat), log2(E_p_sp3) - log2(repeat), log2(E_p_sp4) - log2(repeat)))
  print("- Log E[p_L]: %.4f" % (log2(E_p_L) - log2(repeat)))
  print()

  # Experiment 2: e_11 & e_12 are independent
  E_p_sp1 = 0.0; E_p_sp2 = 0.0; E_p_sp3 = 0.0; E_p_sp4 = 0.0
  E_p_L = 0.0

  for _ in range(repeat):
    e_11 = np.zeros(r); e_12 = np.zeros(r)
    for i in range(r):
      while 1:
        e_11[i] = np.random.uniform(-l1, l1)
        e_12[i] = np.random.uniform(-l1, l1)
        if abs(e_12[i]) <= l1:
          break

    e_21 = np.zeros(r); e_22 = np.zeros(r); e_23 = np.zeros(r); e_24 = np.zeros(r)
    for i in range(r):
      while 1:
        e_21[i] = np.random.uniform(-l2, l2)
        e_22[i] = e_21[i] + e_11[i]
        if abs(e_22[i]) <= l2:
          break
      while 1:
        e_23[i] = np.random.uniform(-l2, l2)
        e_24[i] = e_23[i] + e_12[i]
        if abs(e_24[i]) <= l2:
          break

    p_rep1 = scaler; p_rep2 = scaler; p_rep3 = scaler; p_rep4 = scaler
    for i in range(r):
      p_rep1 *= 1.0 - abs(e_21[i]) / (2*c)
      p_rep2 *= 1.0 - abs(e_22[i]) / (2*c)
      p_rep3 *= 1.0 - abs(e_23[i]) / (2*c)
      p_rep4 *= 1.0 - abs(e_24[i]) / (2*c)

    p_sp1 = 1-(1 - p_rep1)**R; p_sp2 = 1-(1 - p_rep2)**R; p_sp3 = 1-(1 - p_rep3)**R; p_sp4 = 1-(1 - p_rep4)**R

    E_p_sp1 += p_sp1; E_p_sp2 += p_sp2; E_p_sp3 += p_sp3; E_p_sp4 += p_sp4     
    E_p_L += p_sp1*p_sp2*p_sp3*p_sp4
  
  print("* Experiment 2: e_11 & e_12 are independent")
  print("- E[p_sp]: %.4f / %.4f / %.4f / %.4f" % (E_p_sp1/repeat, E_p_sp2/repeat, E_p_sp3/repeat, E_p_sp4/repeat))
  print("- E[p_L]: %.4f" % (E_p_L/repeat))
  print("- Log E[p_sp]: %.4f / %.4f / %.4f / %.4f" % (log2(E_p_sp1) - log2(repeat), log2(E_p_sp2) - log2(repeat), log2(E_p_sp3) - log2(repeat), log2(E_p_sp4) - log2(repeat)))
  print("- Log E[p_L]: %.4f" % (log2(E_p_L) - log2(repeat)))
  print()