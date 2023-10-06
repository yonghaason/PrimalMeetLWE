import numpy as np
from math import log

# from dependence_test import *
# dep_test_lv1(l1=6, l2=6, r=20 error_param=1, error_unif=False)
# dep_test_lv2(l1=6, l2=6, l3=6, r=20, error_param=1, error_unif=False)

def log2(a):
  return log(a, 2)

def p_adm_unif(b, ell):
  if b <= ell:
    return 1.0 - b/(2*ell)
  else:
    return ell/(2*b)

def dep_test_lv1(l1, l2, r, error_param, error_unif=False, C=10, repeat=10000, scaler=1): 
  assert l1 <= l2, "l should increase"
  E_p_rep = p_adm_unif(l1, 2*l2)**r

  print("******* Goal of this experiment **********")
  print("Define p(e) := Product of [1.0 - |e_i| / (2*l2)]")
  print("       & p_sp(e, R) := 1 - (1 - p(e))^R")
  print("Compare E[p_L] = E[p_sp(e_1)*p_sp(e_2)]")
  print()
  print("* Side Info")
  R = int(C/E_p_rep)
  print("Log E[p_rep] = %.3f, R = %d/p_rep = %d" % (log2(E_p_rep), C, R))
  if scaler != 1:
    R = int(R / scaler)
    E_p_rep *= scaler
    print("(R is scaled by %.4f -> R = %d)" % (scaler, R))
  print("******************************")
  print()
  
  print("* Experiment 1: e: sampled from", end = "")
  if error_unif: 
    print(" U[%.2f, %.2f])^%d " % (-error_param, error_param, r))
  else:
    print(" G(%.2f)^%d" % (error_param, r))
  print("              e_1: sampled from U([-%d, %d]^%d)" % (l1, l1, r))  
  print("              e_2 = e_1 + e")

  E_p_sp1 = 0.0; E_p_sp2 = 0.0
  E_p_L = 0.0

  for _ in range(repeat):
    errors = np.zeros(r)
    if error_unif:
      for i in range(r):
        errors[i] = np.random.uniform(-error_param, error_param)
    else:
      for i in range(r):
        errors[i] = np.random.normal(scale=error_param)

    v1 = np.zeros(r)
    v2 = np.zeros(r)
    for i in range(r):
      while 1:
        v1[i] = np.random.uniform(-l1, l1)
        v2[i] = v1[i] + errors[i]
        if abs(v2[i]) <= l1:
          break

    p_rep1 = scaler
    p_rep2 = scaler
    for i in range(r):
      p_rep1 *= 1.0 - abs(v1[i]) / (2*l2)
      p_rep2 *= 1.0 - abs(v2[i]) / (2*l2)

    p_sp1 = 1-(1 - p_rep1)**R
    p_sp2 = 1-(1 - p_rep2)**R

    E_p_sp1 += p_sp1
    E_p_sp2 += p_sp2
    E_p_L += p_sp1*p_sp2

  print("- E[p_sp]: %.4f / %.4f" % (E_p_sp1/repeat, E_p_sp2/repeat))
  print("- E[p_L]=E[p_sp1*p_sp2]: %.4f" % (E_p_L/repeat))
  print("- Log E[p_sp]: %.4f / %.4f" % (log2(E_p_sp1) - log2(repeat), log2(E_p_sp2) - log2(repeat)))
  print("- Log E[p_L]: %.4f" % (log2(E_p_L) - log2(repeat)))
  print()

  Exp_1 = E_p_L/repeat

  print("Experiment 2: e1, e2 independently sampled from U([-%d, %d]^%d)" % (l1, l1, r))
  E_p_sp1 = 0.0; E_p_sp2 = 0.0
  E_p_L = 0.0

  for _ in range(repeat):
    v1 = np.zeros(r)
    v2 = np.zeros(r)
    for i in range(r):
      v1[i] = np.random.uniform(-l1, l1)
      v2[i] = np.random.uniform(-l1, l1)

    p_rep1 = scaler
    p_rep2 = scaler
    for i in range(r):
      p_rep1 *= 1.0 - abs(v1[i]) / (2*l2)
      p_rep2 *= 1.0 - abs(v2[i]) / (2*l2)

    p_sp1 = 1-(1 - p_rep1)**R
    p_sp2 = 1-(1 - p_rep2)**R

    E_p_sp1 += p_sp1
    E_p_sp2 += p_sp2
    E_p_L += p_sp1*p_sp2
  
  print("- E[p_sp]: %.4f / %.4f" % (E_p_sp1/repeat, E_p_sp2/repeat))
  print("- E[p_L]=E[p_sp1*p_sp2]: %.4f" % (E_p_L/repeat))
  print("- Log E[p_sp]: %.4f / %.4f" % (log2(E_p_sp1) - log2(repeat), log2(E_p_sp2) - log2(repeat)))
  print("- Log E[p_L]: %.4f" % (log2(E_p_L) - log2(repeat)))

  Exp_2 = E_p_L/repeat

  return Exp_1, Exp_2

def dep_test_lv2(l1, l2, l3, r, error_param, error_unif=False, C=10, repeat=10000, scaler=1):
  assert l1 <= l2, "l should increase"
  assert l2 <= l3, "l should increase"

  E_p_rep = (p_adm_unif(l2, 2*l3))**r
  
  print("******* Goal of this experiment **********")
  print("Define p(e) := Product of [1.0 - |e_i| / (2*l3)]")
  print("       & p_sp(e, R) := 1 - (1 - p(e))^R")
  print("Compare E[p_L] = E[p_sp(e_21)*p_sp(e_22)*p_sp(e_23)*p_sp(e_24)]")
  print("       while assuming several dependence between e_2*")
  print()
  print("* Side Info")
  R = int(C/E_p_rep)
  print("Log E[p_rep] = %.3f, R = %d/p_rep = %d" % (log2(E_p_rep), C, R))
  if scaler != 1:
    R = int(R / scaler)
    E_p_rep *= scaler
    print("(R is scaled by %.4f -> R = %d)" % (scaler, R))
  print("******************************")
  print()
  
  print("* Experiment 1: e_0: sampled from", end = "")
  if error_unif: 
    print(" U[%.2f, %.2f])^%d " % (-error_param, error_param, r))
  else:
    print(" G(%.2f)^%d" % (error_param, r))
  print("                e_11: sampled from U([-%d, %d]^%d)" % (l1, l1, r))  
  print("                e_12 = e_11 + e_0")
  print("                e_21, e_23: (independently) sampled from U([-%d, %d]^%d)" % (l2, l2, r))
  print("                e_22 = e_21 + e_11 and e_24 = e_23 + e_12")

  E_p_sp1 = 0.0; E_p_sp2 = 0.0; E_p_sp3 = 0.0; E_p_sp4 = 0.0
  E_p_L = 0.0

  for _ in range(repeat):
    e_0 = np.zeros(r)
    if error_unif:
      for i in range(r):
        e_0[i] = np.random.uniform(-error_param, error_param)
    else:
      for i in range(r):
        e_0[i] = np.random.normal(scale=error_param)

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
      p_rep1 *= 1.0 - abs(e_21[i]) / (2*l3)
      p_rep2 *= 1.0 - abs(e_22[i]) / (2*l3)
      p_rep3 *= 1.0 - abs(e_23[i]) / (2*l3)
      p_rep4 *= 1.0 - abs(e_24[i]) / (2*l3)

    p_sp1 = 1-(1 - p_rep1)**R; p_sp2 = 1-(1 - p_rep2)**R; p_sp3 = 1-(1 - p_rep3)**R; p_sp4 = 1-(1 - p_rep4)**R

    E_p_sp1 += p_sp1; E_p_sp2 += p_sp2; E_p_sp3 += p_sp3; E_p_sp4 += p_sp4     
    E_p_L += p_sp1*p_sp2*p_sp3*p_sp4

  print("- E[p_sp]: %.4f / %.4f / %.4f / %.4f" % (E_p_sp1/repeat, E_p_sp2/repeat, E_p_sp3/repeat, E_p_sp4/repeat))
  print("- E[p_L]: %.4f" % (E_p_L/repeat))
  print("- Log E[p_sp]: %.4f / %.4f / %.4f / %.4f" % (log2(E_p_sp1) - log2(repeat), log2(E_p_sp2) - log2(repeat), log2(E_p_sp3) - log2(repeat), log2(E_p_sp4) - log2(repeat)))
  print("- Log E[p_L]: %.4f" % (log2(E_p_L) - log2(repeat)))
  print()

  Exp_1 = E_p_L/repeat
  
  print("* Experiment 2: Same with 1, except e_11 & e_12 are independent")
  print("                i.e., e_11, e_12: independently sampled from U([-%d, %d]^%d)" % (l1, l1, r))  
  E_p_sp1 = 0.0; E_p_sp2 = 0.0; E_p_sp3 = 0.0; E_p_sp4 = 0.0
  E_p_L = 0.0

  for _ in range(repeat):
    e_11 = np.zeros(r); e_12 = np.zeros(r)
    for i in range(r):
      e_11[i] = np.random.uniform(-l1, l1)
      e_12[i] = np.random.uniform(-l1, l1)

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
      p_rep1 *= 1.0 - abs(e_21[i]) / (2*l3)
      p_rep2 *= 1.0 - abs(e_22[i]) / (2*l3)
      p_rep3 *= 1.0 - abs(e_23[i]) / (2*l3)
      p_rep4 *= 1.0 - abs(e_24[i]) / (2*l3)

    p_sp1 = 1-(1 - p_rep1)**R; p_sp2 = 1-(1 - p_rep2)**R; p_sp3 = 1-(1 - p_rep3)**R; p_sp4 = 1-(1 - p_rep4)**R

    E_p_sp1 += p_sp1; E_p_sp2 += p_sp2; E_p_sp3 += p_sp3; E_p_sp4 += p_sp4     
    E_p_L += p_sp1*p_sp2*p_sp3*p_sp4
  
  print("- E[p_sp]: %.4f / %.4f / %.4f / %.4f" % (E_p_sp1/repeat, E_p_sp2/repeat, E_p_sp3/repeat, E_p_sp4/repeat))
  print("- E[p_L]: %.4f" % (E_p_L/repeat))
  print("- Log E[p_sp]: %.4f / %.4f / %.4f / %.4f" % (log2(E_p_sp1) - log2(repeat), log2(E_p_sp2) - log2(repeat), log2(E_p_sp3) - log2(repeat), log2(E_p_sp4) - log2(repeat)))
  print("- Log E[p_L]: %.4f" % (log2(E_p_L) - log2(repeat)))
  print()

  Exp_2 = E_p_L/repeat

  print("* Experiment 3: All e_21 ... e_24 independently sampled from U([-%d, %d]^%d)" % (l2, l2, r))  
  E_p_sp1 = 0.0; E_p_sp2 = 0.0; E_p_sp3 = 0.0; E_p_sp4 = 0.0
  E_p_L = 0.0

  for _ in range(repeat):

    e_21 = np.zeros(r); e_22 = np.zeros(r); e_23 = np.zeros(r); e_24 = np.zeros(r)
    for i in range(r):
      e_21[i] = np.random.uniform(-l2, l2)
      e_22[i] = np.random.uniform(-l2, l2)
      e_23[i] = np.random.uniform(-l2, l2)
      e_24[i] = np.random.uniform(-l2, l2)
        
    p_rep1 = scaler; p_rep2 = scaler; p_rep3 = scaler; p_rep4 = scaler
    for i in range(r):
      p_rep1 *= 1.0 - abs(e_21[i]) / (2*l3)
      p_rep2 *= 1.0 - abs(e_22[i]) / (2*l3)
      p_rep3 *= 1.0 - abs(e_23[i]) / (2*l3)
      p_rep4 *= 1.0 - abs(e_24[i]) / (2*l3)

    p_sp1 = 1-(1 - p_rep1)**R; p_sp2 = 1-(1 - p_rep2)**R; p_sp3 = 1-(1 - p_rep3)**R; p_sp4 = 1-(1 - p_rep4)**R

    E_p_sp1 += p_sp1; E_p_sp2 += p_sp2; E_p_sp3 += p_sp3; E_p_sp4 += p_sp4     
    E_p_L += p_sp1*p_sp2*p_sp3*p_sp4
  
  print("- E[p_sp]: %.4f / %.4f / %.4f / %.4f" % (E_p_sp1/repeat, E_p_sp2/repeat, E_p_sp3/repeat, E_p_sp4/repeat))
  print("- E[p_L]: %.4f" % (E_p_L/repeat))
  print("- Log E[p_sp]: %.4f / %.4f / %.4f / %.4f" % (log2(E_p_sp1) - log2(repeat), log2(E_p_sp2) - log2(repeat), log2(E_p_sp3) - log2(repeat), log2(E_p_sp4) - log2(repeat)))
  print("- Log E[p_L]: %.4f" % (log2(E_p_L) - log2(repeat)))
  print()
  
  Exp_3 = E_p_L/repeat

  return Exp_1, Exp_2, Exp_3