from estimator import *
param_specified = {}
filename = False

n=509; q=2048; h=254
# n=512; q=12289; h=154
# n=512; q=8383489; h=342
# n=653; q=4621; h=288
# n=659; q=2048; h=76
# n=677; q=2048; h=254
# n=761; q=2048; h=84
# n=761; q=4591; h=286
# n=821; q=4096; h=510
# n=857; q=5167; h=322
# n=1087; q=2048; h=126
# n=1499; q=2048; h=158

param = [n, q, 'unif', 1, h, n]

for lv in range(1, 4):
  filename = "0908_log/" + str(n) + "_" + str(q) + "_" + str(h) + "_" + str(lv) + "_MATZOV.txt" 
  _ = primal_may(param, t=lv, param_specified=param_specified, filename=filename)

# param = [n, q, 'unif', 1, h, 100]
# param_specified = {'m': 786, 'd': 225, 'beta': 529, 'w_g': 106}
# _ = primal_may(param, t=2, param_specified=param_specified, filename=filename)
