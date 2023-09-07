from estimator import *
param_specified = {}
filename = False

n = 509; q = 2048; h = 254
param = [n, q, 'unif', 1, h, n]

# for lv in range(1, 4):
#   filename = str(n) + "_" + str(q) + "_" + str(h) + "_" + str(lv) + "_MATZOV.txt" 
#   _ = primal_may(param, t=lv, param_specified=param_specified, filename=filename)

param = [n, q, 'unif', 1, h, 100]
param_specified = {'m': 786, 'd': 225, 'beta': 529, 'w_g': 106}
_ = primal_may(param, t=2, param_specified=param_specified, filename=filename)
