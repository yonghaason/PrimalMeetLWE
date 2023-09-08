from estimator import *
param_specified = {}
filename = False

logn=15; logq=699; h=128
# logn=15; logq=768; h=192
# logn=16; logq=1450; h=64
# logn=16; logq=1553; h=192

param = [2**logn, 2**logq, 'gaussian', 3.2, h, 2**logn]

for lv in range(1, 4):
  filename = "0908_log/" + str(logn) + "_" + str(logq) + "_" + str(h) + "_" + str(lv) + "_MATZOV.txt"
  _ = primal_may(param, t=lv, param_specified=param_specified, filename=filename)

# param_specified = {'d':11008, 'm':111873, 'beta': 295}
# _ = primal_may(param, t=3, param_specified=param_specified, filename=filename)
