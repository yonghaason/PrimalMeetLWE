from estimator_euro import *
import sys

param = None
if len(sys.argv) == 6:
    logn = int(sys.argv[1])
    logq = int(sys.argv[2])
    q = 2**logq
    n = 2**logn
    param = [n, q, float(sys.argv[3]), int(sys.argv[4]), n]
    lv = int(sys.argv[5])
    filename = 'output/n_' + str(n) + '_logq_' + str(logq) + '_hw_' + sys.argv[4] + '_lv_' + str(lv) + '.txt'
else:
    print("Caution: Set Default param (HE n = 2**15)")
    param = HE_param(2**15)

#print("Input Params")
#print("(n, q, stddev, w, m) =", param)
#print("----------------------")
primal_may(param, lv=lv, filename=filename)

print('----------------- lv 1')
best1 = primal_may(param, lv=1)
prettyprint(best1)

print('----------------- lv 2')
best2 = primal_may(param, lv=2)
prettyprint(best2)

print('----------------- lv 3')
best3 = primal_may(param, lv=3)
prettyprint(best3)

