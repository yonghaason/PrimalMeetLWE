from estimator import *
param_specified = {'zeta': 6144, 'd': 51457, 'beta': 295}
param = [2**15, 2**768, 'gaussian', 3.19, 192, 2**15]
_ = primal_may(param, lv=2, param_specified=param_specified)