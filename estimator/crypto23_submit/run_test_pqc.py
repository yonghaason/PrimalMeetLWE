from estimator import *
import sys
import os.path
            
param = None
filename = False
if len(sys.argv) >= 7:
    lv = int(sys.argv[1])
    n = int(sys.argv[2])
    q = int(sys.argv[3])
    error_type = sys.argv[4]
    error_param = float(sys.argv[5])
    secret_weight_param = int(sys.argv[6])
    # is_sec_bal = (sys.argv[6].lower() == 'true')
    is_sec_bal = False
    param = [n, q, error_type, error_param, secret_weight_param, n, is_sec_bal]
    if len(sys.argv) == 8:
        if sys.argv[7] == 'save':
            filename = 'log/' + str(n) + '_' + str(q) + '_' + str(secret_weight_param)
            if is_sec_bal:
                filename += '_b'
            filename +=  '_lv' + str(lv) + '.txt'
            if os.path.isfile(filename):
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                print("         ", filename, "already exists, runs as print mode")
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
                filename = False
            else:
                print("Output will be saved in", filename)
    primal_may(param, lv=lv, filename=filename)
else:
    print("Check whether arguments are complete")
    print("- (lv) (n) (q) (error_type) (error_param) (secret_weight_param)")