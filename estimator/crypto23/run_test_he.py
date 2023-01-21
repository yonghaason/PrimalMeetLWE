from estimate import *
import sys
import os.path
            
param = None
filename = False
if len(sys.argv) >= 8:
    lv = int(sys.argv[1])
    logn = int(sys.argv[2])
    logq = int(sys.argv[3])
    error_type = sys.argv[4]
    error_param = float(sys.argv[5])
    is_sec_bal = (sys.argv[6].lower() == 'true')
    secret_weight_param = int(sys.argv[7])
    q = 2**logq
    n = 2**logn
    param = [n, q, error_type, error_param, is_sec_bal, secret_weight_param, n]
    if len(sys.argv) == 9:
        if sys.argv[8] == 'save':
            filename = 'log/' + str(logn) + '_' + str(logq) + '_' + str(secret_weight_param)
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
    print("- (lv) (logn) (logq) (error_type) (error_param) (is_secret_balanced) (secret_weight_param)") 