from estimator import *

# n=2**15;q=2**768;stddev=3.19;is_sec_bal=False;w=192

# zeta=4864;d=54017;beta=312;nu=41.674
# helv22 = cost_beta(beta, zeta, n, q, d, stddev, is_sec_bal, w, nu, lv=2)

n=2**16;q=2**1450;stddev=3.19;is_sec_bal=False;w=64

zeta=28672;d=70657;beta=176;nu=102.08
lv3 = cost_beta(beta, zeta, n, q, d, stddev, is_sec_bal, w, nu, lv=3)

# zeta=5760;d=53633;beta=298;nu=41.674
# helv2 = cost_beta(beta, zeta, n, q, d, stddev, is_sec_bal, w, nu, lv=2)

# zeta=5760;d=52609;beta=298;nu=41.674
# helv3 = cost_beta(beta, zeta, n, q, d, stddev, is_sec_bal, w, nu, lv=3)

# n=508;q=2048;stddev=0.82;is_sec_bal=False;w=254

# zeta=130;d=847;beta=405;nu=1.155
# ntrulv1=cost_beta(beta, zeta, n, q, d, stddev, is_sec_bal, w, nu, lv=1)   