import numpy as np
import sys
import operator as op
from functools import reduce
from math import log2, pi, e, erf, sqrt, exp, ceil, log, floor

def Log2(a):
    try:
        return round(log2(a), 2)
    except ValueError:
        return -np.inf

def HE_param(n):
    stddev = 3.19; w = 64; 
    if n == 2048:
        q = 2**45
    elif n == 4096:
        q = 2**82
    elif n == 8192:
        q = 2**158
    elif n == 2**15:
        q = 2**768
    elif n == 2**16:
        q = 2**1553; w = 192
    else:
        print("Not in preset params")

    return [n, q, 'gaussian', stddev, True, w, n]

# NTRU Round3
def NTRU_param(n):
    if n == 508:
        q = 2048; w = 254
    elif n == 676:
        q = 2048; w = 254
    elif n == 820:
        q = 4096; w = 510
    elif n == 652:
        q = 4621; w = 288
    elif n == 760:
        q = 4591; w = 286
    elif n == 856:
        q = 5167; w = 322
    else:
        print("Not in preset params")
    return [n, q, 'unif', 1, True, w, n]

def prettyprint(dic, filename = False):
    
    if filename is not False:
        f = open(filename, 'a')
        sys.stdout = f

    print('* Final Cost:', dic['cost'])
    for keyy in dic:
        if keyy != 'cost' and keyy != 'guess_details':
            print('   - \'', keyy, '\': ', dic[keyy], sep = '', end='\n')
    
    guess = dic['guess_details']
    # guess = dic['guess_details'][-1]
    print('Guess details:')
    for keyy in ['cost', 'log_p_suc']:
        print('   - \'', keyy, '\': ', guess[keyy], sep = '', end='\n')
    for keyy in guess:
        if keyy not in ['lsh_stats', 'cost', 'log_p_suc', 'top_reach_fail'] and len(guess[keyy]) != 0:
            print('   - \'', keyy, '\': ', guess[keyy], sep = '', end='\n')
    if 'top_reach_fail' in guess:
        print('   - \'top_reach_fail\': ', guess['top_reach_fail'], sep = '', end='\n')
    
    if filename is not False:
        f.close()

def my_binary_search(left, right, func, cur_depth, key='cost', 
                    max_depth = 5, val_l = None, val_m = None, val_r = None):
    ''' 
    Main idea is to shrink the input range to half length
    by comparing function values on "middle_left, middle, middle_right".
    There would be more elegant way to implement this idea than below, 
    but this logic indeed finds global minima, if exists.
    '''
    if left == right:
        return func(left)
        
    middle = (left + right) // 2
    if middle == left: # right = left + 1 case
        if val_l is None:
            val_l = func(left)
        if val_r is None:
            val_r = func(left)
        val_min = min(val_l[key], val_r[key])
        if val_min == val_l[key]:
            return val_l
        else:
            return val_r
    else: 
        if val_m is None:
            val_m = func(middle)

    middle_left = (left + middle) // 2
    middle_right = (middle + right) // 2
    val_ml = None
    val_mr = None
    if middle_left == left:
        if val_l is None:
            val_ml = func(left)
        else:
            val_ml = val_l
    else:
        val_ml = func(middle_left)
    if middle_right == middle: 
        if val_m is None:
            val_mr = func(middle)
        else:
            val_mr = val_m
    else: 
        val_mr = func(middle_right)

    # Here is the main logic (the above are for extreme cases)
    vals = [val_ml[key], val_m[key], val_mr[key]]
    vals.sort()
    val_min = vals[0]

    if val_min == val_mr[key]:
        if cur_depth == max_depth:
            return val_mr
        else:
            return my_binary_search(middle, right, func, cur_depth+1, key, 
                    val_l = val_m, val_m = val_mr, val_r = val_r, max_depth = max_depth)
    elif val_min == val_m[key]:
        if cur_depth == max_depth:
            return val_m
        else:
            return my_binary_search(middle_left, middle_right, func, cur_depth+1, key, 
                    val_l = val_ml, val_m = val_m, val_r = val_mr, max_depth = max_depth)
    else:
        if cur_depth == max_depth:
            return val_ml
        else:
            return my_binary_search(left, middle, func, cur_depth+1, key, 
                    val_l = val_l, val_m = val_ml, val_r = val_m, max_depth = max_depth)

##################################################
#                Lattices
#   Some functions are taken, or adapted from "lattice-estimator"
##################################################

def delta_0(beta):
    small = (
        (2, 1.02190),  # noqa
        (5, 1.01862),  # noqa
        (10, 1.01616),
        (15, 1.01485),
        (20, 1.01420),
        (25, 1.01342),
        (28, 1.01331),
        (40, 1.01295),
    )

    if beta <= 2:
        return float(1.0219)
    elif beta < 40:
        for i in range(1, len(small)):
            if small[i][0] > beta:
                return float(small[i - 1][1])
    elif beta == 40:
        return float(small[-1][1])
    else:
        return float(beta / (2 * pi * e) * (pi * beta) ** (1 / beta)) ** (1 / (2 * (beta - 1)))


def log_BKZ_cost(d, beta, model="MATZOV22"):    
    beta_ = beta - d4f(beta)
    if model == "ADPS16":
        return 0.292 * beta_
    if model == "MATZOV22":
        svp_calls = 5.46 * max(d - beta, 1)
        gate_count = 5.46 * 2 ** (0.29613500308205365 * beta_ + 20.387885985467914)
        return log2(LLL(d) + svp_calls * gate_count)

def d4f(beta):
    """
    Dimensions "for free" following [EC:Ducas18]_.
    :param beta: Block size ≥ 2.
    If β' is output by this function then sieving is expected to be required up to dimension β-β'.
    EXAMPLE::
        >>> from estimator.reduction import RC
        >>> RC.Kyber.d4f(500)
        42.597...
    """
    return max(float(beta * log2(4 / 3.0) / log2(beta / (2 * pi * e))), 0.0)

def LLL(d, B=None):
    """
    Runtime estimation for LLL algorithm based on [AC:CheNgu11]_.
    :param d: Lattice dimension.
    :param B: Bit-size of entries.
    """
    if B:
        return d ** 3 * B ** 2
    else:
        return d ** 3  # ignoring B for backward compatibility

def GSA(q, d, d1, beta, nu = 1):    
    GSnorm = []
    d2 = d - d1
    log_vol = float(log2(q) * (d2 - 1) + log2(nu) * d1)
    delta = delta_0(beta)
    r_log = [(d - 1 - 2 * i) * float(log2(delta)) + log_vol / d for i in range(d)]
    for r_ in r_log:
        try:
            GSnorm.append(2**r_)
        except OverflowError:
            GSnorm.append(np.inf)
    return GSnorm

def GSA_mod(q, d, d1, beta, nu = 1):    
    GSnorm = []
    d2 = d - d1
    log_vol = float(log2(q) * d2 + log2(nu) * d1)
    delta = delta_0(beta)
    first_R = (d - 1) * float(log2(delta)) + log_vol / d
    while first_R > log2(q):
        d -= 1
        log_vol -= log2(q)
        GSnorm.append(q)
        first_R = (d - 1) * float(log2(delta)) + log_vol / d
    r_log = [(d - 1 - 2*i) * float(log2(delta)) + log_vol / d for i in range(d)]
    for r_ in r_log:
        try:
            GSnorm.append(2**r_)
        except OverflowError:
            GSnorm.append(np.inf)
    return GSnorm

def prob_np(GSnorm, stddev):
    pr_np = 1
    for q in GSnorm:
        bnd = q / (2 * stddev * sqrt(2))
        pr_np *= erf(bnd)
    return pr_np

def babai_cost(d):
    return d**2
    # return d

##################################################
#            Maths & Combinatorics
##################################################

def comb(n, c):
    if n <= 0:
        return 0
    r = min(c, n-c)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer // denom

#This function is taken from https://github.com/ElenaKirshanova/ntru_with_lsh
def multinom(n, c):
    """
    Product of binimial coefficients:
    (n choose c[0]) * (n-c[0] choose c[1])*...*
    """
    assert sum(c) == n, 'bad input to multinom!'
    try:
        res = 1
        n_ = n
        for i in range(len(c)):
            res *= comb(n_, c[i])        
            n_ = n_ - c[i]
        return res
    except ValueError:
        return np.inf

def easy_multinom(n, a, b):
    return multinom(n, [a, b, n - a - b])

def num_ternary_secret(n, w, is_sec_bal=False):
    if is_sec_bal:
        return easy_multinom(n, w, w)
    else:
        return 2**w * comb(n, w)

def prob_hw(n, zeta, w, w0, is_sec_bal=False):
    '''
    If balanced secret: 
        prob of n-dim, w pm 1 (2w nonzero) => zeta-dim, w0 pm 1 (2w0 nonzero)
    else: 
        prob of n-dim, w nonzero => zeta-dim, exactly w0 nonzero
    '''
    if is_sec_bal:
        return comb(n-2*w, zeta-2*w0) * comb(w, w0)**2 / comb(n, zeta)
    else:
        return comb(n-w, zeta-w0) * comb(w, w0) / comb(n, zeta)

def probs_hw_precompute(n, zeta, w, is_sec_bal=False, bound = 31):
    numer = comb(n, zeta)
    if is_sec_bal:
        return [comb(n-2*w, zeta-2*i) * comb(w, i)**2 / numer for i in range(bound)]
    else:
        return [comb(n-w, zeta-i) * comb(w, i) / numer for i in range(bound)]
        

def ambiguity(n, h, w, is_sec_bal=False):
    '''
    The number of representations s = s_1 - s_2
    h: weight parameter of s
    w: weight parameter of s_1 and s_2
    if is_sec_bal:
        s, s_1, s_2 have balanced number of +1 and -1
        exactly, s has h +1s and h -1s ..
    else:
        s has total h +1 or -1.
    '''
    assert h % 2 == 0, 'weight should be even'
    eps = w - h//2
    if is_sec_bal:
        return comb(h, h//2)**2 * easy_multinom(n - 2*h, eps, eps)
    else:
        return comb(h, w - eps) * comb(n - h, eps) * 2**eps

def prob_admissible_gaussian(b, stddev):
    # Prob[x + y in [0, b]], where x <-[0, b] and y <- Gaussian(stddev)
    x = b/(sqrt(2) * stddev)
    try:
        return erf(x) + (exp(-x**2) - 1)/(x * sqrt(pi))
    except OverflowError:
        return 1

def prob_admissible_uniform(b, e):
    # Prob[x + y in [0, b]], where x <-[0, b] and y <- [-e, e]
    if e <= b:
        return 1 - e/(2*b)
    else:
        return b/(2*e)

p_jensen_hardcode_gauss = [0.994, 0.976, 0.953, 0.923, 0.892, 0.86, 0.823, 0.797, 0.763, 0.741, 0.702, 0.677, 0.646, 0.617, 0.595]
p_jensen_hardcode_unif = [0.964, 0.891, 0.815, 0.733, 0.655, 0.562, 0.523, 0.463, 0.421, 0.373, 0.327, 0.303]