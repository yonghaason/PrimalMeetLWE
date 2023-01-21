import numpy as np
import sys
from math import log2, comb, pi, e, erf, sqrt, exp, ceil, log

def Log2(a):
    try:
        return round(log2(a), 2)
    except ValueError:
        return -np.inf

def prettyprint(dic, filename = False):
    
    if filename is not False:
        f = open(filename, 'a')
        sys.stdout = f

    print('* Final Cost:', dic['cost'])
    print('* Lattice:', dic['lat'])
    print('* Guess: {', end='')
    for keyy in dic['guess']:
        if keyy != 'lsh_info':
            print('\'', keyy, '\': ', dic['guess'][keyy], sep = '', end=', ')
    print('}')
    print('---- Further details on guess')
    if 'lsh_info' in dic['guess']:
        for i in range(len(dic['guess']['lsh_info'])):
            print('* lsh info of Lv', i+1, '->', i, ':', dic['guess']['lsh_info'][i])

    if filename is not False:
        f.close()

def my_binary_search(left, right, func, cur_depth, key='cost', 
                    max_depth = 5, val_l = None, val_m = None, val_r = None):
    ''' 
    Main idea is to shrink the input range to half length
    by comparing function values on "middle_left, middle, middle_right".
    I strongly believe there would be more elegant way to implement this idea than below, 
    but also believe that this logic indeed finds global minima.
    '''
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

    if val_min == val_ml[key]:
        if cur_depth == max_depth:
            return val_ml
        else:
            return my_binary_search(left, middle, func, cur_depth+1, key, val_l = val_l, val_m = val_ml, val_r = val_m)
    elif val_min == val_m[key]:
        if cur_depth == max_depth:
            return val_m
        else:
            return my_binary_search(middle_left, middle_right, func, cur_depth+1, key, val_l = val_ml, val_m = val_m, val_r = val_mr)
    else:
        if cur_depth == max_depth:
            return val_mr
        else:
            return my_binary_search(middle, right, func, cur_depth+1, key, val_l = val_m, val_m = val_mr, val_r = val_r)

##################################################
#                Lattices
##################################################

# Taken from "Lattice-estimator"
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

# Taken from "Lattice-estimator"
def log_BKZ_cost(d, beta, C=5.46):
    svp_calls = C * max(d - beta, 1)
    beta_ = beta - d4f(beta)
    # Old one from ADPS16
    # gate_count = C * 2 ** (0.2988026130564745 * beta_ + 26.011121212891872)
    # New one from MATZOV22
    gate_count = C * 2 ** (0.29613500308205365 * beta_ + 20.387885985467914)
    return log2(LLL(d) + svp_calls * gate_count)

# Taken from "Lattice-estimator"
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

# Adapted from "Lattice-estimator"
def ext_GSA(m, q, d, d1, beta, nu = 1):    
    Rfactors = [q] * (m - d)

    d2 = d - d1
    log_vol = float(log2(q) * (d2 - 1) + log2(nu) * d1)
    delta = delta_0(beta)
    r_log = [(d - 1 - 2 * i) * float(log2(delta)) + log_vol / d for i in range(d)]
    for r_ in r_log:
        try:
            Rfactors.append(2**r_)
        except OverflowError:
            Rfactors.append(np.inf)

    return Rfactors

def prob_np(q_list, sigma):
    pr_np = 1
    for q in q_list:
        bnd = q / (2 * sigma * sqrt(2))
        pr_np *= erf(bnd)
    return pr_np

# Probability of n-dim, w pm 1 => d-dim, w0 pm 1
def prob_hw_half(n, d, w, w0):    
    return comb(n-2*w,d-2*w0) * comb(w, w0)**2 / comb(n,d)

# Probability of n-dim, w nonzero => zeta-dim, less than w0 nonzero
def prob_hw(n, zeta, w, w0):
    cases = 0
    for i in range(w0+1):
        cases += comb(n-w, zeta-i) * comb(w, i)
    return cases / comb(n, zeta)

def prob_hw_exact(n, zeta, w, w0):
    return comb(n-w, zeta-w0) * comb(w, w0) / comb(n, zeta)

##################################################
#                Combinatorics
##################################################

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
            res*=comb(n_, c[i])        
            n_ = n_ - c[i]
        return res
    except ValueError:
        return np.inf

def easy_multinom(n, a, b):
    return multinom(n, [a, b, n - a - b])

def ambiguity_half(n, w, w_sp):
    w_even = w
    if w % 2 == 1:
        w_even += 1
    if n < 2*w_even:
        n += 1
    eps = w_sp - w_even // 2
    try:
        return comb(w_even, w_even//2)**2 * multinom(n - 2*w_even, [eps, eps, n - 2*w_even - 2*eps])
    except ValueError:
        return 0

def ambiguity(n, h, w):
    ell = w - (h+1)//2
    return comb(h, w - ell) * comb(n - h, ell) * 2**ell

def prob_admissible(b, k, GSnorm, unif=True):
    
    p = 1
    for i in range(k):
        p *= prob_admissible_axis(b, GSnorm[-i-1], unif=unif)

    return p

def prob_admissible_axis(b, l, unif=True):
    '''
    Consider a random element x in a length l box (NOT a torus!).
    Output the probability that x + e also in the same box,
    when e is b-bounded uniform error (or Gaussian error of std.dev b)
    '''
    if unif:
        if l > 2*b:
            return 1 - b/(2*l)
        elif l > b:
            return 1/4 * (1 + l/b)
        else:
            return l/(2*b) # just in case..
    else:
        x = l/(sqrt(2) * b)
        try:
            prob = erf(x) + (exp(-x**2) - 1)/(x * sqrt(pi))
            if (prob > 1):
                print("prob > 1 ! with prob =", prob)
                return 1
            return prob
        except OverflowError:
            return 1

def amplify(target_success_probability, success_probability):
    try:
        # target_success_probability = 1 - (1-success_probability)^trials
        return ceil(log2(1 - target_success_probability) / log2(1 - success_probability))
    except ValueError:
        return np.inf