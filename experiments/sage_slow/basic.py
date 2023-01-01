from sage.all import shuffle, randint, ceil, next_prime, log, cputime, set_random_seed, sqrt
from copy import copy
from sage.all import GF, ZZ
from sage.all import random_matrix, random_vector, zero_vector, vector, matrix, identity_matrix
from sage.rings.all import QQ, RR, ZZ, RealField, PowerSeriesRing, RDF
from sage.rings.infinity import PlusInfinity, Infinity
from sage.structure.element import parent
from sage.symbolic.all import pi, e
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler \
    as DiscreteGaussian

from fpylll import IntegerMatrix, FPLLL, GSO
from math import log2, comb, pi, e, erf, sqrt, exp, ceil, log

def modq(v, q, bal = True):
    v %= q
    if bal and v > q / 2:
        v -= q
    return v

def sparse_secret_list(M, n, w):
    m = M.nrows()
    if n < w:
        return []
    elif w == 0:
        return [[[0] * n, vector(ZZ, [0] * m)]]
    else:
        S = []
        S1 = sparse_secret_list(M, n-1, w-1)
        S2 = sparse_secret_list(M, n-1, w)
        for v in S1:
            S.append([v[0] + [1], v[1] + M.column(n-1)])
            S.append([v[0] + [-1], v[1] - M.column(n-1)])
        for v in S2:
            S.append([v[0] + [0], v[1]])
        return S

def gen_instance(n, stddev, q, h, m=None):
    if m is None:
        m = n

    A = matrix(ZZ, m, n)
    for i in range(m):
        for j in range(n):
            A[i,j] = randint(0, q-1)
    
    S = [-1, 1]
    s = [S[randint(0, 1)] for i in range(h)]
    s += [0 for _ in range(n-h)]
    shuffle(s)
    s_ = s + [-1]
    s = vector(ZZ, s)
    b = A*s

    D = DiscreteGaussian(stddev)

    for i in range(m):
        b[i] += D()
        b[i] = modq(b[i], q)

    M = matrix(ZZ, m, n+1)
    for i in range(m):
        row = list(A[i]) + [b[i]]
        M[i] = vector(row)

    return M, vector(ZZ, s_)

def ambiguity(n, h, w):
    ell = w - (h+1)//2
    return comb(h, w - ell) * comb(n - h, ell) * 2**ell

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

def constraint_area(R, b, l, GSnorm, unif=True):
    '''
    For a ternary s of weight w0 satisfying ||Bs||_inf < b,
    we have R many ternary pairs (s1, s2) of weight w1 s.t. s = s1 - s2.
    We compute a box A of length l, 
    where there is at least good pair (s1, s2) s.t. both Bs1 & Bs2 lie on A.
    '''
    
    r = 1
    p_adm = 1
    vol_ratio = 1

    while True:
        cur_p_adm = 1
        if GSnorm[-r] < l:
            cur_p_adm = prob_admissible_axis(b, GSnorm[-r], unif = unif)
        else:
            cur_ratio = GSnorm[-r] / l
            cur_p_adm = prob_admissible_axis(b, l, unif=unif) / cur_ratio

        if r == len(GSnorm):
            break

        inv = 1.0 / (p_adm * cur_p_adm)
        if inv > R:
            r -= 1
            break
        else:
            p_adm *= cur_p_adm
            r += 1

    return r, p_adm

def build_list(M, w, r, l, q):
    ''' 
    Output all vectors s such that
    1) weight w
    2) pi_r(Mv) in box of length l
    '''
    m = M.nrows()
    n = M.ncols()

    # S = {(v, Mv): HW(v) = w}
    S = sparse_secret_list(M, n, w)

    L = []    
    for pair in S:
        v = pair[0]; Mv = pair[1]
        Mv = vector(ZZ, [modq(Mv[i], q) for i in range(m)])

        proj = vector(ZZ, Mv[-r:])
        if proj.norm(Infinity) <= l:
            L.append([vector(ZZ, v), Mv])
            
    return L

def near_collision(L, stddev, q, h):
    sol = []
    N = len(L)
    for i in range(N):
        for j in range(i+1, N):
            s_ = L[i][0] - L[j][0]
            if s_.hamming_weight() != h:
                continue
            if s_[-1] == 1:
                continue
            # Ms_ = M*s_
            Ms_ = L[i][1] - L[j][1]
            Ms_ = vector(ZZ, [modq(Ms_[i], q) for i in range(len(Ms_))])

            if Ms_.norm(Infinity) < 6*stddev:
                sol.append(s_)
    return sol