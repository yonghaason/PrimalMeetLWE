from utils import *
from functools import partial
import copy

"""

* Requires Python >= 3.8,
because "math" module supports comb function after this version.
You may implement and replace comb function, to run this code with lower version of python.

* How to use?
In python prompt:

> from estimator import *
> param = [(n), (q), (stddev), (hw), (m)] # Insert the parameters
> primal_may(param, (lv)) # Insert the level

* It takes LONG LONG time (over several hours!) for high level, or large LWE parameters. If the reviewer is curious about the detailed parameters (that provides the numbers in Table 1 of paper), we recommend to see the text files in `logs' directory.

* About the consistency with `lattice-estimator:

Our level-0 attack strategy is covered by `primal_hybrid' with options `mitm=False, babai=True' in lattice-estimator.
We observe that both scripts output similar numbers for this strategy.
Meanwhile, our level-1 attack strategy is almost the same as the primal hybrid-MitM attack [HG07].
This attack is covered by 'primal_hybrid' with options 'mitm=True, babai=True' in 'lattice-estimator'.
However, we found that the computation is over-simplified from the original description in [HG07],
and we believe our script more accurately computes the attack complexity.
This results in quite different estimation of this attack for PQC parameters.
However, it has no actual impact since the both resulting numbers are non-competitive compared to other attacks.
For FHE parameters, two estimations are quite similar
and there is no notable problem from lattice-estimator's over-simplified estimation.

"""

# HE
def HE_param(n, filename=False):

    if filename is not False:
        f = open(filename, 'a')
        sys.stdout = f

    stddev = 3.19; w = 64; 
    if n == 2048:
        q = 2**45
    elif n == 4096:
        q = 2**82
    elif n == 8192:
        q = 2**158
    elif n == 2**15:
        q = 2**768; w = 64
    elif n == 2**16:
        q = 2**1553; w = 192
    else:
        print("Not in preset params")
    m_max = n

    print("Sparse-ternary HE Parameter setting")
    print("n = %d, q = %d, stddev = %.2f, w = %d, m = %d" % (n, q, stddev, w, m_max))

    if filename is not False:
        f.close()

    return [n, q, stddev, w, m_max]

# NTRU Round3
def NTRU_param(n):
    stddev = 0.82
    if n == 508:
        q = 2048; w = 254
    elif n == 676:
        q = 2048; w = 254; stddev = 0.61
    elif n == 820:
        q = 4096; w = 510; stddev = 0.79
    elif n == 652:
        q = 4621; w = 288
    elif n == 760:
        q = 4591; w = 286
    elif n == 856:
        q = 5167; w = 322
    else:
        print("Not in preset params")
    m_max = n
    print("NTRUHPS & NTRUPrime Parameter setting")
    print("n = %d, q = %d, stddev = %.2f, w = %d, m = %d" % (n, q, stddev, w, m_max))
    return [n, q, stddev, w, m_max]

def BLISS_param(n):
    stddev = 0.82
    if n == 512:
        q = 12289; w = 154
    m_max = n
    print("BLISS I + II Parameter setting")
    print("n = %d, q = %d, stddev = %.2f, w = %d, m = %d" % (n, q, stddev, w, m_max))
    return [n, q, stddev, w, m_max]

def GLP_param(n):
    stddev = 0.82
    if n == 512:
        q = 8383489; w = 342
    m_max = n
    print("GLP setting")
    print("n = %d, q = %d, stddev = %.2f, w = %d, m = %d" % (n, q, stddev, w, m_max))
    return [n, q, stddev, w, m_max]

def primal_may(param, lv=0, filename=False):
    return primal_may_inner(param[0], param[1], param[2], param[3], param[4], lv=lv, filename=filename)

def primal_may_inner(n, q, stddev, w, m_max = None,lv=0, filename=False):

    if filename is not False:
        f = open(filename, 'a')
        sys.stdout = f

    print("Input Params")
    print("n = %d, q = %d, stddev = %.2f, w = %d, m = %d" % (n, q, stddev, w, m_max))
    print("----------------------")

    m = m_max
    if m_max is None:
        m = n

    nu = stddev / float(sqrt((w) / n))
    func = partial(cost_zeta, n=n, q=q, stddev=stddev, w=w, m=m, nu=nu, lv=lv)
    best = my_binary_search(0, n, func, cur_depth = 0, max_depth = 3)
    
    print("----------------------")

    prettyprint(best)

    if filename is not False:
        f.close()
    else:
        return best

def cost_zeta(zeta, n, q, stddev, w, m, nu, lv=0, filename=False):

    if filename is not False:
        f = open(filename, 'a')
        sys.stdout = f

    m_ = m + n + 1 - zeta
    d1 = n + 1 - zeta

    numer = comb(n, zeta)
    probs_hw = [comb(n-w, zeta-i) * comb(w, i) / numer for i in range(31)]

    # Find the best submatrix size f < m_
    func = partial(cost_d, zeta=zeta, n=n, q=q, stddev=stddev, w=w, nu=nu, lv=lv, probs_hw = probs_hw)
    best = my_binary_search(d1 + 1, m_, func, cur_depth = 0, max_depth = 3)

    line = ' WIP: zeta = %d : %f' % (zeta, best['cost'])
    print(line)

    if filename is not False:
        f.close()

    return best

def cost_d(d, zeta, n, q, stddev, w, nu, lv, probs_hw = False):

    if d < 40:
        return {'cost': np.inf}

    beta_max = min(d, 500)

    # Find the best blocksize beta
    func = partial(cost_beta, zeta=zeta, n=n, q=q, d=d, stddev=stddev, w=w, nu=nu, lv=lv, probs_hw = probs_hw)
    best = my_binary_search(40, beta_max, func, cur_depth = 0, max_depth = 3)

    print('     zeta = %d / d = %d :' % (zeta, d), best['cost'])

    return best

def cost_beta(beta, zeta, n, q, d, stddev, w, nu, lv, probs_hw = False):

    d1 = n + 1 - zeta
    log_lat, GSnorm = convert_LWE(d, q, d, d1, beta, nu)

    lat = {'cost': round(log_lat, 2)}
    lat['zeta'] = zeta
    lat['d'] = d
    lat['beta'] = beta
    lat['nu'] = round(nu, 3)
    lat['last_R'] = round(GSnorm[-1], 2)

    pr_np = prob_np(GSnorm, stddev)
    if pr_np == 0:
        lat['log_pr_np'] = np.inf
        lat['cost'] = np.inf
        return lat
    log_pr_np = Log2(pr_np)
    lat['log_pr_np'] = log_pr_np

    w_g_min = max((w - (n-zeta)), 0)
    w_g_max = min(zeta, w) # arbitrary bound

    lat['w_g_range'] = [w_g_min, w_g_max]

    # Find the best guessing weight w_g
    func = partial(cost_guess, lat=lat, GSnorm=GSnorm, zeta=zeta, n=n, d=d, stddev=stddev, w=w, lv=lv, probs_hw = probs_hw)
    best = my_binary_search(w_g_min, w_g_max, func, cur_depth = 0, max_depth = 5)

    return best

def cost_guess(w_g, lat, GSnorm, zeta, n, d, stddev, w, lv, probs_hw = False):
    log_lat = lat['cost']
    log_pr_np = lat['log_pr_np']

    current = {'lat': lat.copy()}
    current['guess'] = {}
    current['guess']['w_g'] = w_g

    if w_g == 0:
        pr_hw = probs_hw[0]
        log_pr_hw = Log2(pr_hw)
        current['cost'] = round(log_lat - log_pr_hw - log_pr_np, 2)
    else:
        guess = None
        babai_cost = 2*Log2(d)
        if w_g < len(probs_hw):
            log_pr_hw = Log2(probs_hw[w_g])
        else:
            log_pr_hw = Log2(prob_hw_exact(n, zeta, w, w_g))
        
        current['guess']['log_pr_w'] = log_pr_hw
        if lv == 0: # Exhaustive Search
            guess = {'cost': round(Log2(num_ternary_secret(zeta, w_g)) + babai_cost, 2)}
        else: # MitM or Meet-lwe
            guess = solve_LWE_over_M(zeta, GSnorm, w_g, stddev, 2*stddev, lv=lv)

        current['cost'] = round(max(log_lat, guess['cost']) - log_pr_hw - log_pr_np, 2)
        for key in guess:
            current['guess'][key] = guess[key]

    return current

def convert_LWE(m_, q, d, d1, beta, nu):
    log_lat = log_BKZ_cost(d, beta)
    return log_lat, ext_GSA(m_, q, d, d1, beta, nu=nu)

def solve_LWE_over_M(zeta, GSnorm, w0, b, l1, lv = 1):
    '''
    Given a matrix B (of size d X zeta) and M (of size d X d),
    compute "Noisy-may-mitm cost" to find s (of weight w0) such that [Bs]_M = e
    '''
    best = {'cost': np.inf}
    f = len(GSnorm)
    babai_cost = f**2 # Cost for Babai's nearest plane algorithm

    for w1 in range((w0 + 1) // 2, w0 + 1):
        R1 = ambiguity(zeta, w0, w1)
        if R1 == 0:
            continue

        if lv == 1:
            try:
                p_adm1 = prob_admissible(b, f, GSnorm, unif=False)
                if p_adm1 == 0:
                    continue
                if R1 * p_adm1 < 1:
                    continue
                logL1 = Log2(num_ternary_secret(zeta, w1)) - 0.5 * (Log2(R1) + Log2(p_adm1))
                L1 = 2**logL1

                costlsh1, values1 = lsh_cost(L1, zeta, f, b, l1, w0, w1, GSnorm, unif=False)
                costtop = L1 * babai_cost

                current = {'cost': Log2(costlsh1 + costtop)}
                current['logcost_i'] = [Log2(costlsh1), Log2(costtop)]
                current['logL_i'] = [Log2(L1)]
                current['w_i'] = [w1]
                current['logR_i'] = [Log2(R1)]
                current['logp_adm'] = [Log2(p_adm1)]
                current['lsh_info'] = [values1]
                
                if best is None:
                    best = copy.deepcopy(current)
                elif current['cost'] < best['cost']:
                    best = copy.deepcopy(current)

            except OverflowError:
                continue

        else:
            try:
                k1, p_box1, vol_ratio1 = constraint_area(R1, b, l1, GSnorm, unif=False)
                logL1 = Log2(num_ternary_secret(zeta, w1)) - Log2(vol_ratio1)
                L1 = 2**logL1
                costlsh1, values1 = lsh_cost(L1, zeta, f, b, l1, w0, w1, GSnorm, ks=k1, unif=False)
            except OverflowError:
                continue

            for w2 in range((w1 + 1) // 2, w1 + 1):                
                l2 = 2 * l1
                R2 = ambiguity(zeta, w1, w2)
                if R2 == 0:
                    continue
            
                if lv == 2:
                    try:
                        p_adm2 = prob_admissible(l1, k1, GSnorm)
                        if p_adm2 == 0:
                            continue
                        if R2 * p_adm2 < 1:
                            continue
                        L2 = num_ternary_secret(zeta, w2) / sqrt(R2 * p_adm2)
                        costlsh2, values2 = lsh_cost(L2, zeta, k1, l1, l2, w1, w2, GSnorm)
                        costtop = L2 * babai_cost

                        current = {'cost': Log2(costlsh1 + costlsh2 + costtop)}
                        current['logcost_i'] = [Log2(costlsh1), Log2(costlsh2), Log2(costtop)]
                        current['logL_i'] = [Log2(L1), Log2(L2)]
                        current['w_i'] = [w1, w2]
                        current['logR_i'] = [Log2(R1), Log2(R2)]
                        current['logp_adm'] = [Log2(p_adm2)]
                        current['proj'] = [k1]
                        current['lsh_info'] = [values1, values2]

                        if best is None:
                            best = copy.deepcopy(current)
                        elif current['cost'] < best['cost']:
                            best = copy.deepcopy(current)
                        
                    except OverflowError:
                        continue       

                else:
                    try:
                        k2, p_box2, vol_ratio2 = constraint_area(R2, l1, l2, GSnorm)
                        L2 = num_ternary_secret(zeta, w2) / vol_ratio2
                        costlsh2, values2 = lsh_cost(L2, zeta, k1, l1, l2, w1, w2, GSnorm, ks=k2)

                    except OverflowError:
                        continue    

                    for w3 in range((w2 + 1) // 2, w2 + 1):
                        try:
                            l3 = 2 * l2
                            R3 = ambiguity(zeta, w2, w3)
                            p_adm3 = prob_admissible(l2, k2, GSnorm)
                            if R3 == 0 or p_adm3 == 0:
                                continue
                            if R3 * p_adm3 < 1:
                                continue
                            L3 = num_ternary_secret(zeta, w3) / sqrt(R3 * p_adm3)
                            costlsh3, values3 = lsh_cost(L3, zeta, k2, l2, l3, w2, w3, GSnorm)
                            costtop = L3 * babai_cost

                            current = {'cost': Log2(costlsh1 + costlsh2 + costlsh3 + costtop)}
                            current['logcost_i'] = [Log2(costlsh1), Log2(costlsh2), Log2(costlsh3), Log2(costtop)]
                            current['logL_i'] = [Log2(L1), Log2(L2), Log2(L3)]
                            current['w_i'] = [w1, w2, w3]
                            current['logR_i'] = [Log2(R1), Log2(R2), Log2(R3)]
                            current['logp_adm'] = [Log2(p_adm3)]
                            current['proj'] = [k1, k2]
                            current['lsh_info'] = [values1, values2, values3]

                            if best is None:
                                best = copy.deepcopy(current)
                            elif current['cost'] < best['cost']:
                                best = copy.deepcopy(current)
                                
                        except OverflowError:
                            continue    
                                    
    return best


def constraint_area(R, b, l, GSnorm, unif=True):
    '''
    For a ternary s of weight w0 satisfying ||Bs||_inf < b,
    we have R many ternary pairs (s1, s2) of weight w1 s.t. s = s1 - s2.
    We compute a box A of length l, 
    where there is at least good pair (s1, s2) s.t. both Bs1 & Bs2 lie on A.
    '''
    
    k = 1
    p_adm = 1
    vol_ratio = 1

    while True:
        cur_ratio = 1
        cur_p_adm = 1
        if GSnorm[-k] < l:
            cur_p_adm = prob_admissible_axis(b, GSnorm[-k], unif = unif)
        else:
            cur_ratio = GSnorm[-k] / l
            cur_p_adm = prob_admissible_axis(b, l, unif=unif) / cur_ratio

        if k == len(GSnorm):
            break
        
        if vol_ratio * cur_ratio > R * cur_p_adm:
            k -= 1
            break
        else:
            vol_ratio *= cur_ratio
            p_adm *= cur_p_adm
            k += 1

    return k, p_adm, vol_ratio

def prob_lsh(L, d, kl, b, l, w1, w2, GSnorm, ks = 1, unif=True):
    '''
    Compute the dimension of torus-LSH that makes false-positive probability as small as possible,
    so that near-collision finding takes not so much time.

    b : distance of near-collision
    l : length of box for (torus-)LSH
    * For the axis where GSnorm is small, define b_ = min(b, GSnorm) and l = min(l, GSnorm).

    Except top level, the list has a constraint on the length of the last "ks" axis.
    As we take LSH length = constraint length, any pair trivially collides on the last "ks" axis.
    '''

    dim_torus = ks
    p_torus_good = 1
    p_torus_bad = 1

    while True:
        cur_len = l
        if GSnorm[-dim_torus] < l:
            cur_len = GSnorm[-dim_torus]
        cur_p_torus_good = prob_admissible_axis(b, cur_len, unif=unif)
        cur_p_torus_bad = cur_len / GSnorm[-dim_torus]
        
        if L * p_torus_bad * cur_p_torus_bad < 1 :
            dim_torus -= 1
            break
        else:
            p_torus_good *= cur_p_torus_good
            p_torus_bad *= cur_p_torus_bad
            if dim_torus == kl:
                break
            dim_torus += 1

    exp_bad = L * p_torus_bad

    values = {}
    values['dim_torus'] = dim_torus
    try:
        values['log_p_good'] = Log2(p_torus_good)
    except ValueError:
        values['log_p_good'] = -np.inf
    try:
        values['log_p_bad'] = Log2(p_torus_bad)
    except ValueError:
        values['log_p_bad'] = -np.inf
    
    return values

def num_ternary_secret_half(n, w):
    return multinom(n, [w, w, n - 2*w])

def num_ternary_secret(n, w):
    return 2**w * comb(n, w)

def lsh_cost(L, d, kl, b, l, w1, w2, GSnorm, ks = 1, unif=True):
    values = prob_lsh(L, d, kl, b, l, w1, w2, GSnorm, ks = ks, unif=unif).copy()
    log_cost = 2*Log2(L) + values['log_p_bad'] - values['log_p_good']
    try:
        return 2**log_cost, values
    except OverflowError:
        return np.inf, values

def prob_hd_good(d, k, w):
    return easy_multinom(d - k, w, w)

def prob_hd_bad(d, k, w):
    u = min(w, d)
    prob = 0
    for i in range(u):
        for j in range(i):
            prob += easy_multinom(k, i, j) * easy_multinom(d-k, w-i, w-j)**2
    prob /= easy_multinom(d, w, w)**2
    return prob