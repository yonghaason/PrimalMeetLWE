from utils import *
from functools import partial
import copy

"""

* How to use?
In 'python' shell (NOT 'sage' shell!):
> from estimate import *
> # param = [(n), (q), (error_type), (error_param), (is_secret_balanced), (secret_param), (m)] 
> # primal_may(param, (lv))
> param = [2**15, 2**768, 'gaussian', 3.19, True, 96, 2**15]
> primal_may(param, lv=2)

- It takes LONG LONG time (over several hours!) for high level, or large LWE parameters. 
- If the reviewer is curious about the detailed parameters,
we recommend to see the text files in `logs' directory.
- Checked with python3.8

* About the consistency with `lattice-estimator:
Our level-0 attack strategy is covered by `primal_hybrid' with options `mitm=False, babai=True' in lattice-estimator.
We observe that both scripts output similar numbers for this strategy.
Meanwhile, our level-1 attack strategy is almost the same as the primal hybrid-MitM attack [HG07].
This attack is covered by 'primal_hybrid' with options 'mitm=True, babai=True' in 'lattice-estimator'.
However, we found that the computation is over-simplified from the original description in [HG07],
and we believe our script more accurately computes the attack complexity.
This results in quite different estimation of this attack for PQC parameters.
However, it has no actual impact on the current parameter selections
since resulting numbers are non-competitive compared to other lattice attacks.
For FHE parameters, two estimations are quite similar
and there is no notable problem from lattice-estimator's over-simplified estimation.

"""

def primal_may(param, lv, filename=False):
    
    if filename is not False:
        f = open(filename, 'a')
        sys.stdout = f

    n = param[0]
    q = param[1]
    error_type = param[2]
    error_param = param[3]
    is_secret_balanced = param[4]
    w = param[5]
    m = param[6]

    if error_type == 'unif':
        stddev = sqrt(((2*error_param+1)**2 - 1)/12)
    elif error_type == 'gaussian':
        stddev = error_param
    else:
        print("Please specify error type: `unif' or `gaussian'.")
        return False

    print("----------------------")
    print("Input Params")
    if q > 10000:
        print("n = %d, logq = %d, stddev = %.2f, m = %d" % (n, log(q,2), stddev, m))
    else:
        print("n = %d, q = %d, stddev = %.2f, m = %d" % (n, q, stddev, m))
    if is_secret_balanced:
        print("HW(s) = %d, where %d 1s and %d -1s" % (2*w, w, w))
    else:
        print("HW(s) = %d" % w)
    if error_type == 'unif':
        print("(The input uniform error distribution is considered as a gaussian with the same variance)")
    print("----------------------")

    return primal_may_inner(n, q, stddev, is_secret_balanced, w, m, lv=lv, filename=filename)

def primal_may_inner(n, q, stddev, is_sec_bal, w, m, lv, filename=False):

    if filename is not False:
        f = open(filename, 'a')
        sys.stdout = f

    nu = stddev / float(sqrt((w) / n))
    func = partial(cost_zeta, n=n, q=q, stddev=stddev, 
                    is_sec_bal=is_sec_bal, w=w, m=m, nu=nu, lv=lv)
    best = my_binary_search(0, n, func, cur_depth = 0, max_depth = 4)
    
    print("----------------------")

    prettyprint(best)

    if filename is not False:
        f.close()
    else:
        return best

def cost_zeta(zeta, n, q, stddev, is_sec_bal, w, m, nu, lv, filename = False):

    if filename is not False:
        f = open(filename, 'a')
        sys.stdout = f

    m_ = m + n + 1 - zeta
    d1 = n + 1 - zeta

    probs_hw = probs_hw_precompute(n, zeta, w, is_sec_bal)        

    # Find the best submatrix size d < m_
    func = partial(cost_d, zeta=zeta, n=n, q=q, stddev=stddev, 
                    is_sec_bal=is_sec_bal, w=w, nu=nu, lv=lv, probs_hw = probs_hw)
    best = my_binary_search(d1 + 1, m_, func, cur_depth = 0, max_depth = 2)
    
    line = ' WIP: zeta = %d : %f' % (zeta, best['cost'])
    print(line)

    if filename is not False:
        f.close()

    return best

def cost_d(d, zeta, n, q, stddev, is_sec_bal, w, nu, lv, probs_hw = []):

    if d < 40:
        return {'cost': np.inf}

    beta_max = min(d, 500)

    # Find the best blocksize beta
    func = partial(cost_beta, zeta=zeta, n=n, q=q, d=d, stddev=stddev, 
                    is_sec_bal=is_sec_bal, w=w, nu=nu, lv=lv, probs_hw = probs_hw)
    best = my_binary_search(40, beta_max, func, cur_depth = 0, max_depth = 3)

    print('     zeta = %d / d = %d :' % (zeta, d), best['cost'])
    return best

def cost_beta(beta, zeta, n, q, d, stddev, is_sec_bal, w, nu, lv, probs_hw = []):

    log_lat = log_BKZ_cost(d, beta)
    GSnorm = GSA(d, q, d,  n + 1 - zeta, beta, nu)

    lat = {'cost': round(log_lat, 2)}
    lat['zeta'] = zeta
    lat['d'] = d
    lat['beta'] = beta
    lat['nu'] = round(nu, 3)
    lat['last_R'] = round(GSnorm[-1], 2)

    pr_np = prob_np(GSnorm, stddev)

    # Skip if the last GSnorm is too short
    if pr_np == 0 or stddev > GSnorm[-1]/2:
        lat['log_pr_np'] = np.inf
        lat['cost'] = np.inf
        return lat
    log_pr_np = Log2(pr_np)
    lat['log_pr_np'] = log_pr_np

    # rough bound assuming |S|^0.25 < lat['cost']
    w_g_bound = ceil(4*lat['cost']/(1 + Log2(zeta)))
    if is_sec_bal:
        w_g_bound = (w_g_bound+1)//2

    w_g_min = max((w - (n-zeta)), 0)
    w_g_max = min(w_g_bound, zeta, w)
    
    lat['w_g_range'] = [(w_g_min//2)*2, ((w_g_max+1)//2)*2]

    # Find the best guessing weight w_g
    func = partial(cost_guess_rec, lat=lat, GSnorm=GSnorm, zeta=zeta, n=n, d=d, stddev=stddev, 
                    is_sec_bal=is_sec_bal, w=w, lv=lv, probs_hw = probs_hw)
    best = my_binary_search(w_g_min//2, w_g_max//2, func, cur_depth = 0, max_depth = 4)

    return best

def cost_guess(w_g_half, lat, GSnorm, zeta, n, d, stddev, is_sec_bal, w, lv, probs_hw = []):
    # make w_g even
    w_g = 2*w_g_half

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
        if w_g >= len(probs_hw):
            log_pr_hw = Log2(prob_hw(n, zeta, w, w_g, is_sec_bal))
        else:
            log_pr_hw = Log2(probs_hw[w_g])
        
        current['guess']['log_pr_w_g'] = log_pr_hw
        if lv == 0: # Exhaustive Search
            guess = {'cost': round(Log2(num_ternary_secret(zeta, w_g, is_sec_bal)) \
                + 2*Log2(babai_cost(d)), 2)}
        else: # MitM or Meet-lwe
            guess = solve_LWE_over_M(zeta, GSnorm, is_sec_bal, w_g, stddev, lv=lv, bound=log_lat)

        current['cost'] = round(max(log_lat, guess['cost']) - log_pr_hw - log_pr_np, 2)
        for key in guess:
            current['guess'][key] = guess[key]

    return current

def solve_LWE_over_M(zeta, GSnorm, is_sec_bal, w0, stddev, lv, bound=np.inf):
    '''
    Given a matrix B (of size d X zeta) and M (of size d X d),
    compute "Noisy-may-mitm cost" to find s (of weight w0) such that [Bs]_M = e
    where e is sampled from Gaussian(stddev)

    TODO 1. Change this to perform binary search, rather than current full range search.
    TODO 2. Clean up by refactoring this function with recursive call ..
    '''

    best = {'cost': np.inf}
    d = len(GSnorm)
    
    constraint_bound0 = stddev
    lsh_length1 = 6*stddev

    if lv == 1:
        p_adm = 1
        for i in range(d):
            p_adm *= prob_admissible_gaussian(GSnorm[-i-1], stddev)
        if p_adm == 0:
            return best
        for w1 in range(w0//2, w0 + 1):
            R1 = ambiguity(zeta, w0, w1, is_sec_bal)
            if R1 * p_adm < 1:
                continue
            S1 = num_ternary_secret(zeta, w1, is_sec_bal)
            logL1 = Log2(S1) - 0.5 * (Log2(R1) + Log2(p_adm))
            L1 = 2**logL1
            costlsh1, values1 = \
                lsh_cost(L1, constraint_bound0, lsh_length1, GSnorm, 
                        lsh_start_idx = 0, 
                        constraint_cube_length = 0, 
                        constraint_dim_idx = d, unif=False)
            costtop = L1 * babai_cost(d)

            current = {'cost': Log2(costlsh1 + costtop)}
            current['logcost_i'] = [Log2(costlsh1), Log2(costtop)]
            current['logS_i'] = [Log2(S1)]
            current['logL_i'] = [Log2(L1)]
            current['w_i'] = [w1]
            current['logR_i'] = [Log2(R1)]
            current['logp_adm'] = [Log2(p_adm)]
            current['cons_bound'] = [constraint_bound0]
            current['lsh_info'] = [values1]
            
            if best is None:
                best = copy.deepcopy(current)
            elif current['cost'] < best['cost']:
                best = copy.deepcopy(current)

    else: # lv >= 2
        # Area_1: [-constraint_bound1, constraint_bound1]^proj_dim1
        constraint_bound1 = constraint_bound0 
        lsh_length2 = 2*constraint_bound1

        w1_start = (((w0//2)+1)//2)*2
        for w1 in range(w1_start, w0 + 1, 2):            
            R1 = ambiguity(zeta, w0, w1, is_sec_bal)
            proj_dim1 = 0
            proj_dim1, p_rep1, vol_ratio1 = \
                constraint_area(R1, constraint_bound0, constraint_bound1, GSnorm, unif=False)
            S1 = num_ternary_secret(zeta, w1, is_sec_bal)
            logL1 = Log2(S1) - Log2(vol_ratio1)
            L1 = 2**logL1
            costlsh1, values1 = \
                lsh_cost(L1, constraint_bound0, lsh_length1, GSnorm, 
                    lsh_start_idx = 0, 
                    constraint_cube_length = 2*constraint_bound1, 
                    constraint_dim_idx = d-proj_dim1, unif=False)
            if proj_dim1 == 0 or Log2(costlsh1) > bound:
                continue

            if lv == 2:
                p_adm = 1
                for i in range(proj_dim1):
                    p_adm *= prob_admissible_uniform(GSnorm[-i-1], constraint_bound1)
                if p_adm == 0:
                    continue
                for w2 in range(w1//2, w1 + 1):
                    R2 = ambiguity(zeta, w1, w2, is_sec_bal)
                    if R2 * p_adm < 1:
                        continue
                    S2 = num_ternary_secret(zeta, w2, is_sec_bal)
                    L2 = S2 / sqrt(R2 * p_adm)
                    costlsh2, values2 = \
                        lsh_cost(L2, constraint_bound1, lsh_length2, GSnorm, 
                                lsh_start_idx = d-proj_dim1, 
                                constraint_cube_length = 0, 
                                constraint_dim_idx = d, unif=True)
                    costtop = L2 * babai_cost(d)

                    current = {'cost': Log2(costlsh1 + costlsh2 + costtop)}
                    current['logcost_i'] = [Log2(costlsh1), Log2(costlsh2), Log2(costtop)]
                    current['logS_i'] = [Log2(S1), Log2(S2)]
                    current['logL_i'] = [Log2(L1), Log2(L2)]
                    current['w_i'] = [w1, w2]
                    current['logR_i'] = [Log2(R1), Log2(R2)]
                    current['logp_rep'] = [Log2(p_rep1)]
                    current['logvol_ratio'] = [Log2(vol_ratio1)]
                    current['logp_adm'] = [Log2(p_adm)]
                    current['proj_dim'] = [proj_dim1]
                    current['cons_bound'] = [constraint_bound0, constraint_bound1]
                    current['lsh_info'] = [values1, values2]

                    if best is None:
                        best = copy.deepcopy(current)
                    elif current['cost'] < best['cost']:
                        best = copy.deepcopy(current)

            else:
                constraint_bound2 = constraint_bound1
                lsh_length3 = 2*constraint_bound2
                w2_start = (((w1//2)+1)//2)*2
                for w2 in range(w2_start, w1 + 1, 2):
                    R2 = ambiguity(zeta, w1, w2, is_sec_bal)            
                    proj_dim2, p_rep2, vol_ratio2 = \
                        constraint_area(R2, constraint_bound1, constraint_bound2, GSnorm, unif=True, proj_last = proj_dim1)
                    S2 = num_ternary_secret(zeta, w2, is_sec_bal)
                    logL2 = Log2(S2) - Log2(vol_ratio2)
                    L2 = 2**logL2
                    costlsh2, values2 = \
                        lsh_cost(L2, constraint_bound1, lsh_length2, GSnorm, 
                                lsh_start_idx = d-proj_dim1,
                                constraint_cube_length = 2*constraint_bound2, 
                                constraint_dim_idx = d-proj_dim2, unif=True)
                    if proj_dim2 == 0 or Log2(costlsh2) > bound:
                        continue

                    
                    p_adm = 1
                    for i in range(proj_dim2):
                        p_adm *= prob_admissible_uniform(GSnorm[-i-1], constraint_bound2)
                    if p_adm == 0:
                        continue
                
                    for w3 in range(w2//2, w2 + 1):
                        R3 = ambiguity(zeta, w2, w3, is_sec_bal)            
                        if R3 * p_adm < 1:
                            continue
                        S3 = num_ternary_secret(zeta, w3, is_sec_bal)
                        L3 = S3 / sqrt(R3 * p_adm)
                        costlsh3, values3 = \
                            lsh_cost(L3, constraint_bound2, lsh_length3, GSnorm, 
                                    lsh_start_idx=d-proj_dim2,
                                    constraint_cube_length = 0, 
                                    constraint_dim_idx = d, unif=True)
                        costtop = L3 * babai_cost(d)

                        current = {'cost': Log2(costlsh1 + costlsh2 + costlsh3 + costtop)}
                        current['logcost_i'] = [Log2(costlsh1), Log2(costlsh2), Log2(costlsh3), Log2(costtop)]
                        current['logS_i'] = [Log2(S1), Log2(S2), Log2(S3)]
                        current['logL_i'] = [Log2(L1), Log2(L2), Log2(L3)]
                        current['w_i'] = [w1, w2, w3]
                        current['logR_i'] = [Log2(R1), Log2(R2), Log2(R3)]
                        current['logp_rep'] = [Log2(p_rep1), Log2(p_rep2)]
                        current['logvol_ratio'] = [Log2(vol_ratio1), Log2(vol_ratio2)]
                        current['logp_adm'] = [Log2(p_adm)]
                        current['proj_dim'] = [proj_dim1, proj_dim2]
                        current['cons_bound'] = [constraint_bound0, constraint_bound1, constraint_bound2]
                        current['lsh_info'] = [values1, values2, values3]

                        if best is None:
                            best = copy.deepcopy(current)
                        elif current['cost'] < best['cost']:
                            best = copy.deepcopy(current)
                                    
    return best

def cost_beta_rec(beta, zeta, n, q, d, stddev, is_sec_bal, w, nu, lv, probs_hw = []):

    log_lat = log_BKZ_cost(d, beta)
    GSnorm = GSA(d, q, d,  n + 1 - zeta, beta, nu)

    lat = {'cost': round(log_lat, 2)}
    lat['zeta'] = zeta
    lat['d'] = d
    lat['beta'] = beta
    lat['nu'] = round(nu, 3)
    lat['last_R'] = round(GSnorm[-1], 2)

    pr_np = prob_np(GSnorm, stddev)

    # Skip if the last GSnorm is too short
    if pr_np == 0 or stddev > GSnorm[-1]/2:
        lat['log_pr_np'] = np.inf
        lat['cost'] = np.inf
        return lat
    log_pr_np = Log2(pr_np)
    lat['log_pr_np'] = log_pr_np

    # rough bound assuming |S|^0.25 < lat['cost']
    w_g_bound = ceil(4*lat['cost']/(1 + Log2(zeta)))
    if is_sec_bal:
        w_g_bound = (w_g_bound+1)//2

    w_g_min = max((w - (n-zeta)), 0)
    w_g_max = min(w_g_bound, zeta, w)
    
    lat['w_g_range'] = [(w_g_min//2)*2, ((w_g_max+1)//2)*2]

    # Find the best guessing weight w_g
    func = partial(cost_guess_rec, lat=lat, GSnorm=GSnorm, zeta=zeta, n=n, d=d, stddev=stddev, 
                    is_sec_bal=is_sec_bal, w=w, lv=lv, probs_hw = probs_hw)
    best = my_binary_search(w_g_min//2, w_g_max//2, func, cur_depth = 0, max_depth = 4)

    return best

def cost_guess_rec(w_g_half, lat, GSnorm, zeta, n, d, stddev, is_sec_bal, w, lv, probs_hw = []):
    # make w_g even
    w_g = 2*w_g_half

    log_lat = lat['cost']
    log_pr_np = lat['log_pr_np']

    current = {'lat': lat.copy()}
    current['guess'] = {}
    current['guess']['w_g'] = w_g

    log_pr_hw = None
    if w_g >= len(probs_hw):
        log_pr_hw = Log2(prob_hw(n, zeta, w, w_g, is_sec_bal))
    else:
        log_pr_hw = Log2(probs_hw[w_g])

    if w_g == 0:
        current['guess']['log_pr_w_g'] = log_pr_hw
        current['cost'] = round(log_lat - log_pr_hw - log_pr_np, 2)
    else:
        guess = None
        babai_cost = 2*Log2(d)
        current['guess']['log_pr_w_g'] = log_pr_hw
        if lv == 0: # Exhaustive Search
            guess = {'cost': round(Log2(num_ternary_secret(zeta, w_g, is_sec_bal)) + babai_cost, 2)}
        else: # MitM or Meet-lwe
            current_guess = \
             {'cost_i': [], 
             'S_i': [], 
             'L_i': [], 
             'R_i': [], 
             'w_i': [],
             'p_rep': [], 
             'vol_ratio': [], 
             'proj_dim': [], 
             'cons_bound': [], 
             'lsh_stats': []}
            guess = meet_LWE_cost(1, GSnorm, is_sec_bal, zeta, 
                proj_dim = len(GSnorm), w = w_g, is_error_unif = False, 
                error_param = stddev, lsh_length = 6*stddev, 
                last_lv=lv, current_guess = copy.deepcopy(current_guess), abort_bound=log_lat)

        current['cost'] = round(max(log_lat, guess['cost']) - log_pr_hw - log_pr_np, 2)
        for key in guess:
            current['guess'][key] = guess[key]

    return current
    
def meet_LWE_cost(lv, GSnorm, is_sec_bal, zeta,
    proj_dim, w, is_error_unif, error_param, 
    lsh_length, last_lv, current_guess, abort_bound = np.inf):
    
    d = len(GSnorm)
    cur_best = {'cost': np.inf}

    if lv == last_lv:
        p_adm = 1
        if is_error_unif:
            for i in range(proj_dim):
                p_adm *= prob_admissible_uniform(GSnorm[-1-i], error_param)
        else:
            for i in range(proj_dim):
                p_adm *= prob_admissible_gaussian(GSnorm[-1-i], error_param)
        if p_adm == 0:
            return cur_best

        for w_split in range(w//2, w + 1):
            R = ambiguity(zeta, w, w_split, is_sec_bal)
            if R * p_adm < 1:
                continue
            S = num_ternary_secret(zeta, w_split, is_sec_bal)
            L = S / sqrt(R*p_adm)
            cost_lsh, stat_lsh = lsh_cost(L, error_param, lsh_length, GSnorm,
                                lsh_start_idx = d - proj_dim, unif=is_error_unif)
            cost_top = L * babai_cost(proj_dim)

            total_cost = 0.0
            for logcost_i in current_guess['cost_i']:
                total_cost += 2**logcost_i

            if Log2(total_cost) < cur_best['cost']:
                cur_best = copy.deepcopy(current_guess)
                cur_best['cost'] = Log2(total_cost)
                cur_best['cost_i'].append(Log2(cost_lsh))
                cur_best['cost_i'].append(Log2(cost_top))
                cur_best['S_i'].append(Log2(S))
                cur_best['L_i'].append(Log2(L))
                cur_best['R_i'].append(Log2(R))
                cur_best['w_i'].append(w_split)
                cur_best['p_adm'] = (Log2(p_adm))
                cur_best['lsh_stats'].append(stat_lsh)

        return cur_best
        
    else:
        constraint_bound = error_param # One can further fine-tune this ...
        w_start = w//2
        if w_start % 2 == 1: w_start += 1
        for w_split in range(w_start, w + 1, 2):
            R = ambiguity(zeta, w, w_split, is_sec_bal)
            next_proj_dim, p_rep, vol_ratio = \
                constraint_area(R, error_param, constraint_bound, GSnorm,
                                unif=is_error_unif, 
                                proj_last = proj_dim)
            if next_proj_dim == 0:
                continue
            S = num_ternary_secret(zeta, w_split, is_sec_bal)
            L = S / vol_ratio
            cost_lsh, stat_lsh = lsh_cost(L, error_param, lsh_length, GSnorm,
                                lsh_start_idx = d - proj_dim,
                                constraint_cube_length = 2*constraint_bound,
                                constraint_dim_idx = d - next_proj_dim, unif=is_error_unif)
            if Log2(cost_lsh) > abort_bound:
                continue

            current_guess['cost_i'].append(Log2(cost_lsh))
            current_guess['S_i'].append(Log2(S))
            current_guess['L_i'].append(Log2(L))
            current_guess['R_i'].append(Log2(R))
            current_guess['w_i'].append(w_split)
            current_guess['p_rep'].append(Log2(p_rep))
            current_guess['vol_ratio'].append(Log2(vol_ratio))
            current_guess['proj_dim'].append(next_proj_dim)
            current_guess['cons_bound'].append(error_param)
            current_guess['lsh_stats'].append(stat_lsh)

            comp_best = meet_LWE_cost(lv+1, GSnorm, is_sec_bal, zeta,
            next_proj_dim, w_split, True, constraint_bound, 2*constraint_bound, 
            last_lv, copy.deepcopy(current_guess), abort_bound)

            if comp_best['cost'] < cur_best['cost']:
                cur_best = copy.deepcopy(comp_best)            
        
        return cur_best

def constraint_area(R, error_param, constraint_bound, GSnorm, unif=True, proj_last = None):
    '''
    Consider a secret s satisfying Bs = e where e is determined by 'error_param',
    and assume that there are R many ternary pairs (s1, s2) such that s = s1 - s2.
    This function computes a dimension `proj_dim' r satisfying 
    there is at least one pair (s1, s2) such that 
    C = [-b, b]^r (b = 'constraint_bound') contains both Bs1 & Bs2.
    '''
    
    proj_dim = 1
    p_adm = 1
    vol_ratio = 1
    if proj_last == None:
        proj_last = len(GSnorm)
    
    while True:
        if proj_dim == proj_last:
            break

        cur_ratio = 1
        cur_p_adm = 1
        cur_axis_length = GSnorm[-proj_dim]
        if 2*constraint_bound < GSnorm[-proj_dim]:
            cur_axis_length = 2*constraint_bound
            cur_ratio = GSnorm[-proj_dim] / (2*constraint_bound)
        if unif:
            cur_p_adm = prob_admissible_uniform(cur_axis_length, error_param)
        else:
            cur_p_adm = prob_admissible_gaussian(cur_axis_length, error_param)

        inv = 1.0 / (p_adm * cur_p_adm / cur_ratio)
        if R < inv:
            proj_dim -= 1
            break
        else:
            vol_ratio *= cur_ratio
            p_adm *= cur_p_adm / cur_ratio
            proj_dim += 1

    return proj_dim, p_adm, vol_ratio

def prob_lsh(L, error_param, lsh_length, GSnorm, lsh_start_idx, \
             constraint_cube_length = 0, constraint_dim_idx = np.inf, unif = True):
    '''
    Compute the torus-LSH target dimension which
    minimizes the running time of the torus-LSH based near-collision finding algorithm .

    unif : distance is uniform (true) or gaussian (false)?
    error_param: error bound or stddev of near-collision pair distance
    lsh_length: length of torus-LSH box
    lsh_start_idx: the first coordinate where the lower level list's hypercube is defined

    Except top level list, the list is defined with some hypercube constraints 
    parameterized by `constraint_dim_idx' and `constraint_cube_length',
    which makes each element has small length in last coordinates.

    * For the axis where GSnorm is too small, 
    define b_ = min(b, GSnorm) and l = min(l, GSnorm).
    '''

    lsh_last_idx = lsh_start_idx
    p_torus_good = 1
    p_torus_bad = 1

    while True:
        cur_len = lsh_length
        if lsh_last_idx >= constraint_dim_idx:
            cur_len = constraint_cube_length
        if GSnorm[lsh_last_idx] < cur_len:
            cur_len = GSnorm[lsh_last_idx]
        if unif:
            cur_p_torus_good = prob_admissible_uniform(cur_len, error_param)
        else:
            cur_p_torus_good = prob_admissible_gaussian(cur_len, error_param)
        cur_p_torus_bad = cur_len / GSnorm[lsh_last_idx]
        
        if L * p_torus_bad * cur_p_torus_bad < 1 :
            if lsh_last_idx != lsh_start_idx:
                lsh_last_idx -= 1
            else:
                p_torus_good = cur_p_torus_good
                p_torus_bad = cur_p_torus_bad
            break
        else:
            p_torus_good *= cur_p_torus_good
            p_torus_bad *= cur_p_torus_bad
            if lsh_last_idx == len(GSnorm)-1:
                break
            lsh_last_idx += 1

    values = {}
    values['box_length'] = lsh_length
    values['target_coords'] = [lsh_start_idx, lsh_last_idx]
    try:
        values['log_p_good'] = Log2(p_torus_good)
    except ValueError:
        values['log_p_good'] = -np.inf
    try:
        values['log_p_bad'] = Log2(p_torus_bad)
    except ValueError:
        values['log_p_bad'] = -np.inf
    
    return values

def lsh_cost(L, e, lsh_length, GSnorm, lsh_start_idx, \
             constraint_cube_length = 0, constraint_dim_idx = np.inf, unif=True):
    values = prob_lsh(L, e, lsh_length, GSnorm, lsh_start_idx, \
                      constraint_cube_length, constraint_dim_idx, unif).copy()
    babai_dim = len(GSnorm) - lsh_start_idx
    log_cost = 2*Log2(L) + values['log_p_bad'] - values['log_p_good'] + Log2(babai_cost(babai_dim))
    log_cost = max(Log2(L), log_cost)
    try:
        return 2**log_cost, values
    except OverflowError:
        return np.inf, values