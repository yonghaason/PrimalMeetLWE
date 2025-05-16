from utils import *
from functools import partial
import copy

def primal_may(param, t, filename=False, param_specified = {}):
    
    if filename is not False:
        f = open(filename, 'a')
        sys.stdout = f

    n = param[0]
    q = param[1]
    error_type = param[2]
    error_param = param[3]
    w = param[4]
    m_ = param[5]

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
        print("n = %d, logq = %d, stddev = %.2f, m = %d, HW(s) = %d" % (n, log(q,2), stddev, m_, w))
    else:
        print("n = %d, q = %d, stddev = %.2f, m = %d, HW(s) = %d" % (n, q, stddev, m_, w))
    print("Top level = %d" % t)
    if error_type == 'unif':
        print("(The input uniform error distribution is considered as a gaussian with the same variance)")
    print("----------------------")

    return primal_may_inner(n, q, stddev, w, m_, t=t, param_specified = param_specified)

def primal_may_inner(n, q, stddev, w, m_, t, param_specified):

    func = partial(cost_d, n=n, q=q, stddev=stddev, 
                    w=w, m_=m_, t=t, param_specified = param_specified)
    best = None
    if 'd' in param_specified:
        d = param_specified['d']
        best = my_binary_search(d, d, func, cur_depth = 0, max_depth = 6)
    else:
        best = my_binary_search(0, n, func, cur_depth = 0, max_depth = 6)
    
    print("----------------------")

    prettyprint(best)

    return best

def cost_d(d, n, q, stddev, w, m_, t, param_specified):

    m_max = m_ + n + 1 - d
    m_min = n + 1 - d

    probs_hw = probs_hw_precompute(n, d, w)  

    line = ' Loop with d = %d' % (d) 
    print(line)

    # Find the best submatrix size d < m_
    func = partial(cost_m, d=d, n=n, q=q, stddev=stddev, 
                w=w, t=t, probs_hw = probs_hw, param_specified = param_specified)
    if 'm' in param_specified:
        m = param_specified['m']
        best = my_binary_search(m, m, func, cur_depth = 0, max_depth = 6)
    else:
        best = my_binary_search(m_min + 1, m_max, func, cur_depth = 0, max_depth = 5)
    
    line = ' └──── Best: %f' % best['cost']
    print(line)

    return best

def cost_m(m, d, n, q, stddev, w, t, param_specified, probs_hw = []):

    if m < 40:
        return {'cost': np.inf}

    # print(m)

    beta_max = min(m, 800)

    # Find the best blocksize beta
    func = partial(cost_beta, d=d, n=n, q=q, m=m, stddev=stddev, 
                    w=w, t=t, probs_hw = probs_hw, param_specified = param_specified)
    if 'beta' in param_specified:
        beta = param_specified['beta']
        best = my_binary_search(beta, beta, func, cur_depth = 0, max_depth = 4)
    else:
        best = my_binary_search(40, beta_max, func, cur_depth = 0, max_depth = 4)
    
    print(' │ m = %d:' % m, best['cost'])
    return best

def cost_beta(beta, d, n, q, m, stddev, w, t, param_specified, probs_hw = []):

    log_lat = log_BKZ_cost(m, beta)
    GSnorm = GSA(q, m, n + 1 - d, beta)

    lat = {'cost': round(log_lat, 2)}
    lat['d'] = d
    lat['m'] = m
    lat['beta'] = beta
    lat['last_R'] = round(GSnorm[-1], 2)

    pr_np = prob_np(GSnorm, stddev)

    # Skip if the last GSnorm is too short
    if pr_np == 0 or 6*stddev > GSnorm[-1]/2:
        lat['log_pr_np'] = np.inf
        lat['cost'] = np.inf
        return lat
    log_pr_np = Log2(pr_np)
    lat['log_pr_np'] = log_pr_np

    # rough bound assuming |S|^0.25 < lat['cost']
    w_g_bound = 0
    if t > 1:
        while True:
            cand = comb(d, w_g_bound) * 2**w_g_bound
            if Log2(cand)/4 < log_lat:
                w_g_bound += 1
            else:
                break
            if w_g_bound == w:
                break
    else:
        while True:
            cand = comb(d, w_g_bound) * 2**w_g_bound
            if Log2(cand)/2 < log_lat:
                w_g_bound += 1
            else:
                break
            if w_g_bound == w:
                break
            
    w_g_min = max(w - (n-d), 0)
    w_g_max = min(w_g_bound, d, w)
    
    if w_g_min > w_g_max:
        return {'cost': np.inf}

    cost = {'lat': lat.copy()}

    guess = {'cost': np.inf}
    pr_hw = 0.0
    w_g_range = range(((w_g_min+1)//2)*2, (w_g_max//2)*2, 2)

    if 'w_g' in param_specified:
        w_g_range = range(param_specified['w_g'], param_specified['w_g']+1, 2)
    
    ell0 = round(6*stddev, 2)
    b_lsh0 = round(12*stddev, 2)

    # print(w_g_range)
    
    # Our (Generalized) Meet-LWE
    for w_g in w_g_range:
        if w_g == 0:
            continue

        initial_stat = \
        {'cost': np.inf,
        'log_p_suc': None,
        'w': [w_g], 
        'r': [], 
        'ell': [ell0],
        'b_lsh': ['-'],
        'R_lsh': ['-'],
        'log_cost_list': ['-'], 
        'log_cost_lsh': ['-'],  
        'log_fullS': [Log2(num_ternary_secret(d, w_g))],
        'log_S': ['-'], 
        'log_R': ['-'], 
        'log_p_sp': ['-'], 
        'log_p_ncf': ['-'],
        'lsh_dim': ['-'],
        'log_p_col': ['-']}
        cur_guess = meet_LWE_cost(0, GSnorm, d, 
            r = len(GSnorm), w = w_g, e = stddev, ell = ell0, b_lsh = b_lsh0,
            C_proj = 10.0, C_lsh = 10.0, t=t, 
            current_guess = copy.deepcopy(initial_stat), abort_bound=log_lat)

        if cur_guess['cost'] > log_lat:
            if guess['cost'] is not np.inf:
                break

        else:
            if w_g >= len(probs_hw):
                pr_hw = prob_hw(n, d, w, w_g)
            else:
                pr_hw = probs_hw[w_g]
            guess = cur_guess
            cost['w_guess'] = w_g
            
    cost['guess_cost'] = round(guess['cost'], 2)
    cost['log_p_hw'] = Log2(pr_hw)    
    cost['cost'] = round(max(log_lat, guess['cost']) - Log2(pr_hw) - log_pr_np, 2)
    cost['guess_details'] = guess

    return cost
    
def meet_LWE_cost(cur_lv, GSnorm, d,
    r, w, e, ell, b_lsh, C_proj, C_lsh, t,
    current_guess, abort_bound = np.inf, verbose=False):
    
    cur_best = copy.deepcopy(current_guess)

    w_start = w//2
    w_range = None
    
    if cur_lv == t-1:
        # For a quick estimation, we fix w[t] = w[t-1]/2, 
        # rather than try every w[t] in [w[t-1]/2, w[t-1]].
        w_range = range(w_start, w_start + 1)
        # w_range = range(w_start, w)
    else:
        if w_start % 2 == 1: 
            w_start += 1
        w_range = range(w_start, w, 2)

    if verbose:
        if cur_lv == 0: 
            print()
        print(" " * cur_lv, "LV = ",  cur_lv, " (t=", t, "), w", cur_lv, "=", w, ", w_range=", w_range, sep ="")
        print(" " * cur_lv, "error_param=", e, ", norm_bound=", ell, ", lsh_length=", b_lsh, sep ="")
    
    for w_next in w_range:
        running_stat = copy.deepcopy(current_guess)

        R = ambiguity(d, w, w_next)
        S = num_ternary_secret(d, w_next)
        cur_r = r
        next_r = 0

        if cur_lv == t-1:
            # Top-level Optimization
            S *= sqrt(3) / sqrt(R)
        else:
            next_r, p_rep, vol_ratio = \
                set_proj_dim(R, e, ell, GSnorm, cur_lv, cur_r, C_proj)
            
            if next_r == 0:
                continue
            S *= vol_ratio

        if cur_lv == 0:
            # Bottom-level Optimization
            if t != 1:
                cur_r = set_bottom_proj_dim(S, next_r, GSnorm, e, b_lsh, ell, C_lsh)
                if cur_r == 0:
                    continue
            running_stat['r'].append(cur_r)

        lsh_domain = copy.deepcopy(GSnorm[-cur_r:])
        for i in range(next_r):
            lsh_domain[-i-1] = 2*ell
        
        log_cost_list = Log2(S * babai_cost(cur_r))
        log_cost_lsh, log_p_col, R_lsh, lsh_dim = lsh_cost(S, e, b_lsh, lsh_domain, cur_lv, C_lsh)
        
        if Log2(2**log_cost_list + 2**log_cost_lsh) > abort_bound:
            if verbose:
                print(" " * cur_lv, '- r', cur_lv, '=', cur_r, 
                ', w', cur_lv+1, '=', w_next, 
                ', cost_lsh', cur_lv+1, '->', cur_lv, '=', log_cost_lsh, 
                ', cost_list', cur_lv+1, '=', log_cost_list, sep ="")
                print(" " * cur_lv, '- R', cur_lv+1, '=2^', Log2(R), 
                ', (log_p_col=', log_p_col, ', R_lsh=2**', Log2(R_lsh), ', lsh_dim=', lsh_dim, ') => abort', sep ="")
            # if cur_lv < t-1:
            #     continue
            # else:
            #     break
            break
            
        running_stat['w'].append(w_next)
        if cur_lv < t-1:
            running_stat['r'].append(next_r)
            running_stat['ell'].append(ell)
        running_stat['b_lsh'].append(b_lsh)
        running_stat['R_lsh'].append(R_lsh)

        running_stat['log_cost_lsh'].append(log_cost_lsh)
        running_stat['log_cost_list'].append(log_cost_list)
        running_stat['log_fullS'].append(Log2(num_ternary_secret(d, w_next)))
        running_stat['log_S'].append(Log2(S))

        p_sp = 1.0; p_ncf = 1.0
        if cur_lv < t-1:
            if cur_lv == 0:
                if cur_r >= 800:
                    continue
                p_sp = p_jensen_hardcode_gauss[next_r // 100]
                p_ncf = p_jensen_hardcode_gauss[cur_r // 100]
                for i in range(cur_r - next_r):
                    p_ncf *= prob_admissible_gaussian(GSnorm[-cur_r+i], e)
            else:
                if cur_r >= 220:
                    continue
                p_sp = p_jensen_hardcode_unif[next_r // 20]
                p_ncf = p_jensen_hardcode_unif[cur_r // 20]
                for i in range(cur_r - next_r):
                    p_ncf *= prob_admissible_uniform(GSnorm[-cur_r+i], e)     
        else:
            p_sp = 0.95                    
            for i in range(cur_r):
                p_ncf *= prob_admissible_gaussian(GSnorm[i], e)


        running_stat['log_R'].append(Log2(R))
        running_stat['log_p_sp'].append(2**cur_lv * Log2(p_sp))
        running_stat['log_p_ncf'].append(2**cur_lv * Log2(p_ncf))
        running_stat['lsh_dim'].append(lsh_dim)
        running_stat['log_p_col'].append(log_p_col)
        
        if cur_lv < t-1:
            if verbose:
                print(" " * cur_lv, '- r', cur_lv, '=', cur_r, 
                ', w', cur_lv+1, '=', w_next, 
                ', cost_lsh', cur_lv+1, '->', cur_lv, '=', log_cost_lsh, 
                ', cost_list', cur_lv+1, '=', log_cost_list, sep ="")
                print(" " * cur_lv, '- R', cur_lv+1, '=2^', Log2(R), 
                ', (log_p_col=', log_p_col, ', R_lsh=2**', Log2(R_lsh), ', lsh_dim=', lsh_dim, ')', sep ="")

            comp = meet_LWE_cost(cur_lv+1, GSnorm, d,
                            next_r, w_next, ell, ell, 2*ell, C_proj, C_lsh,
                            t, copy.deepcopy(running_stat), abort_bound)
            if comp['cost'] < cur_best['cost']:
                cur_best = copy.deepcopy(comp)
            
        if cur_lv == t-1:
            total_cost = 0
            logprob = 0
            for i in range(1, t+1):
                total_cost += 2**running_stat['log_cost_lsh'][i]
                total_cost += 2**running_stat['log_cost_list'][i]
                logprob += running_stat['log_p_ncf'][i]
                logprob += running_stat['log_p_sp'][i]
            log_total_cost = Log2(total_cost)

            running_stat['log_p_suc'] = round(logprob, 2)
            running_stat['cost'] = round(log_total_cost - logprob, 2)

            if running_stat['cost'] < cur_best['cost']:
                cur_best = copy.deepcopy(running_stat)

            if verbose:
                print(" " * cur_lv, '- r', cur_lv, '=', cur_r, 
                ', w', cur_lv+1, '=', w_next, 
                ', cost_lsh', cur_lv+1, '->', cur_lv, '=', log_cost_lsh, 
                ', cost_list', cur_lv+1, '=', log_cost_list,
                ', total_cost=', round(log_total_cost - logprob, 2), sep ="")
                print(" " * cur_lv, '- R', cur_lv+1, '=2^', Log2(R), 
                ', (log_p_col=', log_p_col, ', R_lsh=2**', Log2(R_lsh), ', lsh_dim=', lsh_dim, ')', sep ="")
        
    # Fail to reduce to higher level -> Rerun this level with [top level = this level]
    # if cur_best['cost'] == np.inf and cur_lv < t-1:
    #     current_guess['top_reach_fail'] = [True]
    #     cur_best = meet_LWE_cost(cur_lv, GSnorm, d,
    #                             r, w, e, ell, b_lsh, C_proj, C_lsh,
    #                             cur_lv+1, copy.deepcopy(current_guess), abort_bound)

    # if cur_best['cost'] == np.inf and cur_lv < t-1:
    #     cur_best['top_reach_fail'] = True
        
    return cur_best

def set_proj_dim(R, e, ell, GSnorm, cur_lv, r_max, C_proj):    
    r = 0
    p_rep = 1
    vol_ratio = 1
    
    while True:
        cur_p_rep = 1
        cur_axis_length = min(GSnorm[-r-1], 2*ell)
        cur_ratio = cur_axis_length / GSnorm[-r-1]
        if cur_lv != 0:
            cur_p_rep = cur_ratio * prob_admissible_uniform(cur_axis_length, e)
        else:
            cur_p_rep = cur_ratio * prob_admissible_gaussian(cur_axis_length, e)

        if R * p_rep * cur_p_rep > 2 * C_proj:
            vol_ratio *= cur_ratio
            p_rep *= cur_p_rep
            r += 1
        else:
            break
        if r == r_max:
            break

    return r, p_rep, vol_ratio

def lsh_cost(S, e, b_lsh, lsh_domain, cur_lv, C_lsh):
    p_lsh = 1.0
    log_p_col = 0.0
    lsh_dim = 1
    r = len(lsh_domain)

    while True:
        logn = max(0, Log2(lsh_domain[lsh_dim - 1]) - Log2(b_lsh))
        log_p_col -= logn
        if logn > 0:
            if cur_lv != 0:
                p_lsh *= prob_admissible_uniform(lsh_domain[lsh_dim - 1], e)
            else:
                p_lsh *= prob_admissible_gaussian(lsh_domain[lsh_dim - 1], e)
        if Log2(S) + log_p_col < 1:
            break
        if lsh_dim == r:
            break
        lsh_dim += 1
    R_lsh = ceil(C_lsh / p_lsh)

    log_N_col = 2*Log2(S) + log_p_col - 1
    log_cost_col = Log2(R_lsh) + Log2(lsh_dim) + Log2(S + 2**log_N_col)

    return round(log_cost_col, 2), round(log_p_col, 2), R_lsh, lsh_dim

def set_bottom_proj_dim(S, r1, GSnorm, e, b_lsh, ell, C_lsh):
    m = len(GSnorm)
    r0 = 1
    p_lsh = 1.0
    log_p_col = 1.0

    while True:
        cur_axis_length = 2*ell
        if r0 >= r1:
            cur_axis_length = GSnorm[-r0]
        logn = max(0, Log2(cur_axis_length) - Log2(b_lsh))
        log_p_col -= logn
        if logn > 0:
            p_lsh *= prob_admissible_gaussian(cur_axis_length, e)
        R_lsh = ceil(C_lsh / p_lsh)

        cost_list = Log2(S * babai_cost(r0))
        log_N_col = 2*Log2(S) + log_p_col - 1
        cost_col = Log2(R_lsh) + Log2(r0) + log_N_col
        if log_p_col + Log2(S) <= 1:
            cost_col = Log2(R_lsh) + Log2(r0) + Log2(S)
        if r0 >= r1 and cost_list >= cost_col:
            break
        if r0 == m:
            break
        r0 += 1
    return r0