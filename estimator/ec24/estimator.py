from utils import *
from functools import partial
import copy

def primal_may(param, lv, is_secret_balanced = False, filename=False, param_specified = {}):
    
    if filename is not False:
        f = open(filename, 'a')
        sys.stdout = f

    n = param[0]
    q = param[1]
    error_type = param[2]
    error_param = param[3]
    w = param[4]
    m = param[5]
    if len(param) > 6:
        is_secret_balanced = param[6]

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

    return primal_may_inner(n, q, stddev, is_secret_balanced, w, m, lv=lv, param_specified = param_specified)

def primal_may_inner(n, q, stddev, is_sec_bal, w, m, lv, param_specified):

    # nu = stddev / float(sqrt(w / n))
    # if is_sec_bal:
    #     nu = stddev / float(sqrt(2*w / n))
    nu = 1
    func = partial(cost_zeta, n=n, q=q, stddev=stddev, 
                    is_sec_bal=is_sec_bal, w=w, m=m, nu=nu, lv=lv, 
                    param_specified = param_specified)
    best = None
    if 'zeta' in param_specified:
        zeta = param_specified['zeta']
        best = my_binary_search(zeta, zeta, func, cur_depth = 0, max_depth = 6)
    else:
        best = my_binary_search(0, n, func, cur_depth = 0, max_depth = 6)
    
    print("----------------------")

    prettyprint(best)

    return best

def cost_zeta(zeta, n, q, stddev, is_sec_bal, w, m, nu, lv, param_specified):

    m_ = m + n + 1 - zeta
    d1 = n + 1 - zeta

    probs_hw = probs_hw_precompute(n, zeta, w, is_sec_bal)  

    line = ' Loop with zeta = %d' % (zeta) 
    print(line)     

    # Find the best submatrix size d < m_
    func = partial(cost_d, zeta=zeta, n=n, q=q, stddev=stddev, 
                    is_sec_bal=is_sec_bal, w=w, nu=nu, lv=lv, probs_hw = probs_hw, param_specified = param_specified)
    if 'd' in param_specified:
        d = param_specified['d']
        best = my_binary_search(d, d, func, cur_depth = 0, max_depth = 6)
    else:
        best = my_binary_search(d1 + 1, m_, func, cur_depth = 0, max_depth = 5)
    
    line = ' └──── Best: %f' % best['cost']
    print(line)

    return best

def cost_d(d, zeta, n, q, stddev, is_sec_bal, w, nu, lv, param_specified, probs_hw = []):

    if d < 40:
        return {'cost': np.inf}

    beta_max = min(d, 1000)

    # Find the best blocksize beta
    func = partial(cost_beta, zeta=zeta, n=n, q=q, d=d, stddev=stddev, 
                    is_sec_bal=is_sec_bal, w=w, nu=nu, lv=lv, probs_hw = probs_hw, param_specified = param_specified)
    if 'beta' in param_specified:
        beta = param_specified['beta']
        best = my_binary_search(beta, beta, func, cur_depth = 0, max_depth = 4)
    else:
        best = my_binary_search(40, beta_max, func, cur_depth = 0, max_depth = 4)
    
    print(' │ d = %d:' % d, best['cost'])
    return best

def cost_beta(beta, zeta, n, q, d, stddev, is_sec_bal, w, nu, lv, param_specified, probs_hw = []):

    log_lat = log_BKZ_cost(d, beta)
    GSnorm = GSA(d, q, d,  n + 1 - zeta, beta, nu)

    # print('         zeta = %d / d = %d / beta = %d :' % (zeta, d, beta))    

    lat = {'cost': round(log_lat, 2)}
    lat['zeta'] = zeta
    lat['d'] = d
    lat['beta'] = beta
    lat['nu'] = round(nu, 3)
    lat['last_R'] = round(GSnorm[-1], 2)

    pr_np = prob_np(GSnorm, stddev)

    # Skip if the last GSnorm is too short
    if pr_np == 0 or 6*stddev*(3**(lv-1)) > GSnorm[-1]/4:
        lat['log_pr_np'] = np.inf
        lat['cost'] = np.inf
        return lat
    log_pr_np = Log2(pr_np)
    lat['log_pr_np'] = log_pr_np

    # rough bound assuming |S|^0.25 < lat['cost']
    w_g_bound = ceil(4*lat['cost']/(1 + Log2(zeta)))
    w_g_min = max(w - (n-zeta), 0)
    w_g_max = min(w_g_bound, zeta, w)

    if is_sec_bal:
        w_g_min = max((w - (n-zeta)//2), 0)
        w_g_max = min(w_g_bound//2, zeta//2, w)

    if w_g_min > w_g_max:
        return {'cost': np.inf}

    cost = {'lat': lat.copy()}

    guesses = []
    guess_cost = 0.0
    pr_hw = 0.0
    w_g_range = range(((w_g_min+1)//2)*2, (w_g_max//2)*2, 2) if (lv >= 1) \
                else range(w_g_min, w_g_max)

    print('log_lat:', log_lat)

    if 'w_g' in param_specified:
        w_g_range = range(param_specified['w_g'], param_specified['w_g'])

    # for w_g in w_g_range:
    for w_g in range(10, 20, 2):
        if w_g == 0:
            continue

        cur_guess = None
        if lv == 0: # Exhaustive Search
            cur_guess = {'cost': round(Log2(num_ternary_secret(zeta, w_g, is_sec_bal) * babai_cost(len(GSnorm))), 2)}
        else: # Our (Generalized) Meet-LWE
            initial_stat = \
            {'cost': np.inf,
            'w_i': [w_g], 
            'proj_dim': [], 
            'norm_bound': [6*stddev],
            'log_cost_list': ['-'], 
            'log_cost_lsh': ['-'],  
            'log_S_i': ['-'], 
            'log_R_i': ['-'], 
            'log_p_sp': ['-'], 
            'log_p_lsh': ['-']}
            cur_guess = meet_LWE_cost(0, GSnorm, is_sec_bal, zeta, 
                proj_dim = len(GSnorm), w = w_g, error_param = stddev,
                norm_bound = 6*stddev, lsh_length = 9*stddev,
                top_lv=lv, current_guess = copy.deepcopy(initial_stat), abort_bound=log_lat)
        
        if cur_guess['cost'] is np.inf or cur_guess['cost'] > log_lat:
            break
        else:
            if w_g >= len(probs_hw):
                pr_hw += prob_hw(n, zeta, w, w_g, is_sec_bal)
            else:
                pr_hw += probs_hw[w_g]
            guesses.append(cur_guess)
            guess_cost += 2**cur_guess['cost']
            cost['w_guess'] = w_g
            
    cost['guess_cost'] = Log2(guess_cost)
    cost['log_p_hw'] = Log2(pr_hw)    
    cost['cost'] = round(max(log_lat, Log2(guess_cost)) - Log2(pr_hw) - log_pr_np, 2)
    cost['guess_details'] = guesses

    return cost
    
def meet_LWE_cost(cur_lv, GSnorm, is_sec_bal, zeta,
    proj_dim, w, error_param, norm_bound, lsh_length, 
    top_lv, current_guess, abort_bound = np.inf):
    
    d = len(GSnorm)
    cur_best = copy.deepcopy(current_guess)

    # trivial_bound = Log2(num_ternary_secret(zeta, w, is_sec_bal) * babai_cost(proj_dim))

    # if cur_lv == top_lv:
    #     cur_best['cost'] = trivial_bound

    w_start = w//2
    w_range = None

    if cur_lv == top_lv-1:
        w_range = range(w_start, w)
    else:
        if w_start % 2 == 1: 
            w_start += 1
        w_range = range(w_start, w, 2)

    if cur_lv == 0: print()
    print(" " * cur_lv, "LV = ",  cur_lv, " (t=", top_lv, "), w", cur_lv, "=", w, ", w_range=", w_range, sep ="")
    
    for w_split in w_range:
        running_stat = copy.deepcopy(current_guess)

        R = ambiguity(zeta, w, w_split, is_sec_bal)
        S = num_ternary_secret(zeta, w_split, is_sec_bal)
        next_proj_dim = 0

        if cur_lv == top_lv-1:
            # Top-level Optimization
            S *= 3.0 / sqrt(R)
        else:
            next_proj_dim, p_rep, vol_ratio = \
                set_proj_dim(R, error_param, norm_bound, GSnorm, cur_lv, proj_last = proj_dim)
            if next_proj_dim == 0:
                continue
            S *= vol_ratio

        lsh_domain = GSnorm      
        for i in range(next_proj_dim):
            lsh_domain[-i-1] = 2*norm_bound

        if cur_lv == 0:
            # Bottom-level Optimization
            proj_dim = set_bottom_proj_dim(S, next_proj_dim, lsh_domain, error_param, lsh_length)
            running_stat['proj_dim'].append(proj_dim)
        
        cost_list = S * babai_cost(proj_dim)
        cost_lsh, p_col, p_lsh = lsh_cost(S, error_param, lsh_length, lsh_domain, cur_lv)
        
        if Log2(cost_lsh + cost_list) > abort_bound:
            print(" " * cur_lv, '- r', cur_lv, '=', proj_dim, 
            ', w', cur_lv+1, '=', w_split, 
            ', cost_lsh', cur_lv+1, '->', cur_lv, '=', Log2(cost_lsh), 
            ', cost_list', cur_lv+1, '=', Log2(cost_list), sep ="")
            print(" " * cur_lv, '- R', cur_lv+1, '=2^', Log2(R), ', (p_col=', p_col, ', p_lsh=', p_lsh, ') => abort', sep ="")
            continue

        running_stat['w_i'].append(w_split)
        if cur_lv < top_lv:
            running_stat['proj_dim'].append(next_proj_dim)
            running_stat['norm_bound'].append(3*norm_bound)

        running_stat['log_cost_lsh'].append(Log2(cost_lsh))
        running_stat['log_cost_list'].append(Log2(cost_list))
        running_stat['log_S_i'].append(Log2(S))
        
        running_stat['log_R_i'].append(Log2(R))
        running_stat['log_p_sp'].append(2**cur_lv * Log2(0.85))
        running_stat['log_p_lsh'].append(2**(cur_lv-1) * Log2(0.85))

        if cur_lv < top_lv-1:
            print(" " * cur_lv, '- r', cur_lv, '=', proj_dim, 
            ', w', cur_lv+1, '=', w_split, 
            ', cost_lsh', cur_lv+1, '->', cur_lv, '=', Log2(cost_lsh), 
            ', cost_list', cur_lv+1, '->', cur_lv, '=', Log2(cost_list), sep ="")
            print(" " * cur_lv, '- R', cur_lv+1, '=2^', Log2(R), ', (p_col=', p_col, ', p_lsh=', p_lsh, ')', sep ="")

            comp = meet_LWE_cost(cur_lv+1, GSnorm, is_sec_bal, zeta,
                            next_proj_dim, w_split, norm_bound, 3*norm_bound, 3*lsh_length, 
                            top_lv, copy.deepcopy(running_stat), abort_bound)
            if comp['cost'] < cur_best['cost']:
                cur_best = copy.deepcopy(comp)
            
        if cur_lv == top_lv-1:            
            total_cost = 0
            logprob = 0
            for i in range(1, top_lv+1):
                total_cost += 2**running_stat['log_cost_lsh'][i]
                total_cost += 2**running_stat['log_cost_list'][i]
                logprob += running_stat['log_p_lsh'][i]
                logprob += running_stat['log_p_sp'][i]
            log_total_cost = Log2(total_cost)

            running_stat['prob'] = logprob
            running_stat['cost'] = log_total_cost - logprob

            if running_stat['cost'] < cur_best['cost']:
                cur_best = copy.deepcopy(running_stat)

            print(" " * cur_lv, '- r', cur_lv, '=', proj_dim, 
            ', w', cur_lv+1, '=', w_split, 
            ', cost_lsh', cur_lv+1, '->', cur_lv, '=', Log2(cost_lsh), 
            ', cost_list', cur_lv+1, '->', cur_lv, '=', Log2(cost_list),
            ', total_cost=', round(log_total_cost - logprob, 2), sep ="")
            print(" " * cur_lv, '- R', cur_lv+1, '=2^', Log2(R), ', (p_col=', p_col, ', p_lsh=', p_lsh, ')', sep ="")
        
    if cur_best['cost'] == np.inf and cur_lv < top_lv-1:
        # Fail to reduce to higher level -> Rerun this level with [top level = this level]
        current_guess['level_fail'] = [True]
        cur_best = meet_LWE_cost(cur_lv, GSnorm, is_sec_bal, zeta,
                                proj_dim, w, error_param, norm_bound, lsh_length,
                                cur_lv+1, copy.deepcopy(current_guess), abort_bound)

    # if cur_best['cost'] == trivial_bound:
    #     current_guess['level_fail'] = [True]
        
    return cur_best

def set_proj_dim(R, error_param, norm_bound, GSnorm, cur_lv, proj_last = None):    
    proj_dim = 0
    p_rep = 1
    vol_ratio = 1
    
    while True:
        cur_p_rep = 1
        cur_axis_length = min(GSnorm[-proj_dim-1], 2*norm_bound)
        cur_ratio = cur_axis_length / GSnorm[-proj_dim-1]
        if cur_lv != 0:
            cur_p_rep = cur_ratio * prob_admissible_uniform(cur_axis_length, error_param)
        else:
            cur_p_rep = cur_ratio * prob_admissible_gaussian(cur_axis_length, error_param)

        if R * p_rep * cur_p_rep > 2 * 3.0:
            vol_ratio *= cur_ratio
            p_rep *= cur_p_rep
            proj_dim += 1
        else:
            break
        if proj_dim == proj_last:
            break

    return proj_dim, p_rep, vol_ratio

def lsh_cost(S, error_param, lsh_length, lsh_domain, cur_lv):
    p_lsh = 1.0
    p_col = 1.0
    r = len(lsh_domain)

    for i in range(len(lsh_domain)):
        n = max(1, floor(lsh_domain[i] / lsh_length))
        p_col /= n
        if n > 1:
            if cur_lv != 0:
                p_lsh *= prob_admissible_uniform(lsh_domain[i], error_param)
            else:
                p_lsh *= prob_admissible_gaussian(lsh_domain[i], error_param)
    R_lsh = ceil(3.0 / p_lsh)

    return R_lsh*(S + 0.5*(S**2)*p_col)*r, p_col, p_lsh

def set_bottom_proj_dim(S, start_dim, lsh_domain, error_param, lsh_length):
    m = len(lsh_domain)
    proj_dim = 0
    p_lsh = 1.0
    p_col = 1.0

    for i in range(1, m+1):
        n = max(1, floor(lsh_domain[-i] / lsh_length))
        p_col /= n
        if n > 1:
            p_lsh *= prob_admissible_gaussian(lsh_domain[-i], error_param)
        cost_list = S * babai_cost(i)
        R_lsh = ceil(3.0 / p_lsh)
        cost_lsh = R_lsh*(S + 0.5*(S**2)*p_col)*i
        if i >= start_dim and cost_list > cost_lsh:
            proj_dim = i
            break
    return proj_dim