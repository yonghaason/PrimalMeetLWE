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

    nu = stddev / float(sqrt(w / n))
    if is_sec_bal:
        nu = stddev / float(sqrt(2*w / n))
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
    if pr_np == 0 or stddev > GSnorm[-1]/2:
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
    if 'w_g' in param_specified:
        w_g_range = range(param_specified['w_g'], param_specified['w_g'])

    for w_g in w_g_range:
        if w_g == 0:
            continue

        cur_guess = None
        if lv == 0: # Exhaustive Search
            cur_guess = {'cost': round(Log2(num_ternary_secret(zeta, w_g, is_sec_bal) * babai_cost(len(GSnorm))), 2)}
        else: # Our (Generalized) Meet-LWE
            initial_stat = \
            {'log_S_i': [], 'w_i': [w_g],'log_R_i': [], 'log_L_i': [1], 
            'log_cost_babai': [], 'log_cost_lsh': [], 'proj_dim': [d], 
            'log_p_rep': [], 'log_vol_ratio': [], 'cons_bound': [], 'lsh_stats': []}
            cur_guess = meet_LWE_cost(1, GSnorm, is_sec_bal, zeta, 
                proj_dim = len(GSnorm), w = w_g, is_error_unif=False, 
                error_param = stddev, lsh_length = 6*stddev,
                last_lv=lv, current_guess = copy.deepcopy(initial_stat), abort_bound=log_lat)
        
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
    
def meet_LWE_cost(lv, GSnorm, is_sec_bal, zeta,
    proj_dim, w, is_error_unif, error_param, 
    lsh_length, last_lv, current_guess, abort_bound = np.inf):
    
    d = len(GSnorm)
    cur_best = copy.deepcopy(current_guess)
    # Set upper bound by brute-force search of candidates
    
    if lv == last_lv:
        cur_best['cost'] = Log2(num_ternary_secret(zeta, w, is_sec_bal) * babai_cost(proj_dim))
        p_adm = p_admissible(GSnorm, range(d-proj_dim, d), is_error_unif, error_param)
        
        for w_split in range(w//2, w):
            running_stat = copy.deepcopy(current_guess)
            R = ambiguity(zeta, w, w_split, is_sec_bal)
            if R * p_adm < 10.0:
                continue
            S = num_ternary_secret(zeta, w_split, is_sec_bal)
            L = S / sqrt(R*p_adm)

            cost_lsh = None
            stat_lsh = None

            if lv == 1:
                cost_lsh, stat_lsh = lsh_cost_bottom(L, error_param, lsh_length, GSnorm,
                                        lsh_last_idx = d - 1,
                                        is_error_unif = is_error_unif)
                running_stat['proj_dim'].append(d - stat_lsh['lsh_coords'][0])
                S0 = num_ternary_secret(zeta, w, is_sec_bal)
                vol_ratio0 = 1
                for i in range(stat_lsh['lsh_coords'][0], stat_lsh['lsh_coords'][1]+1):
                    vol_ratio0 *= GSnorm[i] / lsh_length
                running_stat['log_L_i'] = [max(0, Log2(S0/vol_ratio0))]
                
            else:
                cost_lsh, stat_lsh = lsh_cost(L, error_param, lsh_length, GSnorm,
                                lsh_start_idx = d-proj_dim, 
                                is_error_unif = is_error_unif)
            
            log_cost_babai = running_stat['log_L_i'][-1] + Log2(babai_cost(running_stat['proj_dim'][-2]))
            cost_top_babai = L * babai_cost(running_stat['proj_dim'][-1])
            
            total_cost = cost_lsh + cost_top_babai + 2**log_cost_babai
            for log_lsh in running_stat['log_cost_lsh']:
                total_cost += 2**log_lsh
            for log_babai in running_stat['log_cost_babai']:
                total_cost += 2**log_babai
            log_total_cost = Log2(total_cost)

            running_stat['cost'] = log_total_cost
            running_stat['log_cost_lsh'].append(Log2(cost_lsh))
            running_stat['log_cost_babai'].append(round(log_cost_babai, 2))
            running_stat['log_cost_babai'].append(Log2(cost_top_babai))
            running_stat['log_S_i'].append(Log2(S))
            running_stat['log_L_i'].append(Log2(L))
            running_stat['log_R_i'].append(Log2(R))
            running_stat['w_i'].append(w_split)
            running_stat['log_p_adm'] = [Log2(p_adm)]
            running_stat['lsh_stats'].append(stat_lsh)

            if running_stat['cost'] < cur_best['cost']:
                cur_best = copy.deepcopy(running_stat)
        
    else:
        cur_best['cost'] = np.inf
        constraint_bound = error_param # One can further fine-tune this ...
        w_start = w//2
        if w_start % 2 == 1: w_start += 1
        
        for w_split in range(w_start, w, 2):
            running_stat = copy.deepcopy(current_guess)
            R = ambiguity(zeta, w, w_split, is_sec_bal)
            next_proj_dim, p_rep, vol_ratio = \
                constraint_area(R, error_param, constraint_bound, GSnorm,
                                is_error_unif = is_error_unif, 
                                proj_last = proj_dim)

            if next_proj_dim == 0:
                continue
            S = num_ternary_secret(zeta, w_split, is_sec_bal)
            L = S / vol_ratio
                
            cost_lsh = None
            stat_lsh = None

            if lv == 1:
                cost_lsh, stat_lsh = lsh_cost_bottom(L, error_param, lsh_length, GSnorm,
                                        lsh_last_idx = d - next_proj_dim,
                                        is_error_unif = is_error_unif)
                running_stat['proj_dim'].append(d - stat_lsh['lsh_coords'][0])
                S0 = num_ternary_secret(zeta, w, is_sec_bal)
                vol_ratio0 = 1
                for i in range(stat_lsh['lsh_coords'][0], stat_lsh['lsh_coords'][1]+1):
                    vol_ratio0 *= GSnorm[i] / lsh_length
                running_stat['log_L_i'] = [max(0, Log2(S0/vol_ratio0))]
                
            else:
                cost_lsh, stat_lsh = lsh_cost(L, error_param, lsh_length, GSnorm,
                                lsh_start_idx = d - proj_dim,
                                constraint_cube_length = 2*constraint_bound,
                                constraint_dim_idx = d - next_proj_dim, 
                                is_error_unif = is_error_unif)
            
            log_cost_babai = running_stat['log_L_i'][-1] + Log2(babai_cost(running_stat['proj_dim'][-2]))
            
            if Log2(cost_lsh) > abort_bound:
                continue

            running_stat['log_cost_lsh'].append(Log2(cost_lsh))
            running_stat['log_cost_babai'].append(round(log_cost_babai, 2))
             
            running_stat['log_S_i'].append(Log2(S))
            running_stat['log_L_i'].append(Log2(L))
            running_stat['log_R_i'].append(Log2(R))
            running_stat['w_i'].append(w_split)
            running_stat['log_p_rep'].append(Log2(p_rep))
            running_stat['log_vol_ratio'].append(Log2(vol_ratio))
            running_stat['proj_dim'].append(next_proj_dim)
            running_stat['cons_bound'].append(error_param)
            running_stat['lsh_stats'].append(stat_lsh)

            comp = meet_LWE_cost(lv+1, GSnorm, is_sec_bal, zeta,
            next_proj_dim, w_split, True, constraint_bound, 2*constraint_bound, 
            last_lv, copy.deepcopy(running_stat), abort_bound)

            if comp['cost'] < cur_best['cost']:
                cur_best = copy.deepcopy(comp)

        if cur_best['cost'] == np.inf:
            # Fail to reduce to higher level, so set the current level as the top level
            cur_best = meet_LWE_cost(lv, GSnorm, is_sec_bal, zeta,
                proj_dim, w, True, constraint_bound, 2*constraint_bound, 
                lv, copy.deepcopy(current_guess), abort_bound)
    
    return cur_best

def constraint_area(R, error_param, constraint_bound, GSnorm, is_error_unif, proj_last = None):
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
        if is_error_unif:
            cur_p_adm = prob_admissible_uniform(cur_axis_length, error_param)
        else:
            cur_p_adm = prob_admissible_gaussian(cur_axis_length, error_param)

        inv = 10.0 / (p_adm * cur_p_adm / cur_ratio)
        if R < inv:
            proj_dim -= 1
            break
        else:
            vol_ratio *= cur_ratio
            p_adm *= cur_p_adm / cur_ratio
            proj_dim += 1

    return proj_dim, p_adm, vol_ratio

def prob_lsh(L, error_param, lsh_length, GSnorm, lsh_start_idx, \
             constraint_cube_length = 0, constraint_dim_idx = np.inf, is_error_unif = True):
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
        axis_len = GSnorm[lsh_last_idx]
        cur_p_torus_good = 1
        cur_p_torus_bad = 1
        if lsh_last_idx >= constraint_dim_idx:
            axis_len = constraint_cube_length
        if axis_len > lsh_length:
            if is_error_unif:
                cur_p_torus_good = prob_admissible_uniform(lsh_length, error_param)
            else:
                cur_p_torus_good = prob_admissible_gaussian(lsh_length, error_param)
            cur_p_torus_bad = lsh_length / axis_len
        
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
    values['lsh_coords'] = [lsh_start_idx, lsh_last_idx]
    try:
        values['log_p_good'] = Log2(p_torus_good)
    except ValueError:
        values['log_p_good'] = -np.inf
    try:
        values['log_p_bad'] = Log2(p_torus_bad)
    except ValueError:
        values['log_p_bad'] = -np.inf
    
    return values

def lsh_cost(L, error_param, lsh_length, GSnorm, lsh_start_idx, \
             constraint_cube_length = 0, constraint_dim_idx = np.inf, is_error_unif=True):
    values = prob_lsh(L, error_param, lsh_length, GSnorm, lsh_start_idx, \
                      constraint_cube_length, constraint_dim_idx, is_error_unif).copy()
    log_num_collision = 2*Log2(L) + values['log_p_bad'] - values['log_p_good']
    log_cost = max(Log2(L), log_num_collision)
    try:
        return 2**log_cost, values
    except OverflowError:
        return np.inf, values

def lsh_cost_bottom(L, error_param, lsh_length, GSnorm, lsh_last_idx, is_error_unif):
    lsh_start_idx = lsh_last_idx
    p_torus_good = 1
    p_torus_bad = 1

    while True:
        cur_p_torus_good = 1
        cur_p_torus_bad = 1
        if GSnorm[lsh_start_idx] > lsh_length:
            if is_error_unif:
                cur_p_torus_good = prob_admissible_uniform(lsh_length, error_param)
            else:
                cur_p_torus_good = prob_admissible_gaussian(lsh_length, error_param)
            cur_p_torus_bad = lsh_length / GSnorm[lsh_start_idx]
        
        if L**2 * p_torus_bad * cur_p_torus_bad < 1 :
            if lsh_last_idx != lsh_start_idx:
                lsh_start_idx += 1
            else:
                p_torus_good = cur_p_torus_good
                p_torus_bad = cur_p_torus_bad
            break
        else:
            p_torus_good *= cur_p_torus_good
            p_torus_bad *= cur_p_torus_bad
            if lsh_start_idx == 0:
                break
            lsh_start_idx -= 1

    values = {}
    values['box_length'] = lsh_length
    values['lsh_coords'] = [lsh_start_idx, lsh_last_idx]
    values['log_p_good'] = Log2(p_torus_good)
    values['log_p_bad'] = Log2(p_torus_bad)
    
    log_num_collision = 2*Log2(L) + values['log_p_bad'] - values['log_p_good']
    log_cost = max(Log2(L), log_num_collision)
    try:
        return 2**log_cost, values
    except OverflowError:
        return np.inf, values
