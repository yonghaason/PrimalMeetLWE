#include "functions.h"
#include <iostream>
#include <random>
#include <algorithm>
#include <map>

using namespace std;

UniformTest::UniformTest(
  uint32_t n, double e, uint64_t q,
  uint32_t h, uint32_t m, bool unif)
  : n(n), e(e), h(h), m(m), unif(unif)
{
  random_device rd;
  mt19937 gen(rd());
  vector<uniform_real_distribution<double>> coord_sampler(m);

  GSnorm.resize(m);
  for (size_t i = 0; i < m; i++) {
    GSnorm[i] = pow(q, (1.0 - (double) i/(2*m)));
  }

  B.resize(m);
  for (size_t i = 0; i < m; i++) {
    B[i].resize(m);
  }
  for (size_t j = 0; j < m; j++) {
    B[j][j] = GSnorm[j];
    coord_sampler[j] = uniform_real_distribution<double>(-GSnorm[j]/2, GSnorm[j]/2);
    for (size_t i = 0; i < j; i++) {
      B[i][j] = coord_sampler[j](gen);
    }
  }

  cout << "- Dim n = " << n << endl;
  cout << "- Coordinate Length: [" << GSnorm[0] << ", ... , " 
      << GSnorm[m-1] << "] (GSA)" << endl;
  cout << "- # Sample m = " << m << endl;
  cout << "- Error distribution : ";
  if (unif) cout << "Uniform over [-" << e << ", " << e << "]" << endl;
  else cout << "Gaussian of stddev = " << e << endl;
  cout << endl;
}

vector<pair<secret, vector<double>>> UniformTest::sparse_secret_list(
  matrix& Mtrans, uint32_t d, uint32_t w) 
{
  vector<pair<secret, vector<double>>> R;
  if (d == 0 | d < w) {
    return R;
  }
  else if (w == 0) {
    R.push_back(make_pair(vector<int64_t>(d), vector<double>(m)));
    return R;
  }
  else if (d == 1 && w == 1) {
    R.push_back(make_pair(vector<int64_t>{1}, Mtrans[d-1]));
    vector<double> zerovec(m);
    R.push_back(make_pair(vector<int64_t>{-1}, subvec(zerovec, Mtrans[d-1])));
    return R;
  }
  else {
    auto R1 = sparse_secret_list(Mtrans, d-1, w-1);
    auto R2 = sparse_secret_list(Mtrans, d-1, w);

    for (size_t i = 0; i < R1.size(); i++) {
      auto v = R1[i].first;
      v.push_back(1);
      R.push_back(make_pair(v, addvec(R1[i].second, Mtrans[d-1])));
      v[d-1] = -1;
      R.push_back(make_pair(v, subvec(R1[i].second, Mtrans[d-1])));
    }
    for (size_t i = 0; i < R2.size(); i++) {
      auto v = R2[i].first;
      v.push_back(0);
      R.push_back(make_pair(v, R2[i].second));
    }
    return R;
  }
}

vector<secret> UniformTest::enumerate_secrets(uint32_t d, uint32_t w)
{
  vector<secret> R;
  if (d == 0 | d < w) {
    return R;
  }
  else if (w == 0) {
    R.push_back(vector<int64_t>(d));
    return R;
  }
  else if (d == 1 && w == 1) {
    R.push_back(vector<int64_t>{1});
    R.push_back(vector<int64_t>{-1});
    return R;
  }
  else {
    auto R1 = enumerate_secrets(d-1, w-1);
    auto R2 = enumerate_secrets(d-1, w);

    for (size_t i = 0; i < R1.size(); i++) {
      auto v = R1[i];
      v.push_back(1);
      R.push_back(v);
      v[d-1] = -1;
      R.push_back(v);
    }
    for (size_t i = 0; i < R2.size(); i++) {
      auto v = R2[i];
      v.push_back(0);
      R.push_back(v);
    }
    return R;
  }
}

uint32_t ConstraintTest::new_check_near_collision(
  matrix& M, secret& s, vector<double>& e, 
  uint32_t guess_weight, double constraint_bound) 
{
  matrix Mtrans = transpose(M);
  uint32_t count = 0;
  vector<double> Ms1(m);
  vector<int32_t> nonzeroidx(guess_weight);
  vector<int32_t> sign(guess_weight);
  recur_check(count, Mtrans, e, s, constraint_bound, 
              nonzeroidx, sign, Ms1, 0, n-1, 0, guess_weight);
  return count;
}

void ConstraintTest::recur_check(
  uint32_t& count, matrix& Mtrans, vector<double>& e, secret& s, double& constraint_bound, 
  vector<int32_t>& nonzeroidx, vector<int32_t>& sign, vector<double>& Ms1, 
  uint32_t start, uint32_t end, uint32_t cur_hw, uint32_t& k)
{
  if (cur_hw == k) {
    secret s1(n);
    for (size_t i = 0; i < nonzeroidx.size(); i++) {
      if (sign[i] == 1) s1[nonzeroidx[i]] = 1;
      else s1[nonzeroidx[i]] = -1;
    }
    auto s2 = sub(s1, s);
    if (weight_ternary_check(s2, k)) {
      Ms1 = babaiNP(Ms1, B);
      if (inf_norm(Ms1, m-proj_dim) <= constraint_bound) {
        auto Ms2 = subvec(Ms1, e);
        Ms2 = babaiNP(Ms2, B);
        if (inf_norm(Ms2, m-proj_dim) <= constraint_bound) {
          count++;
        }
      }
    }
  }
  else {
    vector<double> Ms1_next(m);
    for (int i = start; i <= end && end - i + 1 >= k - cur_hw; i++) {
        nonzeroidx[cur_hw] = i;
        sign[cur_hw] = 1;
        Ms1_next = addvec(Ms1, Mtrans[i]);
        recur_check(count, Mtrans, e, s, constraint_bound,
                    nonzeroidx, sign, Ms1_next, i+1, end, cur_hw+1, k);
        sign[cur_hw] = -1;  
        Ms1_next = subvec(Ms1, Mtrans[i]);
        recur_check(count, Mtrans, e, s, constraint_bound,
                    nonzeroidx, sign, Ms1_next, i+1, end, cur_hw+1, k);
    }
  }
}

void UniformTest::gen_noisy_instance(
  matrix& M, vector<int64_t>& s, vector<double>& error)
{
  random_device rd;
  mt19937 gen(rd());
  vector<uniform_real_distribution<double>> coord_sampler(m);
  
  matrix A;
  A.resize(m);
  for (size_t i = 0 ; i < m; i++) {
    A[i].resize(n-1);
    coord_sampler[i] = uniform_real_distribution<double>(-GSnorm[i]/2, GSnorm[i]/2);
    for (size_t j = 0; j < n-1; j++) {
      A[i][j] = coord_sampler[i](gen);
    }
  }

  s.resize(n-1);
  uniform_int_distribution<int64_t> binrand(0, 1);
  for (size_t i = 0; i < h-1; i++) {
    s[i] = 2*binrand(gen)-1;
  }
  for (size_t i = h-1; i < n-1; i++) {
    s[i] = 0;
  }
  shuffle(begin(s), end(s), gen);
  s.push_back(1);

  vector<double> b(m);
  for (size_t k = 0; k < n-1; k++) {
    if (s[k] != 0) {
      for (size_t i = 0; i < m; i++) {
        b[i] += s[k] * A[i][k];
      }
    }
  }

  error.resize(m);
  normal_distribution<double> gaussian_noise(0, e);
  uniform_real_distribution<double> uniform_noise(-e, e);
  if (unif) {
    for (size_t k = 0; k < m; k++) {
      error[k] = uniform_noise(gen);
      b[k] += error[k];
    }
  }
  else {
    for (size_t k = 0; k < m; k++) {
      error[k] = gaussian_noise(gen);
      b[k] += error[k];
    }
  }
  b = babaiNP(b, B);
  
  M.resize(m);
  for (size_t i = 0; i < m; i++) {
    M[i].resize(n);
    for (size_t j = 0; j < n-1; j++) {
      M[i][j] = -A[i][j];
    }
    M[i][n-1] = b[i];
  }
};

void UniformTest::gen_noisy_instance_fixed_error(
  matrix& M, vector<int64_t>& s, const vector<double>& error)
{
  random_device rd;
  mt19937 gen(rd());
  vector<uniform_real_distribution<double>> coord_sampler(m);
  
  matrix A;
  A.resize(m);
  for (size_t i = 0 ; i < m; i++) 
  {
    A[i].resize(n-1);
    coord_sampler[i] = uniform_real_distribution<double>(-GSnorm[i]/2, GSnorm[i]/2);
    for (size_t j = 0; j < n-1; j++) 
    {
      A[i][j] = coord_sampler[i](gen);
    }
  }

  s.resize(n-1);
  uniform_int_distribution<int64_t> binrand(0, 1);
  for (size_t i = 0; i < h-1; i++) 
  {
    s[i] = 2*binrand(gen)-1;
  }
  for (size_t i = h-1; i < n-1; i++) 
  {
    s[i] = 0;
  }
  shuffle(begin(s), end(s), gen);
  s.push_back(1);

  vector<double> b(m);
  for (size_t k = 0; k < n-1; k++) 
  {
    if (s[k] != 0) {
      for (size_t i = 0; i < m; i++) 
      {
        b[i] += s[k] * A[i][k];
      }
    }
  }

  for (size_t k = 0; k < m; k++) 
  {
    b[k] += error[k];
  }
  b = babaiNP(b, B);

  M.resize(m);
  for (size_t i = 0; i < m; i++) 
  {
    M[i].resize(n);
    for (size_t j = 0; j < n-1; j++) 
    {
      M[i][j] = -A[i][j];
    }
    M[i][n-1] = b[i];
  }
};


vector<pair<secret, vector<double>>> ConstraintTest::build_list(
  matrix& M, uint32_t guess_weight, double constraint_bound)
{
  // S = {(v, Mv): HW(v) = w}
  matrix Mtrans = transpose(M);
  auto S = sparse_secret_list(Mtrans, n, guess_weight);
  
  vector<pair<secret, vector<double>>> L;
  for (auto& p: S) {
    auto v = p.first;
    auto Mv = p.second;
    Mv = babaiNP(Mv, B);

    if (inf_norm(Mv, m-proj_dim) <= constraint_bound)
      L.push_back(make_pair(v, Mv));
  }
  return L;
};

void ConstraintTest::set_constraint_dim(uint32_t guess_weight, double constraint_bound, uint32_t target_pair_num)
{
  auto R = ambiguity(n, h, guess_weight);
  proj_dim = 1;
  vol_ratio = 1;
  p_rep_avg = 1;
  p_rep_wst = 1.0;
  auto m = GSnorm.size();

  double error_bound = unif ? e : 3*e;

  while (true) {
    if (proj_dim == m) break;
    
    double cur_ratio = 1;
    double cur_axis_length = (2*constraint_bound < GSnorm[m-proj_dim]) ?
                              2*constraint_bound : GSnorm[m-proj_dim];
    double cur_p_rep_avg = unif? 
                      prob_admissible_uniform(e, cur_axis_length): 
                      prob_admissible_gaussian(e, cur_axis_length);
    double cur_p_rep_wst = (1.0 - error_bound / cur_axis_length);

    if (2*constraint_bound < GSnorm[m-proj_dim]) {
      cur_ratio = GSnorm[m-proj_dim] / (2*constraint_bound);
    }

    auto inv = (double) target_pair_num / (p_rep_avg * cur_p_rep_avg / cur_ratio);
    if (R < inv) {
      proj_dim -= 1;
      break;
    }
    else {
      p_rep_avg *= cur_p_rep_avg / cur_ratio;
      p_rep_wst *= cur_p_rep_wst / cur_ratio;
      vol_ratio *= cur_ratio;
      proj_dim += 1;
    }
  }

  cout << "**** Theory & Expectation ****"  << endl;
  cout << "- # of reps (Ambiguity R): " << R << endl;
  cout << "- Projection dim (r): " << proj_dim << endl;
  cout << "- # pairs = " << R * p_rep_avg << " (R*p_rep_avg)" << endl;
  cout << "- p_rep_avg (log): " << p_rep_avg << " (" << log2(p_rep_avg) << ")" << endl;  
  cout << "- Pr[#pairs = 0] = " << pow(1 - p_rep_avg, R/2) << endl;
  // cout << "- p_rep_wst (log): " << p_rep_wst << " (" << log2(p_rep_wst) << ")" << endl;
  // cout << "- Pr[#pairs = 0] < " << pow(1 - p_rep_wst, R/2) << endl;  
  cout << "- vol_ratio (log): " << vol_ratio << " (" << log2(vol_ratio) << ")" << endl;

  cout << endl;
}

set<pair<secret, secret>> ConstraintTest::fast_check_near_collision(
  matrix& M, secret& s, vector<double>& e, 
  vector<secret>& partial_list, uint32_t guess_weight, double constraint_bound) 
{
  matrix Mtrans = transpose(M);
  set<pair<secret, secret>> sol;
  for (auto& s1: partial_list) {
    auto s2 = sub(s1, s);
    if (weight_ternary_check(s2, guess_weight)) {
      auto Ms1 = matmul(Mtrans, s1);
      auto Ms2 = subvec(Ms1, e);
      Ms1 = babaiNP(Ms1, B);
      Ms2 = babaiNP(Ms2, B);
      if ((inf_norm(Ms1, m-proj_dim) <= constraint_bound) && (inf_norm(Ms2, m-proj_dim) <= constraint_bound)) {
        sol.insert(make_pair(s1, s2));
      }
    }
  }
  return sol;
}

set<pair<secret, secret>> ConstraintTest::check_near_collision(
  vector<pair<secret, vector<double>>>& L, secret& s, uint32_t guess_weight)
{
  set<pair<secret, secret>> sol;
  set<secret> secrets;
  for (auto& p: L) {
    secrets.insert(p.first);
  }

  // If s1 in S is a valid parital key s.t. s = s1 - s2, 
  // s2 = s1 - s should be also in S
  for (auto& p: L) {
    secret s1 = p.first;
    secret s2 = sub(s1, s);
    if (weight_ternary_check(s2, guess_weight)) {
      if (secrets.find(s2) != secrets.end()) {
        sol.insert(make_pair(s1, s2));
      }
    }
  }
  return sol;
}

NearCollisionTest::NearCollisionTest(uint32_t r, uint64_t q, double e, bool unif, uint32_t lsh_dim, double lsh_length, uint32_t iter_mult)
  : r(r), e(e), unif(unif), lsh_dim(lsh_dim), lsh_length(lsh_length), iter_mult(iter_mult)
{
  dom.resize(r);
  for (size_t i = 0; i < r; i++) {
    dom[i] = pow(q, (1.0 - (double) i/(2*r)));
    // dom[i] = q;
  }

  cout << "------- LSH near-collision test -------" << endl;  
  cout << "- Dimension r = " << r << endl;
  cout << "- Domain: Simulation of GSA from " << q << " to " << dom[r-1] << endl;
  cout << "- Near-collision distance : ";
  if (unif) cout << "Uniform over [-" << e << ", " << e << "]" << endl;
  else cout << "Gaussian of stddev = " << e << endl;
  cout << endl;
  
  box_num.resize(r);
  lsh_lengths.resize(r);

  double p_bad = 1;  

  for (size_t i = 0; i < r; i++) {
    if (dom[i] < lsh_length) box_num[i] = 1;
    else box_num[i] = ceil(dom[i] / lsh_length);
    lsh_lengths[i] = (double) dom[i] / box_num[i];
    // p_bad /= box_num[i];
  }
  
  for (size_t i = 0; i < lsh_dim; i++) {
    if (box_num[i] > 1) {
      if (unif) p_good *= prob_admissible_uniform(e, lsh_lengths[i]);
      else p_good *= prob_admissible_gaussian(e, lsh_lengths[i]);
    }
  }

  uint64_t iteration = ceil(iter_mult / p_good);
  cout << "LSH length: " << lsh_length << ", LSH dim = " << lsh_dim 
  << ", E(p_good) (log) = " << p_good << " (" << log2(p_good) << ")"
  << ", Pr[LSH fails] = " << pow((1 - p_good), iteration) << endl;
  cout << endl;
}

matrix NearCollisionTest::gen_instance(
  size_t near_collision_num, size_t output_size) 
{
  matrix result;

  random_device rd;
  mt19937 gen(rd());
  vector<uniform_real_distribution<double>> sampler;
  for (size_t i = 0; i < r; i++) {
    sampler.push_back(uniform_real_distribution<double>(0, dom[i]));
  }
  
  normal_distribution<double> gaussian_sampler(0, e);
  uniform_real_distribution<double> unif_sampler(-e, e);

  double p_good_max = 0, p_good_min = 1.0, p_good_mean = 0, p_good_var = 0;

  for (size_t i = 0; i < near_collision_num; i++) {
    vector<double> pair1(r);
    vector<double> pair2(r);
    vector<double> error(r);
    for (size_t j = 0; j < r; j++) {
      pair1[j] = sampler[j](gen);
      error[j] = unif? unif_sampler(gen): gaussian_sampler(gen);
      while ((pair1[j] + error[j] > dom[j]) | (pair1[j] + error[j] <= 0)) {
        error[j] = unif? unif_sampler(gen): gaussian_sampler(gen);
      }
      pair2[j] = pair1[j] + error[j];
    }
    double p_good_cur = prob_admissible_fixed(error, lsh_lengths, lsh_dim);
    p_good_max = max(p_good_cur, p_good_max);
    p_good_min = min(p_good_cur, p_good_min);
    p_good_mean += p_good_cur;
    p_good_var += p_good_cur * p_good_cur;

    result.push_back(pair1);
    result.push_back(pair2);
  }

  cout << fixed;
	cout.precision(4);

  p_good_mean = p_good_mean / near_collision_num;
  p_good_var = p_good_var / near_collision_num - p_good_mean * p_good_mean;

  cout << "p_good in [" << p_good_min << ", " << p_good_max << "], E(p_good): " << p_good_mean << ", V(p_good): " << p_good_var << endl;

  while (result.size() < output_size) {
    vector<double> randvec(r);
    for (size_t j = 0; j < r; j++) {
      randvec[j] = sampler[j](gen);
    }
    result.push_back(randvec);
  }

  return result;
}

set<vector<double>> NearCollisionTest::lsh_based_search(matrix& L) 
{
  set<vector<double>> sol;
 
  random_device rd;
  mt19937 gen(rd());

  uint64_t iteration = ceil(iter_mult / p_good);
  vector<uint64_t> iter_checks{0};
  for (size_t i = 0; i < iter_mult; i++) {
    iter_checks.push_back(ceil((i+1) / p_good));
  }

  double near_collision_condition = unif? e: 6*e;  

  size_t hash_collision_num = 0;

  for (size_t cp = 0; cp < iter_mult; cp++) {
    for (size_t i = iter_checks[cp]; i < iter_checks[cp+1]; i++) {
      // Pick torus-LSH by starting points
      vector<double> lsh_starts(lsh_dim);
      for (size_t i = 0; i < lsh_dim; i++) {
        uniform_real_distribution<double> lsh_random_starts(0, lsh_lengths[i]);
        lsh_starts[i] = lsh_random_starts(gen);
      }

      // Fill LSH table
      map<vector<uint32_t>, matrix> hash_table;
      for (auto& p: L) {
        std::vector<uint32_t> address(lsh_dim);
        vector<pair<secret, vector<double>>> bin;
        for (size_t i = 0; i < lsh_dim; i++) {
          int64_t a = floor((p[i] - lsh_starts[i]) / lsh_lengths[i]);
          address[i] = a % box_num[i];
        }
        if (hash_table.count(address) > 0) {
          hash_table[address].push_back(p);
        }
        else {
          hash_table[address] = matrix{p};
        }
      }

      for (auto &iterator: hash_table) {
        auto bin = iterator.second;
        auto N = bin.size();
        for (size_t i = 0; i < N; i++) {
          for (size_t j = i+1; j < N; j++) {
            auto check = subvec(bin[i], bin[j]);            
            if (inf_norm(check) <= near_collision_condition) {
              sol.insert(check);
              hash_collision_num++;
            } 
          }
        }
      }
    }
    // cout << "   " << iter_checks[cp+1] << " LSH iterations (" << cp+1 << "-multiple): " << sol.size() << " pairs " << endl;
  }

  // cout << "# Total hash collision = " << hash_collision_num << endl;
  // cout << "(Expected 1.5 * |L|^2 * p_bad / p_good = " << 1.5 * L.size() * p_bad / p_good << endl;

  return sol;
};
