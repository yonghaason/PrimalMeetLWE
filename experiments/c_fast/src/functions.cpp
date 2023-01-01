#include "functions.h"
#include <iostream>
#include <random>
#include <algorithm>
#include <map>

using namespace std;

my_list sparse_secret_list(
  vector<vector<int64_t>>& Mtrans, uint32_t n, uint32_t w) 
{
  auto m = Mtrans[0].size();
  my_list R;
  if (n == 0) return R;
  else if (w == 0) {
    R.push_back(make_pair(vector<int64_t>(n), vector<int64_t>(m)));
    return R;
  }
  else if (n == 1 && w == 1) {
    R.push_back(make_pair(vector<int64_t>{1}, Mtrans[n-1]));
    vector<int64_t> zerovec(m);
    R.push_back(make_pair(vector<int64_t>{-1}, subvec(zerovec, Mtrans[n-1])));
    return R;
  }
  else if (n >= w) {
    auto R1 = sparse_secret_list(Mtrans, n-1, w-1);
    auto R2 = sparse_secret_list(Mtrans, n-1, w);

    for (size_t i = 0; i < R1.size(); i++) {
      auto v = R1[i].first;
      v.push_back(1);
      R.push_back(make_pair(v, addvec(R1[i].second, Mtrans[n-1])));
      v[n-1] = -1;
      R.push_back(make_pair(v, subvec(R1[i].second, Mtrans[n-1])));
    }
    for (size_t i = 0; i < R2.size(); i++) {
      auto v = R2[i].first;
      v.push_back(0);
      R.push_back(make_pair(v, R2[i].second));
    }
    return R;
  }
}

void gen_instance(
  vector<vector<int64_t>>& M, vector<int64_t>& s,
  uint32_t n, double stddev, uint64_t q, uint32_t h, uint32_t m)
{
  if (m == 0) m=n;

  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<int64_t> qrand(0, q);
  
  vector<vector<int64_t>> A;
  A.resize(m);
  for (size_t i = 0 ; i < m; i++) {
    A[i].resize(n);
    for (size_t j = 0; j < n; j++) {
      A[i][j] = qrand(gen);
    }
  }

  uniform_int_distribution<int64_t> binrand(0, 1);
  for (size_t i = 0; i < h; i++) {
    s.push_back(2*binrand(gen)-1);
  }
  for (size_t i = h; i < n; i++) {
    s.push_back(0);
  }
  shuffle(begin(s), end(s), gen);
  s.push_back(-1);

  vector<int64_t> b(m);
  for (size_t k = 0; k < n; k++) {
    if (s[k] != 0) {
      for (size_t i = 0; i < m; i++) {
        b[i] += s[k] * A[i][k];
      }
    }
  }

  normal_distribution<double> noise(0, stddev);
  for (size_t k = 0; k < m; k++) {
    b[k] += round(noise(gen));
  }
  modq(b, q);

  M.resize(m);
  for (size_t i = 0; i < m; i++) {
    M[i].resize(n + 1);
    for (size_t j = 0; j < n; j++) {
      M[i][j] = A[i][j];
    }
    M[i][n] = b[i];
  }
};

double prob_admissible_axis(double b, double ell, bool unif)
{
  if (unif) {
    assert(ell >= b);
    if (ell >= 2*b)
      return 1 - b/(2*ell);
    else // ell >= b
      return 1/4 * (1 + ell/b);
  }
  else {
    auto x = ell/(sqrt(2) * b);
    return erf(x) + (exp(-pow(x, 2)) - 1)/(x * sqrt(M_PI));
  }
}

uint32_t constraint_dim(
  uint64_t ambiguity, double b, double ell,
  vector<uint64_t>& GSnorm, double& p_adm, bool unif)
{
  uint32_t r = 1;
  p_adm = 1;
  auto m = GSnorm.size();

  while (true) {
    if (r == m) break;
    
    double cur_p_adm = 1;
    if (GSnorm[m-r] < ell) {
      cur_p_adm = prob_admissible_axis(b, GSnorm[m-r], unif);
    }
    else {
      double cur_ratio = GSnorm[m-r] / ell;
      cur_p_adm = prob_admissible_axis(b, ell, unif) / cur_ratio;
    }

    auto inv = 1.0 / (p_adm * cur_p_adm);
    if ((double) ambiguity < inv) {
      r -= 1;
      break;
    }
    else {
      p_adm *= cur_p_adm;
      r += 1;
    }
  }
  return r;
}

my_list build_list(vector<vector<int64_t>>& M,
                uint32_t w, uint32_t r,
                double l, uint64_t q)
{
  auto m = M.size();
  auto n = M[0].size();

  vector<vector<int64_t>> Mtrans;
  Mtrans.resize(n);
  for (size_t i = 0; i < n; i++) {
    Mtrans[i].resize(m);
    for (size_t j = 0; j < m; j++) {
      Mtrans[i][j] = M[j][i];
    }
  }

  // S = {(v, Mv): HW(v) = w}
  auto S = sparse_secret_list(Mtrans, n, w);

  my_list L;
  for (auto& p: S) {
    auto v = p.first;
    auto Mv = p.second;
    modq(Mv, q);

    if (inf_norm(Mv, m-r) <= l)
      L.push_back(make_pair(v, Mv));
  }
  return L;
};

vector<vector<int64_t>> near_collision_naive(
  my_list& L, double stddev, uint64_t q, uint32_t h)
{
  vector<vector<int64_t>> sol;
  auto N = L.size();
  auto n = L[0].first.size();
  for (size_t i = 0; i < N; i++) {
    for (size_t j = i+1; j < N; j++) {
      auto s_ = subvec(L[i].first, L[j].first);
      if (hamming_weight(s_) != h)
        continue;
      if (s_[n-1] == 1)
        continue;
      
      auto Ms_ = subvec(L[i].second, L[j].second);
      modq(Ms_, q);

      if (inf_norm(Ms_) < 6*stddev)
        sol.push_back(s_);
    }
  }
  return sol;
};

vector<vector<int64_t>> near_collision_lsh(
  my_list& L, double stddev, uint64_t q, uint32_t h, 
  double& lsh_length, vector<uint64_t>& domain) 
{
  vector<vector<int64_t>> sol;
  
  /**
   * 0. box_num = ceil((double) q / box_len), box_len_real = (double) q / box_num;
   * 1. Take a random center in box_len_real
   * 2. Insert each point in corresponding box
   * - How to implement?
   * map<vector<uint64_t>, list>
   *       Box Address  (s, [Ms]_q)
   * 3. For each box, check collision.
   */

  size_t box_num = ceil((double) q / lsh_length);
  lsh_length = (double) q / box_num;

  map<vector<uint32_t>, my_list> hash_table;

  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<double> lsh_random_center(0, lsh_length);

  for (auto& p: L) {
    auto Mv = p.second;
    std::vector<uint32_t> address;
    for (size_t i = 0; i < Mv.size(); i++) {
      
    }
  }
  
    
  for (size_t i = 0; i < N; i++) {
    for (size_t j = i+1; j < N; j++) {
      auto s_ = subvec(L[i].first, L[j].first);
      if (hamming_weight(s_) != h)
        continue;
      if (s_[n-1] == 1)
        continue;
      
      auto Ms_ = subvec(L[i].second, L[j].second);
      modq(Ms_, q);

      if (inf_norm(Ms_) < 6*stddev)
        sol.push_back(s_);
    }
  }
  return sol;
};