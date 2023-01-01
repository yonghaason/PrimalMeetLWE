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
  if (n == 0 | n < w) {
    return R;
  }
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
  else {
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

set<pair<vector<int64_t>, vector<int64_t>>> near_collision_naive(
  my_list& L, double stddev, uint64_t q, uint32_t h)
{
  set<pair<vector<int64_t>, vector<int64_t>>> sol;
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
        sol.insert(make_pair(L[i].first, L[j].first));
    }
  }
  return sol;
};

set<pair<vector<int64_t>, vector<int64_t>>> near_collision_lsh(
  my_list& L, double stddev, uint64_t q, uint32_t h, 
  double lsh_length, vector<double>& domain) 
{
  set<pair<vector<int64_t>, vector<int64_t>>> sol;

  random_device rd;
  mt19937 gen(rd());

  auto n = L[0].first.size();
  auto m = L[0].second.size();

  auto iteration = 1; // TODO: Should be 1/p_good .. like this

  vector<size_t> box_num(m);
  vector<size_t> lsh_length_real(m);

  for (size_t i = 0; i < m; i++) {
    if (domain[i] < lsh_length) box_num[i] = 1;
    else box_num[i] = ceil(domain[i] / lsh_length);
    lsh_length_real[i] = domain[i] / box_num[i];
  }

  for (size_t iter = 0; iter < iteration; iter++) {    
    // Pick torus-LSH by starting points
    vector<double> lsh_starts(m);
    for (size_t i = 0; i < m; i++) {
      uniform_real_distribution<double> lsh_random_starts(0, lsh_length_real[i]);
      lsh_starts[i] = lsh_random_starts(gen);
    }

    // Fill LSH table
    map<vector<uint32_t>, my_list> hash_table;
    for (auto& p: L) {
      auto Mv = p.second;
      std::vector<uint32_t> address(m);
      my_list bin;
      for (size_t i = 0; i < m; i++) {
        int64_t a = floor((Mv[i] - lsh_starts[i]) / lsh_length_real[i]);
        address[i] = a % box_num[i];
      }
      if (hash_table.count(address) > 0) {
        hash_table[address].push_back(p);
      }
      else {
        hash_table[address] = my_list{p};
      }
    }
      
    for (auto &iterator: hash_table) {
      auto bin = iterator.second;
      auto N = bin.size();
      auto n = bin[0].first.size();
      for (size_t i = 0; i < N; i++) {
        for (size_t j = i+1; j < N; j++) {
          auto s_ = subvec(bin[i].first, bin[j].first);
          if (hamming_weight(s_) != h)
            continue;
          if (s_[n-1] == 1)
            continue;
          
          auto Ms_ = subvec(bin[i].second, bin[j].second);
          modq(Ms_, q);

          if (inf_norm(Ms_) < 6*stddev)
            sol.insert(make_pair(bin[i].first, bin[j].first));
        }
      }
    }
}
  
  return sol;
};