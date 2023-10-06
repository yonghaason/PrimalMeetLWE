#include "functions.h"
#include <iostream>
#include <random>
#include <algorithm>
#include <map>

using namespace std;

void gen_uniform_matrix(
  const int m, const int d, const domain& q,
  matrix& A)
{
  random_device rd;
  mt19937 gen(rd());
  vector<uniform_real_distribution<double>> coord_sampler(m);
  
  A.resize(m);
  for (size_t i = 0 ; i < m; i++) {
    A[i].resize(d);
    coord_sampler[i] = uniform_real_distribution<double>(-q[i]/2, q[i]/2);
    for (size_t j = 0; j < d; j++) {
      A[i][j] = coord_sampler[i](gen);
    }
  }
}

void gen_noisy_instance(
  const int m, const int d, const int h, const double stddev, const domain& q,
  matrix& M, matrix& B, vector<int64_t>& s, vector<double>& e)
{
  random_device rd;
  mt19937 gen(rd());
  vector<uniform_real_distribution<double>> coord_sampler(m);

  B.resize(m);
  for (size_t i = 0; i < m; i++) {
    B[i].resize(m);
  }
  for (size_t j = 0; j < m; j++) {
    B[j][j] = q[j];
    coord_sampler[j] = uniform_real_distribution<double>(-q[j]/2, q[j]/2);
    for (size_t i = 0; i < j; i++) {
      B[i][j] = coord_sampler[j](gen);
    }
  }
  
  matrix A;
  A.resize(m);
  for (size_t i = 0 ; i < m; i++) {
    A[i].resize(d-1);
    coord_sampler[i] = uniform_real_distribution<double>(-q[i]/2, q[i]/2);
    for (size_t j = 0; j < d-1; j++) {
      A[i][j] = coord_sampler[i](gen);
    }
  }

  s.resize(d-1);
  uniform_int_distribution<int64_t> binrand(0, 1);
  for (size_t i = 0; i < h-1; i++) {
    s[i] = 2*binrand(gen)-1;
  }
  for (size_t i = h-1; i < d-1; i++) {
    s[i] = 0;
  }
  shuffle(begin(s), end(s), gen);
  s.push_back(1);

  vector<double> b(m);
  for (size_t k = 0; k < d-1; k++) {
    if (s[k] != 0) {
      for (size_t i = 0; i < m; i++) {
        b[i] += s[k] * A[i][k];
      }
    }
  }

  e.resize(m);
  normal_distribution<double> gaussian_noise(0, stddev);
  for (size_t k = 0; k < m; k++) {
    e[k] = gaussian_noise(gen);
    b[k] += e[k];
  }
  b = babaiNP(b, B);
  
  M.resize(m);
  for (size_t i = 0; i < m; i++) {
    M[i].resize(d);
    for (size_t j = 0; j < d-1; j++) {
      M[i][j] = -A[i][j];
    }
    M[i][d-1] = b[i];
  }
};

list sparse_secret_list(
  matrix& Mtrans, uint32_t m, uint32_t d, uint32_t w) 
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
    R.push_back(make_pair(vector<int64_t>{1}, Mtrans[0]));
    vector<double> zerovec(m);
    R.push_back(make_pair(vector<int64_t>{-1}, subvec(zerovec, Mtrans[0])));
    return R;
  }
  else {
    auto R1 = sparse_secret_list(Mtrans, m, d-1, w-1);
    auto R2 = sparse_secret_list(Mtrans, m, d-1, w);

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
};

set<secret> enumerate_secrets(uint32_t d, uint32_t w) 
{
  set<secret> R;
  if (d == 0 | d < w) {
    return R;
  }
  else if (w == 0) {
    R.insert(vector<int64_t>(d));
    return R;
  }
  else if (d == 1 && w == 1) {
    R.insert(vector<int64_t>{1});
    R.insert(vector<int64_t>{-1});
    return R;
  }
  else {
    auto R1 = enumerate_secrets(d-1, w-1);
    auto R2 = enumerate_secrets(d-1, w);

    for (auto v: R1) {
      v.push_back(1);
      R.insert(v);
      v[d-1] = -1;
      R.insert(v);
    }
    for (auto v: R2) {
      v.push_back(0);
      R.insert(v);
    }
    return R;
  }
};

set<secret> NCF_lsh(
  domain dom, list& L, double ell, int h,
  double block_length, uint64_t R_lsh,
  size_t& collision_nums)
{
  auto r = dom.size();

  vector<int> n(r);
  vector<double> b(r);  
  for (size_t i = 0; i < r; i++) {
    n[i] = max(1, (int) floor(dom[i] / block_length));
    b[i] = (double) dom[i] / n[i];
  }

  random_device rd;
  mt19937 gen(rd());

  set<secret> sol;
  
  for (size_t iter = 0; iter < R_lsh; iter++)
  {
    // Pick torus-LSH by starting points
    vector<double> lsh_starts(r);
    for (size_t i = 0; i < r; i++) {
      uniform_real_distribution<double> lsh_random_starts(0, b[i]);
      lsh_starts[i] = lsh_random_starts(gen);
    }

    // Fill LSH table
    map<vector<int32_t>, list> hash_table;
    for (auto s_Ms: L) {
      auto Ms = s_Ms.second;
      std::vector<int32_t> address(r);
      for (size_t i = 0; i < r; i++) {
        int32_t a = floor((Ms[i] + dom[i] - lsh_starts[i]) / b[i]);
        address[i] = a % n[i];
      }
      if (hash_table.count(address) != 0) {
        hash_table[address].push_back(s_Ms);
      }
      else {
        list s_Ms_copy{s_Ms};
        hash_table[address] = s_Ms_copy;
      }
    }

    for (auto &iterator: hash_table) 
    {
      auto bin = iterator.second;
      auto N = bin.size();
      for (size_t i = 0; i < N; i++) 
      {
        for (size_t j = i+1; j < N; j++) {
          collision_nums++;
          auto norm_check = subvec(bin[i].second, bin[j].second);
          if (inf_norm(norm_check) <= ell) 
          {
            auto weight_check = sub(bin[i].first, bin[j].first);
            if (weight_ternary_check(weight_check, h))
            {
              sol.insert(weight_check);
              // Symmetry
              secret zero(weight_check.size());
              weight_check = sub(zero, weight_check);
              sol.insert(weight_check);
            }
          }
        }
      }
    }
  }
  return sol;
};