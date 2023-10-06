#pragma once

#include "utils.h"
#include <set>

using namespace std;

void gen_uniform_matrix(
  const int m, const int d, const domain& q,
  matrix& A);

void gen_noisy_instance(
  const int m, const int d, const int h, const double stddev, const domain& q,
  matrix& M, matrix& B, vector<int64_t>& s, vector<double>& e);

list sparse_secret_list(matrix& Mtrans, uint32_t m, uint32_t d, uint32_t w);

set<secret> enumerate_secrets(uint32_t d, uint32_t w);

set<secret> NCF_lsh(
  domain dom, list& L, double ell, int h,
  double block_length, uint64_t C,
  size_t& collision_nums);