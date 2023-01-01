#pragma once

#include "utils.h"

using namespace std;

double prob_admissible_axis(
  double b, double ell, bool unif = true);

uint32_t constraint_dim(
  uint64_t ambiguity, double b, double ell,
  vector<uint64_t>& GSnorm, double& p_adm,
  bool unif = true);

my_list sparse_secret_list(
  vector<vector<int64_t>>& Mtrans, uint32_t n, uint32_t w);

void gen_instance(
  vector<vector<int64_t>>& M, vector<int64_t>& s,
  uint32_t n, double stddev, uint64_t q, uint32_t h, uint32_t m = 0);

my_list build_list(
  vector<vector<int64_t>>& M, uint32_t w, uint32_t r, double l, uint64_t q);

vector<vector<int64_t>> near_collision_naive(
  my_list& L, double stddev, uint64_t q, uint32_t h);
vector<vector<int64_t>> near_collision_lsh(
  my_list& L, double stddev, uint64_t q, uint32_t h, 
  double& lsh_length, vector<uint64_t>& domain);