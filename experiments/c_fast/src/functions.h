#pragma once

#include "utils.h"
#include <set>

using namespace std;

///////////////////////////////////////////////////////////////////////
/* Tests for Torus LSH-based near-collision finding */
///////////////////////////////////////////////////////////////////////

class NearCollisionTest {
public:  
  NearCollisionTest(uint32_t m, uint64_t q, double e, bool unif);

  vector<vector<double>> gen_instance(
    size_t near_collision_num, size_t output_size = 0);

  set<vector<double>> lsh_based_search(
    vector<vector<double>>& L, double lsh_length, size_t lsh_dim, uint64_t iter_multiple);

  domain dom;

private:
  uint32_t m; // dimension
  double e; // stddev if gaussian, error bound if unif
  bool unif; // error is unif or gaussian
};

///////////////////////////////////////////////////////////////////////
/* Tests for the correctness of constraint setting */
///////////////////////////////////////////////////////////////////////

class ConstraintTest {
public:
  ConstraintTest(
    uint32_t n, double e, uint64_t q,
    uint32_t h, uint32_t m, bool unif);

  void gen_noisy_instance(
    matrix& M, vector<int64_t>& s, vector<double>& error);

  /* Output the largest projection dimension r of given domain D,
  such that "ambiguity * p_adm > 1" where "p_adm = Pr[x and x + e both in pi_r(D)]" */
  void set_constraint_dim(uint32_t guess_weight, double box_length);

  /* Output {(s, Ms): HW(s) = w} */
  vector<pair<secret, vector<double>>> sparse_secret_list(
    matrix& Mtrans, uint32_t d, uint32_t w);

  /* Output {s: HW(s) = w} */
  vector<secret> enumerate_secrets(uint32_t d, uint32_t w);

  /* Output {(s, Ms): HW(s) = w, ||pi_r(Ms)||_inf <= l} */
  vector<pair<secret, vector<double>>> build_list(
    matrix& M, uint32_t guess_weight, double box_length);

  set<pair<secret, secret>> check_near_collision(
    vector<pair<secret, vector<double>>>& L, secret& s, uint32_t guess_weight);

  set<pair<secret, secret>> fast_check_near_collision(
    matrix& M, secret& s, vector<double>& e, 
    vector<secret>& partial_list, uint32_t guess_weight, double box_length);
    
  // /* Torus LSH-based search method */
  // set<pair<vector<int64_t>, vector<int64_t>>> near_collision_lsh(
  //   my_list& L, double e, vector<double>& GSnorm, uint32_t h, 
  //   vector<size_t>& box_num, vector<double>& lsh_lengths, vector<double>& lsh_domain, 
  //   uint64_t iteration_multiple, bool unif);

  uint64_t proj_dim;
  uint64_t full_list_size;
  double p_adm;
  double vol_ratio;  
  domain GSnorm;

private:
  uint32_t m; // # Sample
  uint32_t n; // LWE dim
  uint32_t h; // secret weight
  double e; // stddev if gaussian, error bound if unif
  bool unif; // error is unif or gaussian
};