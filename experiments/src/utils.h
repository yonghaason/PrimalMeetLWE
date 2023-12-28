#pragma once

#include <stdint.h>
#include <cmath>
#include <vector>
#include <cassert>

#define _USE_MATH_DEFINES

using namespace std;

using secret = vector<int64_t>;
using domain = vector<double>;
using matrix = vector<vector<double>>;
using list = vector<pair<secret, vector<double>>>;

uint64_t ambiguity(uint32_t n, uint32_t h, uint32_t w);

// Prob[x + e in [0, ell]], 
// where x <-[0, ell] and e <- gaussian of stddev, or unif[-b, b]
double prob_admissible_gaussian(double stddev, double ell);
double prob_admissible_uniform(double b, double ell);
double prob_admissible_fixed(vector<double> error, double ell, size_t proj_dim = 0);
double prob_admissible_fixed(vector<double> error, vector<double> ell, size_t lsh_dim);

bool weight_ternary_check(secret& s, uint32_t weight);

// Basic Maths
uint64_t binom(uint32_t n, uint32_t k);
vector<double> addvec(vector<double>& a, vector<double>& b);
vector<double> subvec(vector<double>& a, vector<double>& b);
secret add(secret& a, secret& b);
secret sub(secret& a, secret& b);
void fmodvec(vector<double>& v, domain& GSnorm, bool balanced = true);
vector<double> babaiNP(const vector<double>& v, matrix B, int r = -1);
double inf_norm(domain& v, size_t start = 0);
matrix transpose(matrix& M);
vector<double> matmul(matrix& Mtrans, secret& s);
uint32_t hamming_weight(secret& v);

void print(secret& v);
void print(vector<double>& v);