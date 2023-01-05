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

uint64_t ambiguity(uint32_t n, uint32_t h, uint32_t w);

// Prob[x + e in [-b, b]], 
// where x <-[-b, b] and e <- gaussian of stddev, or unif[-e, e]
double prob_admissible_gaussian(double stddev, double b);
double prob_admissible_uniform(double e, double b);

bool weight_ternary_check(secret& s, uint32_t weight);

// Basic Maths
uint64_t binom(uint32_t n, uint32_t k);
vector<double> addvec(vector<double>& a, vector<double>& b);
vector<double> subvec(vector<double>& a, vector<double>& b);
secret add(secret& a, secret& b);
secret sub(secret& a, secret& b);
void fmodvec(vector<double>& v, domain& GSnorm, bool balanced = true);
double inf_norm(domain& v, size_t start = 0);
matrix transpose(matrix& M);
vector<double> matmul(matrix& Mtrans, secret& s);
uint32_t hamming_weight(secret& v);

void print(secret& v);
void print(vector<double>& v);
// void print_list(my_list L);
