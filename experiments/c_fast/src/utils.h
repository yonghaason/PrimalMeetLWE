#pragma once

#include <stdint.h>
#include <cmath>
#include <vector>
#include <cassert>

#define _USE_MATH_DEFINES

using namespace std;

using my_list = vector<pair<vector<int64_t>, vector<int64_t>>>;

void modq(vector<int64_t>& v, uint64_t q, bool balanced = true);
uint64_t inf_norm(vector<int64_t>& v, size_t start = 0);
uint32_t hamming_weight(vector<int64_t>& v);

uint64_t binom(uint32_t n, uint32_t k);
uint64_t ambiguity(uint32_t n, uint32_t h, uint32_t w);

vector<int64_t> addvec(vector<int64_t>& a, vector<int64_t>& b);
vector<int64_t> subvec(vector<int64_t>& a, vector<int64_t>& b);

void print_list(my_list L);
