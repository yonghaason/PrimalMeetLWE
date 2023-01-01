#include "utils.h"
#include <iostream>

uint64_t inf_norm(vector<int64_t>& v, size_t start) {
  int64_t norm = 0;
  for (size_t i = start; i < v.size(); i++) {
    norm = max(norm, abs(v[i]));
  }
  return norm;
}

uint64_t binom(uint32_t n, uint32_t k) {
  if (k == n || k == 0) return 1;
  return binom(n - 1, k - 1) + binom(n - 1, k);
};

uint64_t ambiguity(uint32_t n, uint32_t h, uint32_t w) {  
  auto eps = w - (h+1) / 2;
  return binom(h, w - eps) * binom(n - h, eps) * (1ull << eps);
};

vector<int64_t> addvec(vector<int64_t>& a, vector<int64_t>& b) {
  assert(a.size() == b.size());
  vector<int64_t> c(a.size());
  for (size_t i = 0; i < a.size(); i++) {
    c[i] = a[i] + b[i];
  }
  return c;
}

vector<int64_t> subvec(vector<int64_t>& a, vector<int64_t>& b) {
  assert(a.size() == b.size());
  vector<int64_t> c(a.size());
  for (size_t i = 0; i < a.size(); i++) {
    c[i] = a[i] - b[i];
  }
  return c;
}

void modq(vector<int64_t>& v, uint64_t q, bool balanced) {
  for (size_t i = 0; i < v.size(); i++)
    v[i] %= q;
  if (balanced) {
    for (size_t i = 0; i < v.size(); i++) {
      if (v[i] > q/2) v[i] -= q;
    }
  }    
}

uint32_t hamming_weight(vector<int64_t>& v) {
  uint32_t weight = 0;
  for (size_t i = 0; i < v.size(); i++)
    if (v[i] != 0) weight++;
  return weight;
}

void print_list(my_list L) {
  for (size_t i = 0; i < L.size(); i++) {
    auto v = L[i].first;
    cout << "[";
    for (size_t k = 0; k < v.size(); k++){
      cout << v[k] << ", ";
    } 
    cout << "]" << endl;
  }
  cout << endl;
}