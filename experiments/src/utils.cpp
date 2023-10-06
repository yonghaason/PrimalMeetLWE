#include "utils.h"
#include <iostream>

double prob_admissible_gaussian(double stddev, double ell)
{
  auto x = ell/(sqrt(2)*stddev);
  auto p = erf(x) + (exp(-pow(x, 2)) - 1)/(x * sqrt(M_PI));
  return p;
}

double prob_admissible_uniform(double b, double ell)
{ 
  if (b <= ell) return 1.0 - b/(2*ell);
  else return ell/(2*b);
}

double prob_admissible_fixed(vector<double> error, double ell, size_t proj_dim)
{
  double prob = 1.0;
  if (proj_dim == 0) proj_dim = error.size();
  for (size_t i = error.size() - proj_dim; i < error.size(); i++)
  {
    prob *= (1.0 - abs(error[i]) / ell);
  }
  return prob;
}

double prob_admissible_fixed(vector<double> error, vector<double> lengths, size_t lsh_dim)
{
  double prob = 1.0;
  for (size_t i = 0; i < lsh_dim; i++)
  {
    prob *= (1.0 - abs(error[i]) / lengths[i]);
  }
  return prob;
}

bool weight_ternary_check(secret& s, uint32_t weight)
{
  if (hamming_weight(s) != weight) return false;
  for (size_t i = 0; i < s.size(); i++) {
    if (abs(s[i]) > 1) return false;
  }
  return true;
}

double inf_norm(domain& v, size_t start) {
  double norm = 0;
  for (size_t i = start; i < v.size(); i++) {
    norm = max(norm, abs(v[i]));
  }
  return norm;
}

matrix transpose(matrix& M) {
  matrix A;
  A.resize(M[0].size());
  for (size_t i = 0; i < A.size(); i++) {
    A[i].resize(M.size());
    for (size_t j = 0; j < A[i].size(); j++) {
      A[i][j] = M[j][i];
    }
  }
  return A;
}

// uint64_t binom(uint32_t n, uint32_t k) {
//   if (n == 0) return 1;
//   if (k == n || k == 0) return 1;
//   return binom(n - 1, k - 1) + binom(n - 1, k);
// };

uint64_t binom(uint32_t n, uint32_t k) {
  if (k == 0) return 1;
  std::vector<uint64_t> aSolutions(k);
  aSolutions[0] = n - k + 1;

  for (int i = 1; i < k; ++i) {
    aSolutions[i] = aSolutions[i - 1] * (n - k + 1 + i) / (i + 1);
  }

  return aSolutions[k - 1];
}

uint64_t ambiguity(uint32_t n, uint32_t h, uint32_t w) {
  assert(h % 2 == 0);
  auto eps = w - h/2;
  return binom(h, w - eps) * binom(n - h, eps) * (1ull << eps);
};

vector<double> addvec(vector<double>& a, vector<double>& b) {
  assert(a.size() == b.size());
  vector<double> c(a.size());
  for (size_t i = 0; i < a.size(); i++) {
    c[i] = a[i] + b[i];
  }
  return c;
}

vector<double> subvec(vector<double>& a, vector<double>& b) {
  assert(a.size() == b.size());
  vector<double> c(a.size());
  for (size_t i = 0; i < a.size(); i++) {
    c[i] = a[i] - b[i];
  }
  return c;
}

secret add(secret& a, secret& b) {
  assert(a.size() == b.size());
  secret c(a.size());
  for (size_t i = 0; i < a.size(); i++) {
    c[i] = a[i] + b[i];
  }
  return c;
}

secret sub(secret& a, secret& b) {
  assert(a.size() == b.size());
  secret c(a.size());
  for (size_t i = 0; i < a.size(); i++) {
    c[i] = a[i] - b[i];
  }
  return c;
}

void fmodvec(vector<double>& v, domain& GSnorm, bool balanced) {
  for (size_t i = 0; i < v.size(); i++) {
    v[i] = fmod(v[i], GSnorm[i]);
  }
  if (balanced) {
    for (size_t i = 0; i < v.size(); i++) {
      if (v[i] > GSnorm[i]/2) {v[i] -= GSnorm[i];}
      else if (v[i] < -GSnorm[i]/2) {v[i] += GSnorm[i];}
    }
  }
}

vector<double> babaiNP(const vector<double>& v, matrix B, int r) {
  // NOT optimized yet
  assert(v.size() == B.size());
  int m = B.size();
  if (r == -1) r = m;
  vector<double> v_reduced(v);
  for (int i = m-1; i >= m-r; i--) 
  {
    int coeff = round(v_reduced[i] / B[i][i]);
    for (int j = i; j >= m-r; j--)
    {
      v_reduced[j] -= coeff * B[j][i];
    }
  }
  vector<double> result(v_reduced.begin() + m - r, v_reduced.end());
  return result;
}

uint32_t hamming_weight(secret& v) {
  uint32_t weight = 0;
  for (size_t i = 0; i < v.size(); i++)
    if (v[i] != 0) weight++;
  return weight;
}

vector<double> matmul(matrix& Mtrans, secret& s) {
  vector<double> res(Mtrans[0].size());
  for (size_t i = 0; i < s.size(); i++) {
    if (s[i] == 1) res = addvec(res, Mtrans[i]);
    else if (s[i] == -1) res = subvec(res, Mtrans[i]);
  }
  return res;
}

void print(secret& v) {
  cout << "[";
  for (size_t i = 0; i < v.size()-1; i++) {
    cout << v[i] << ", ";
  }
  cout << v[v.size()-1] << "]" << endl;
}

void print(vector<double>& v) {
  cout << "[";
  for (size_t i = 0; i < v.size()-1; i++) {
    cout << v[i] << ", ";
  }
  cout << v[v.size()-1] << "]" << endl;
}
