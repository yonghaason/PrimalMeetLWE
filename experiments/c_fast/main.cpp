#include "src/functions.h"
#include <iostream>
#include <iomanip>
#include <map>

using namespace std;

int main(int argc, char* argv[]) {
  uint32_t n, m, h; 
  uint64_t q;
  double stddev; 
  
  n = 20; m = 20; stddev = 3.2; q = 1024; h = 8;

  double ell = 4*stddev;
  uint32_t w = 5;
  vector<uint64_t> GSnorm(m);
  for (size_t i = 0; i < m; i++) GSnorm[i] = q;

  auto R = ambiguity(n, h, w);
  double p_adm = 1.0;

  auto r = constraint_dim(R, stddev, ell, GSnorm, p_adm, false);

  cout << "------- Params -------"  << endl;
  cout << "- # of reps (Ambiguity): " << R << endl;
  cout << "- Projection dim: " << r << endl;
  cout << "- p_admissible (log): " << p_adm << " (" << log2(p_adm) << ")" << endl;
  cout << "- Expected reps (R * p_adm): " << R * p_adm << endl;
  
  // Actual List Construction & Collision Finding
  cout << "------- Actual Test -------" << endl;
  auto iteration = 50;
  double avg_sol = 0.0;
  double avg_list = 0.0;
  map<size_t, size_t> stats_sol;
  for (size_t i = 0; i < iteration; i++) {
    vector<vector<int64_t>> M;
    vector<int64_t> s;
    gen_instance(M, s, n-1, stddev, q, h-1, m);
    auto L = build_list(M, w, r, ell, q);
    auto sol = near_collision_naive(L, stddev, q, h);
    avg_sol = (avg_sol * i + (double) sol.size()) / (i+1);
    avg_list = (avg_list * i + (double) L.size()) / (i+1);
    cout << "- |L_0| = " << sol.size() 
         << " / |L_1| = " << L.size() 
         << std::setw(5) << " / E(|L_0|) = " << avg_sol
         << std::setw(5) << " / E(|L_1|) = " << avg_list << endl;
    auto filled = stats_sol.count(sol.size());
    if (filled == 1) {stats_sol[sol.size()]++;}
    else {stats_sol.insert({sol.size(), 1});}
  }

  cout << "\nStats" << endl;
  for (auto stat: stats_sol) {
    std::cout << "- \'|L_0| = " << stat.first << "\' occurs " << stat.second << " times" << endl;
  }
  return 0;
}