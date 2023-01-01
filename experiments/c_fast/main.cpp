#include "src/functions.h"
#include <iostream>
#include <iomanip>
#include <chrono>
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


  chrono::system_clock::time_point start, end;

  cout << "------- Actual Test -------" << endl; 
  vector<vector<int64_t>> M;
  vector<int64_t> s;
  gen_instance(M, s, n-1, stddev, q, h-1, m);
  auto L = build_list(M, w, r, ell, q);
  vector<double> lsh_domain(m);
  for (size_t i = 0; i < m - r; i++) {
    lsh_domain[i] = GSnorm[i];
  }
  for (size_t i = m - r; i < m; i++) {
    lsh_domain[i] = ell;
  }
  start = chrono::system_clock::now();
  auto sol = near_collision_lsh(L, stddev, q, h, 6*stddev, lsh_domain);
  end = chrono::system_clock::now();
  cout << "LSH takes " << chrono::duration_cast<chrono::milliseconds>(end-start).count() << " ms"<< endl;

  start = chrono::system_clock::now();
  auto sol2 = near_collision_naive(L, stddev, q, h);
  end = chrono::system_clock::now();
  cout << "Naive takes " << chrono::duration_cast<chrono::milliseconds>(end-start).count() << " ms"<< endl;

  cout << "- |L_0|_naive = " << sol2.size() 
        << " / |L_0|_LSH = " << sol.size()
        << " / |L_1| = " << L.size() << endl;
  
  // Actual List Construction & Collision Finding
  // cout << "------- Actual Test -------" << endl;
  // auto iteration = 50;
  // double avg_sol = 0.0;
  // double avg_list = 0.0;
  // map<size_t, size_t> stats_sol;
  // for (size_t i = 0; i < iteration; i++) {
  //   vector<vector<int64_t>> M;
  //   vector<int64_t> s;
  //   gen_instance(M, s, n-1, stddev, q, h-1, m);
  //   auto L = build_list(M, w, r, ell, q);
  //   auto sol = near_collision_naive(L, stddev, q, h);
  //   avg_sol = (avg_sol * i + (double) sol.size()) / (i+1);
  //   avg_list = (avg_list * i + (double) L.size()) / (i+1);
  //   cout << "- |L_0| = " << sol.size() 
  //        << " / |L_1| = " << L.size() 
  //        << setw(5) << " / E(|L_0|) = " << avg_sol
  //        << setw(5) << " / E(|L_1|) = " << avg_list << endl;
  //   auto filled = stats_sol.count(sol.size());
  //   if (filled == 1) {stats_sol[sol.size()]++;}
  //   else {stats_sol.insert({sol.size(), 1});}
  // }

  // cout << "\nStats" << endl;
  // for (auto stat: stats_sol) {
  //   cout << "- \'|L_0| = " << stat.first << "\' occurs " << stat.second << " times" << endl;
  // }

  return 0;
}