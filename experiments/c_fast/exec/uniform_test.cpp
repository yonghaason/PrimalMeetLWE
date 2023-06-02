#include "functions.h"
#include "cmdline.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <map>
#include <random>

// ./unif_test -d 20 -h 8 -l 0.5 -q 4

using namespace std;

int main(int argc, char* argv[]) {
  
  chrono::system_clock::time_point start, end;

  uint32_t d, h; 
  uint64_t q;
  double l;

  cmdline::parser parser;

  parser.add<uint32_t>("d", 'd', "Dimension (r^(i-1) in paper)", true, 20);
  parser.add<uint32_t>("h", 'h', "Secret weight (w^(i-1) in paper)", true, 0);
  parser.add<double>("l", 'l', "Constraint bound (ell^(i) in paper)", true, 0);
  parser.add<uint64_t>("q", 'q', "Modulus param that determines the domain", false, 2048);

  parser.parse_check(argc, argv);

  d = parser.get<uint32_t>("d");
  h = parser.get<uint32_t>("h");
  q = parser.get<uint64_t>("q");
  l = parser.get<double>("l");

  auto m = d;
    
  UniformTest unif_test(d, 1, q, h, m, true);
  matrix B = unif_test.B;
  auto GSnorm = unif_test.GSnorm;
  matrix M; secret s; vector<double> e;
  unif_test.gen_noisy_instance(M, s, e);
  auto Mtrans = transpose(M);

  auto Ms_list = unif_test.sparse_secret_list(Mtrans, d, h);

  cout << "Ms_list.size() = " << Ms_list.size() << endl;

  size_t count = 0;
  for (auto &pair : Ms_list)
  {
    auto Ms = pair.second;
    Ms = babaiNP(Ms, B);
    double norm = inf_norm(Ms);
    if (norm <= l) 
    {
      count++;
    }
  }

  double vol_ratio = 1;
  for (size_t i = 0; i < m; i++)
  {
    if (GSnorm[i] > 2*l) 
    {
      vol_ratio *= (2*l / GSnorm[i]);
    }
  }

  double Ms_ratio = (double) count / (double) Ms_list.size();

  cout << "Volume Ratio: " << vol_ratio << endl;
  cout << "Ms in Box Ratio: " << Ms_ratio << endl;

  return 0;
}