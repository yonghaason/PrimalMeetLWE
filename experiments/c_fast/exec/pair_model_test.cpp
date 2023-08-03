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

  uint32_t d, h, w;
  uint64_t q;
  double l;

  cmdline::parser parser;

  parser.add<uint32_t>("d", 'd', "Dimension (r^(i-1) in paper)", true, 20);
  parser.add<uint32_t>("h", 'h', "Secret weight (w^(i-1) in paper)", true, 0);
  parser.add<uint32_t>("w", 'w', "Guess weight (w^i in paper)", true, 0);
  parser.add<double>("l", 'l', "Constraint bound (ell^(i) in paper)", true, 0);
  parser.add<uint64_t>("q", 'q', "Modulus param that determines the domain", false, 2048);

  parser.parse_check(argc, argv);

  d = parser.get<uint32_t>("d");
  h = parser.get<uint32_t>("h");
  w = parser.get<uint32_t>("w");
  q = parser.get<uint64_t>("q");
  l = parser.get<double>("l");

  auto m = d;
  
  UniformTest unif_test(d, 1, q, h, m, true);
  matrix B = unif_test.B;
  auto GSnorm = unif_test.GSnorm;
  matrix M; secret s; vector<double> e;
  unif_test.gen_noisy_instance(M, s, e);
  auto Mtrans = transpose(M);

  auto candidates = unif_test.enumerate_secrets(d, w);
  set<vector<double>> Ms1_set;
  
  for (auto &s1 : candidates)
  {
    auto s2 = sub(s, s1);
    if (hamming_weight(s2) == w)
    {
      auto Ms1 = matmul(Mtrans, s1);
      Ms1 = babaiNP(Ms1, B);
      Ms1_set.insert(Ms1);
    }
  }

  // auto Ms1_set = unif_test.sparse_secret_list(Mtrans, d, h);

  cout << "Ms1_set.size() = " << Ms1_set.size() << endl;

  size_t count = 0;
  for (auto Ms1 : Ms1_set)
  {
    double norm = inf_norm(Ms1);
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

  double count_ratio = (double) count / (double) Ms1_set.size();

  cout << "Volume Ratio: " << vol_ratio << endl;
  cout << "Ms in Box Ratio: " << count_ratio << endl;

  return 0;
}