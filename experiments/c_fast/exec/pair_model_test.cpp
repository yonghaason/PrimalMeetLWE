#include "functions.h"
#include "cmdline.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <map>
#include <random>

// ./pair_model_test -d 20 -h 8 -w 5 -l 0.5 -q 4

using namespace std;

int main(int argc, char* argv[]) {
  
  chrono::system_clock::time_point start, end;

  cmdline::parser parser;

  parser.add<uint32_t>("m", 'm', "matrix dimension", true, 0);
  parser.add<uint32_t>("d", 'd', "secret dimension", true, 20);
  parser.add<uint32_t>("h", 'h', "secret weight", true, 0);
  parser.add<uint32_t>("w", 'w', "split weight", true, 0);
  parser.add<double>("l", 'l', "box length [-l, l]", true, 0);
  parser.add<double>("stddev", '\0', "error stddev", false, 1);
  parser.add<double>("q", 'q', "modulus param", false, 2048);
  parser.add<double>("rhf", '\0', "root-hermite-factor", false, 1.05);

  parser.parse_check(argc, argv);

  auto m = parser.get<uint32_t>("m");
  auto d = parser.get<uint32_t>("d");
  auto h = parser.get<uint32_t>("h");  
  auto w = parser.get<uint32_t>("w");  
  auto l = parser.get<double>("l");
  auto stddev = parser.get<double>("stddev");
  auto q0 = parser.get<double>("q");
  auto rhf = parser.get<double>("rhf");

  vector<double> q(m);
  for (int i = 0; i < m; i++) {
    q[m-i-1] = q0 * pow(rhf, i);
  }

  std::cout << "Coordinate length (q) = " << q[0] << " -> " << q[m-1] << endl;

  matrix M; matrix B; secret s; vector<double> e;
  gen_noisy_instance(
    m, d, h, stddev, q, 
    M, B, s, e);
  auto Mtrans = transpose(M);

  auto candidates = enumerate_secrets(d, w);
  set<vector<double>> Ms1_set;
  
  for (auto s1 : candidates)
  {
    auto s2 = sub(s, s1);
    if (hamming_weight(s2) == w)
    {
      auto Ms1 = matmul(Mtrans, s1);
      Ms1 = babaiNP(Ms1, B);
      Ms1_set.insert(Ms1);
    }
  }

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
    if (q[i] > 2*l) 
    {
      vol_ratio *= (2*l / q[i]);
    }
  }

  double count_ratio = (double) count / (double) Ms1_set.size();

  cout << "Volume Ratio: " << vol_ratio << endl;
  cout << "Ms in Box Ratio: " << count_ratio << endl;

  return 0;
}