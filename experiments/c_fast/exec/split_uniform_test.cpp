#include "functions.h"
#include "cmdline.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <map>
#include <random>

// ./split_uniform_test -m 10 -d 20 -h 10 -w 6 -q 4

using namespace std;

int main(int argc, char* argv[]) {
  
  chrono::system_clock::time_point start, end;

  cmdline::parser parser;

  parser.add<uint32_t>("m", 'm', "matrix dimension", true, 0);
  parser.add<uint32_t>("d", 'd', "secret dimension", true, 20);
  parser.add<uint32_t>("h", 'h', "secret weight", true, 0);
  parser.add<uint32_t>("w", 'w', "split weight", true, 0);
  parser.add<double>("stddev", '\0', "error stddev", false, 1);
  parser.add<double>("q", 'q', "modulus param", false, 2048);
  parser.add<double>("rhf", '\0', "root-hermite-factor", false, 1.05);
  parser.add<double>("step", '\0', "step of 1-norm box radius", false, 0.1);
  
  parser.parse_check(argc, argv);

  auto m = parser.get<uint32_t>("m");
  auto d = parser.get<uint32_t>("d");
  auto h = parser.get<uint32_t>("h");  
  auto w = parser.get<uint32_t>("w");
  auto stddev = parser.get<double>("stddev");
  auto q0 = parser.get<double>("q");
  auto rhf = parser.get<double>("rhf");
  auto step = parser.get<double>("step");

  vector<double> q(m);
  for (int i = 0; i < m; i++) {
    q[m-i-1] = q0 * pow(rhf, i);
  }

  cout << "*********** Experiments for Heuristic 2 *********" << endl; 
  cout << "* Settings " << endl;
  cout << "- B = " << m << " X " << m << " matrix of Gram-Schmidt norm " 
            << q[0] << " -> " << q[m-1] 
            << " (GSA with rhf = " << rhf << ")" << endl;
  cout << "- s = ternary vector of HW(s) = " << h << endl;
  cout << "- M = " << m << " X " << d << " matrix, s.t. [Ms]_B = short" << endl;
  cout << "- L = { [Ms_1]_B: (s_1, s_2) is a " << w << "-rep pair of s }" << endl;
  cout << "* Goal: Compare the number of points in a ball of radius r"  << endl;
  cout << endl;

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

  cout << "* Results" << endl;
  cout << "r: |L ∩ Ball(r)|/|L| v.s. vol(Ball(r))/vol(P(B))" << endl;

  double radius = 0.0;
  while (true) 
  {
    radius += step;
    if (radius >= q[0]/2) break;

    double vol_ratio = 1;
    for (size_t i = 0; i < m; i++)
    {
      if (q[i] > 2*radius) 
      {
        vol_ratio *= (2*radius / q[i]);
      }
    }

    if (vol_ratio >= 0.01)
    {
      size_t count = 0;
      for (auto Ms1 : Ms1_set)
      {
        double norm = inf_norm(Ms1);
        if (norm <= radius) 
        {
          count++;
        }
      }

      double count_ratio = (double) count / (double) Ms1_set.size();
      cout << radius << ": " << count_ratio << " v.s. " << vol_ratio << endl;
    }
  }

  return 0;
}