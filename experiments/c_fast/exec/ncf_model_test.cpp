#include "functions.h"
#include "cmdline.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <map>
#include <random>

// ./ncf_model_test -d 20 -h 8 -l 0.5 -q 4

using namespace std;

int main(int argc, char* argv[]) {
  
  chrono::system_clock::time_point start, end;

  cmdline::parser parser;

  parser.add<uint32_t>("m", 'm', "matrix dimension", true, 0);
  parser.add<uint32_t>("d", 'd', "secret dimension", true, 20);
  parser.add<uint32_t>("h", 'h', "secret weight", true, 0);
  parser.add<double>("l", 'l', "box length [-l, l]", true, 0);
  parser.add<double>("stddev", '\0', "error stddev", false, 1);
  parser.add<double>("q", 'q', "modulus param", false, 2048);
  parser.add<double>("rhf", '\0', "root-hermite-factor", false, 1.05);

  parser.parse_check(argc, argv);

  auto m = parser.get<uint32_t>("m");
  auto d = parser.get<uint32_t>("d");
  auto h = parser.get<uint32_t>("h");  
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

  auto Ms_list = sparse_secret_list(Mtrans, m, d, h);

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
    if (q[i] > 2*l) 
    {
      vol_ratio *= (2*l / q[i]);
    }
  }

  double Ms_ratio = (double) count / (double) Ms_list.size();

  cout << "Volume Ratio: " << vol_ratio << endl;
  cout << "Ms in Box Ratio: " << Ms_ratio << endl;

  return 0;
}