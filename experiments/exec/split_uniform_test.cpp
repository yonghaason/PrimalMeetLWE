#include "functions.h"
#include "cmdline.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <map>
#include <random>

using namespace std;

int main(int argc, char* argv[]) {
  
  chrono::system_clock::time_point start, end;

  cmdline::parser parser;

  parser.add<uint32_t>("m", 'm', "matrix dimension", true, 0);
  parser.add<uint32_t>("d", 'd', "secret dimension", true, 20);
  parser.add<uint32_t>("h", 'h', "secret weight", true, 0);
  parser.add<uint32_t>("w", 'w', "split weight", true, 0);
  parser.add<double>("stddev", '\0', "error stddev", false, 3);
  parser.add<double>("lastGS", '\0', "modulus param", false, 100);
  parser.add<double>("rhf", '\0', "root-hermite-factor", false, 1.05);

  parser.add<uint64_t>("repeat", '\0', "# of experiments", false, 1000);
  
  parser.parse_check(argc, argv);

  auto m = parser.get<uint32_t>("m");
  auto d = parser.get<uint32_t>("d");
  auto h = parser.get<uint32_t>("h");  
  auto w = parser.get<uint32_t>("w");
  auto stddev = parser.get<double>("stddev");
  auto lastGS = parser.get<double>("lastGS");
  auto rhf = parser.get<double>("rhf");
  auto repeat = parser.get<uint64_t>("repeat");

  vector<double> q(m);
  for (int i = 0; i < m; i++) {
    q[m-i-1] = lastGS * pow(rhf, i);
  }

  std::cout << "*********** Experiments for Heuristic 2 *********" << endl; 
  std::cout << "* Settings " << endl;
  std::cout << "- B = " << m << " X " << m << " matrix of Gram-Schmidt norm " 
            << q[0] << " -> " << q[m-1] 
            << " (GSA with rhf = " << rhf << ")" << endl;
  std::cout << "- s = ternary vector of HW(s) = " << h << endl;
  std::cout << "- M = " << m << " X " << d << " matrix, s.t. [Ms]_B = short" << endl;
  std::cout << "- L ⊂ { [Ms_1]_B: (s_1, s_2) is a " << w << "-rep pair of s }" << endl;
  std::cout << "* Goal: Compare the number of points in a ball of radius r with random center"  << endl;
  std::cout << endl;

  random_device rd;
  mt19937 gen(rd());

  matrix M; matrix B; secret s; vector<double> e;
  gen_noisy_instance(
    m, d, h, stddev, q, 
    M, B, s, e);
  auto Mtrans = transpose(M);

  vector<uint32_t> secret_nonzeros;
  vector<int> secret_parity;
  vector<uint32_t> secret_zeros;
  for (size_t i = 0; i < d; i++) 
  {
    if (s[i] != 0) secret_nonzeros.push_back(i);
    else secret_zeros.push_back(i);
  }

  set<vector<double>> Ms1_set;
  uniform_int_distribution<uint64_t> binarysampler(0, 1);

  while (Ms1_set.size() < 10000)
  {
    vector<size_t> nonzero_indices(secret_nonzeros.size());
    iota(nonzero_indices.begin(), nonzero_indices.end(), 0);
    random_shuffle(nonzero_indices.begin(), nonzero_indices.end());

    vector<size_t> zero_indices(secret_zeros.size());
    iota(zero_indices.begin(), zero_indices.end(), 0);
    random_shuffle(zero_indices.begin(), zero_indices.end());

    secret s1(d);
    for (size_t i = 0; i < h/2; i++) 
    {
      auto idx = secret_nonzeros[nonzero_indices[i]];
      s1[idx] = s[idx];
    }
    for (size_t i = 0; i < w-h/2; i++) 
    {
      auto idx = secret_zeros[zero_indices[i]];
      s1[idx] = 2*binarysampler(gen) - 1;
    }
    
    auto s2 = sub(s, s1);
    assert(hamming_weight(s2) == w);
    
    auto Ms1 = matmul(Mtrans, s1);
    Ms1 = babaiNP(Ms1, B);
    Ms1_set.insert(Ms1);
  }
  // std::cout << ".... List done" << endl;

  vector<uniform_real_distribution<double>> coord_sampler(m);
  for (size_t i = 0; i < m; i++) 
  {
    coord_sampler[i] = uniform_real_distribution<double>(-q[i]/2, q[i]/2);
  }
  uniform_real_distribution<double> r_sampler(q[0], q[m-1]/2);

  vector<double> result(101);
  vector<int> result_count(101);

  std::cout << "* Results" << endl;
  std::cout << "vol(P(B) ∩ Ball(r))/vol(P(B)) v.s. |L ∩ Ball(r)|/|L|" << endl;

  for (size_t iter = 0; iter < repeat; iter++)
  {
    auto radius = r_sampler(gen);
    vector<double> center(m);

    for (size_t i = 0; i < m; i++)
    {
      center[i] = coord_sampler[i](gen);
    }

    double vol_ratio = 1;
    for (size_t i = 0; i < m; i++)
    {
      double length = 0.0;
      length += min(radius, q[i]/2 - center[i]);
      length += min(radius, q[i]/2 + center[i]);
      vol_ratio *= length / q[i];
    }

    if (vol_ratio <= 1)
    {
      size_t count = 0;
      for (auto Ms1 : Ms1_set)
      {
        auto check = subvec(Ms1, center);
        double norm = inf_norm(check);
        if (norm <= radius) 
        {
          count++;
        }
      }

      double count_ratio = (double) count / (double) Ms1_set.size();
      result[int(vol_ratio * 100)] += count_ratio;
      result_count[int(vol_ratio * 100)] += 1;
    }
  }

  for (size_t i = 0; i < 101; i++) 
  {
    if (i % 3 == 0)
    {if (result_count[i] != 0)
      std::cout << "(" << (double) i / 100 << // << "~" << (double) (i+1) / 100 << 
      ", " << result[i] / result_count[i] << ")" << endl;
    }
  }
  
  return 0;
}