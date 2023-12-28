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
  parser.add<uint32_t>("w", 'w', "secret weight", true, 0);
  parser.add<uint32_t>("r", 'r', "projection dimension", true, 0);
  parser.add<double>("b", 'b', "box length [-b, b]", true, 0);
  parser.add<double>("stddev", '\0', "error stddev", false, 3);
  parser.add<double>("lastGS", '\0', "last Gram-Schmidt norm", false, 100);
  parser.add<double>("rhf", '\0', "root-hermite-factor", false, 1.05);

  parser.add<uint64_t>("repeat", '\0', "# of experiments", false, 1000);
  
  parser.parse_check(argc, argv);

  auto m = parser.get<uint32_t>("m");
  auto d = parser.get<uint32_t>("d");
  auto w = parser.get<uint32_t>("w");
  auto r = parser.get<uint32_t>("r");
  auto b = parser.get<double>("b");
  auto stddev = parser.get<double>("stddev");
  auto lastGS = parser.get<double>("lastGS");
  auto rhf = parser.get<double>("rhf");
  auto repeat = parser.get<uint64_t>("repeat");

  vector<double> q(m);
  for (int i = 0; i < m; i++) {
    q[m-i-1] = lastGS * pow(rhf, i);
  }

  domain D = q;
  for (size_t i = m-r; i++; i < m)
  {
    D[i] = 2*b;
  }

  cout << "*********** Experiments for Heuristic 4 *********" << endl; 
  cout << "* Settings " << endl;
  cout << "- B = " << m << " X " << m << " matrix of Gram-Schmidt norm " 
            << q[0] << " -> " << q[m-1] 
            << " (GSA with rhf = " << rhf << ")" << endl;
  cout << "- D = " << m << "-dimensional cube of i-th coord length " 
       << "q[i] if i < " << m-r << " and " << 2*b << " if i >= " << m-r << endl; 
  cout << "- M = " << m << " X " << d << " matrix, s.t. [Ms]_B = short" << endl;
  cout << "- L ⊂ { [Ms]_B: HW(s) = " << w << " & [Ms]_B ∈ D }" << endl;
  std::cout << "* Goal: Compare the number of points in a ball of radius r with random center"  << endl;
  std::cout << endl;
    
  random_device rd;
  mt19937 gen(rd());

  matrix M; matrix B; secret s; vector<double> e;
  gen_noisy_instance(
    m, d, d-1, stddev, q, 
    M, B, s, e);
  auto Mtrans = transpose(M);

  set<vector<double>> Ms_list;
  uniform_int_distribution<uint64_t> binarysampler(0, 1);

  while (Ms_list.size() < 10000)
  {
    secret stmp(w);
    for (size_t i = 0; i < w; i++) {stmp[i] = 2*binarysampler(gen) - 1;}
    stmp.resize(d);
    random_shuffle(stmp.begin(), stmp.end());
    auto Ms = matmul(Mtrans, stmp);
    
    Ms = babaiNP(Ms, B);
    Ms_list.insert(Ms);
    // double norm = inf_norm(Ms, Ms.size() - r);
    // if (norm <= b) Ms_list.insert(Ms);
  }

  vector<uniform_real_distribution<double>> coord_sampler(m);
  for (size_t i = 0; i < m; i++) 
  {
    coord_sampler[i] = uniform_real_distribution<double>(-D[i]/2, D[i]/2);
  }
  uniform_real_distribution<double> r_sampler(D[0], D[m-1]/2);

  vector<double> result(101);
  vector<int> result_count(101);

  std::cout << "* Results" << endl;
  std::cout << "vol(P(D) ∩ Ball(r))/vol(P(D)), |L ∩ Ball(r)|/|L|" << endl;

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
      length += min(radius, D[i]/2 - center[i]);
      length += min(radius, D[i]/2 + center[i]);
      vol_ratio *= length / D[i];
    }

    if (vol_ratio <= 1)
    {
      size_t count = 0;
      for (auto Ms : Ms_list)
      {
        auto check = subvec(Ms, center);
        double norm = inf_norm(check);
        if (norm <= radius) 
        {
          count++;
        }
      }

      double count_ratio = (double) count / (double) Ms_list.size();
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