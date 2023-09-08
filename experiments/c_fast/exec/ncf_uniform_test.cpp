#include "functions.h"
#include "cmdline.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <map>
#include <random>

using namespace std;

// ./ncf_uniform_test -m 15 -d 10 -w 7 -r 5 -b 4 -q 8

list experiment(
  matrix& M, uint32_t m, uint32_t d, uint32_t h, matrix& B, uint32_t b, uint32_t r);

int main(int argc, char* argv[]) {
  
  chrono::system_clock::time_point start, end;

  cmdline::parser parser;

  parser.add<uint32_t>("m", 'm', "matrix dimension", true, 0);
  parser.add<uint32_t>("d", 'd', "secret dimension", true, 20);
  parser.add<uint32_t>("w", 'w', "secret weight", true, 0);
  parser.add<uint32_t>("r", 'r', "projection dimension", true, 0);
  parser.add<double>("b", 'b', "box length [-b, b]", true, 0);
  parser.add<double>("stddev", '\0', "error stddev", false, 1);
  parser.add<double>("q", 'q', "modulus param", false, 2048);
  parser.add<double>("rhf", '\0', "root-hermite-factor", false, 1.05);
  parser.add<double>("step", '\0', "step of 1-norm box radius", false, 0.1);
  
  parser.parse_check(argc, argv);

  auto m = parser.get<uint32_t>("m");
  auto d = parser.get<uint32_t>("d");
  auto w = parser.get<uint32_t>("w");
  auto r = parser.get<uint32_t>("r");
  auto b = parser.get<double>("b");
  auto stddev = parser.get<double>("stddev");
  auto q0 = parser.get<double>("q");
  auto rhf = parser.get<double>("rhf");
  auto step = parser.get<double>("step");


  vector<double> q(m);
  for (int i = 0; i < m; i++) {
    q[m-i-1] = q0 * pow(rhf, i);
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
  cout << "- L = { [Ms]_B: HW(s) = " << w << " & [Ms]_B ∈ D }" << endl;
  cout << "* Goal: Compare the number of points in a ball of radius r"  << endl;
  cout << endl;
  
  size_t full_size = binom(d, w) * (1ull << w); 
  
  matrix M; matrix B; secret s; vector<double> e;
  gen_noisy_instance(
    m, d, d-1, stddev, q, 
    M, B, s, e);

  auto Ms_list = experiment(
    M, m, d, w, B, b, r);

  cout << "* Results" << endl;
  cout << "r: |L ∩ Ball(r)|/|L| v.s. vol(Ball(r))/vol(D)" << endl;
  
  double radius = 0.0;
  while (true) 
  {
    radius += step;
    if (radius >= q[0]/2) break;
    
    double vol_ratio = 1;
    for (size_t i = 0; i < m; i++)
    {
      if (D[i] > 2*radius)
      {
        vol_ratio *= (2*radius / D[i]);
      }
    }

    if (vol_ratio >= 0.01)
    {
      size_t count = 0;
      for (auto Ms : Ms_list)
      {
        double norm = inf_norm(Ms.second);
        if (norm <= radius) 
        {
          count++;
        }
      }

      double count_ratio = (double) count / (double) Ms_list.size();
      cout << radius << ": " << count_ratio << " v.s. " << vol_ratio << endl;
    }
  }

  return 0;
}

list experiment(
  matrix& M, uint32_t m, uint32_t d, uint32_t h, matrix& B, uint32_t b, uint32_t r)
{
  auto Mtrans = transpose(M);

  auto Ms_full_list = sparse_secret_list(Mtrans, m, d, h);
  list Ms_list;

  for (auto &pair : Ms_full_list)
  {
    auto Ms = pair.second;
    Ms = babaiNP(Ms, B);
    double norm = inf_norm(Ms, Ms.size() - r);
    if (norm <= b) 
    {
      Ms_list.push_back(make_pair(pair.first, Ms));
    }
  }

  return Ms_list;
}