#include "functions.h"
#include "cmdline.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <map>
#include <random>

using namespace std;

// ./ncf_near_col_test -m 15 -d 10 -h 4 -r 6 -b 4 -q 8 --repeat 1000

void experiment(
  matrix& M, uint32_t m, uint32_t d, uint32_t h, matrix& B, uint32_t b, uint32_t l, uint32_t r,
  double& pt_ratio, double& pt_ratio_sq, double& nc_count, double& nc_count_sq);

int main(int argc, char* argv[]) {
  
  chrono::system_clock::time_point start, end;

  cmdline::parser parser;

  parser.add<uint32_t>("m", 'm', "matrix dimension", true, 0);
  parser.add<uint32_t>("d", 'd', "secret dimension", true, 20);
  parser.add<uint32_t>("h", 'h', "secret weight", true, 0);
  parser.add<uint32_t>("r", 'r', "projection dimension", true, 0);
  parser.add<double>("b", 'b', "box length [-b, b]", true, 0);
  parser.add<double>("l", 'l', "near-collision bound ell", false, 0);
  parser.add<double>("stddev", '\0', "error stddev", false, 1);
  parser.add<double>("q", 'q', "modulus param", false, 2048);
  parser.add<double>("rhf", '\0', "root-hermite-factor", false, 1.05);

  parser.add<uint64_t>("repeat", '\0', "# of experiments", false, 1000);

  parser.parse_check(argc, argv);

  auto m = parser.get<uint32_t>("m");
  auto d = parser.get<uint32_t>("d");
  auto h = parser.get<uint32_t>("h");
  auto r = parser.get<uint32_t>("r");
  auto b = parser.get<double>("b");
  auto l = parser.get<double>("l");
  auto stddev = parser.get<double>("stddev");
  auto q0 = parser.get<double>("q");
  auto rhf = parser.get<double>("rhf");
  auto repeat = parser.get<uint64_t>("repeat");

  if (l == 0) l = b;

  vector<double> q(m);
  for (int i = 0; i < m; i++) {
    q[m-i-1] = q0 * pow(rhf, i);
  }

  cout << "*********** Experiments for Heuristic 4 *********" << endl; 
  cout << "* Settings " << endl;
  cout << "- B = " << m << " X " << m << " matrix of Gram-Schmidt norm " 
            << q[0] << " -> " << q[m-1] 
            << " (GSA with rhf = " << rhf << ")" << endl;
  cout << "- D = " << m << "-dimensional cube of i-th coord length " 
       << "q[i] if i < " << m-r << " and " << 2*b << " if i >= " << m-r << endl;
  cout << endl;  
  cout << "* Goal: Compare # of " << l << "-near-collisions in L X L defined as L = { [Ms]_B: HW(s) = " << h << " & [Ms]_B ∈ D }" << endl; 
  cout << " where M = our definition v.s. truly uniform" << endl;
  cout << endl;

  double p_vol = 1;
  for (int i = m-1; i >= m-r; i--)
  {
    if (q[i] > 2*b) 
    {
      p_vol *= (2*b / q[i]);
    }
  }  

  double p_nc = 1.0;
  for (int i = 0; i < m; i++)
  {
    if (i < m-r)
    {
      if (q[i] > l) 
        p_nc *= l/q[i]*(2 - l/q[i]);
    }
    else {
      auto len = min(2*b, q[i]);
      if (len > l)
        p_nc *= l/len*(2 - l/len);
    }
  }

  double pt_ratio_noisy = 0.0;
  double pt_ratio_sq_noisy = 0.0;
  double nc_count_noisy = 0;
  double nc_count_sq_noisy = 0;
  double pt_ratio_unif = 0.0;
  double pt_ratio_sq_unif = 0.0;
  double nc_count_unif = 0;
  double nc_count_sq_unif = 0;

  size_t full_size = binom(d, h) * (1ull << h);
  double nc_expect = full_size * full_size * p_vol * p_vol * p_nc / 2;    
  
  for (size_t iter = 0; iter < repeat; iter++) 
  {
    matrix M; matrix B; secret s; vector<double> e;
    gen_noisy_instance(
      m, d, d-1, stddev, q, 
      M, B, s, e);

    experiment(
      M, m, d, h, B, b, l, r,
      pt_ratio_noisy, pt_ratio_sq_noisy, nc_count_noisy, nc_count_sq_noisy);

    matrix A;
    gen_uniform_matrix(m, d, q, A);

    experiment(
      A, m, d, h, B, b, l, r,
      pt_ratio_unif, pt_ratio_sq_unif, nc_count_unif, nc_count_sq_unif);
    
    cout << "... Done " << iter << "-th exec. Current avg: " << nc_count_noisy / (iter+1) 
          << " v.s. " << nc_count_unif / (iter+1) << " (theory = " << nc_expect << ")"
          // << endl;
          << "\r" << std::flush;    
  }
  cout << endl;

  cout << "* Results"  << endl;
  auto nc_count_avg_noisy = nc_count_noisy / repeat;
  auto nc_count_avg_unif = nc_count_unif / repeat;
  auto nc_count_var_noisy = nc_count_sq_noisy / repeat - nc_count_avg_noisy*nc_count_avg_noisy;
  auto nc_count_var_unif = nc_count_sq_unif / repeat - nc_count_avg_unif*nc_count_avg_unif;
  cout << "# " << l << "-Near-Cols. noisy, unif, expec: "
        << nc_count_avg_noisy << ", " << nc_count_avg_unif << ", " << nc_expect << endl;
  cout << "             (std.dev) "
        << sqrt(nc_count_var_noisy) << ", " << sqrt(nc_count_var_unif) << endl;

  return 0;
}

void experiment(
  matrix& M, uint32_t m, uint32_t d, uint32_t h, matrix& B, uint32_t b, uint32_t l, uint32_t r,
  double& pt_ratio, double& pt_ratio_sq, double& nc_count, double& nc_count_sq)
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

  auto cur_ratio = (double) Ms_list.size() / (double) Ms_full_list.size();
  pt_ratio += cur_ratio;
  pt_ratio_sq += cur_ratio*cur_ratio;

  size_t cur_count = 0;
  for (size_t i = 0; i < Ms_list.size(); i++)
  {
    auto e1 = Ms_list[i].second;
    for (size_t j = i+1; j < Ms_list.size(); j++)
    {
      auto e2 = Ms_list[j].second;
      auto ss = subvec(e1, e2);
      if (inf_norm(ss) <= l) cur_count++;
    }
  }

  nc_count += cur_count;
  nc_count_sq += cur_count*cur_count;
}