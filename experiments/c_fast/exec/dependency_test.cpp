#include "functions.h"
#include "cmdline.h"
#include <iostream>
#include <iomanip>
#include <random>

using namespace std;

int main(int argc, char* argv[]) {

  random_device rd;
  mt19937 gen(rd());

  cmdline::parser parser;
  // parser.add<double>("coord-length", 'q', "Last coordinate length", true, 0);
  // parser.add<double>("rhf", '\0', "Coordinate length ratio (GSA)", false, 1.12);  
  parser.add<uint32_t>("domain-dim", 'd', "Dimension d of the domain [-l, l]^d", true, 0);
  parser.add<double>("domain-len", 'l', "Length of the domain [-l, l]^d", true, 0);
  parser.add<uint32_t>("cons-dim", 'r', "Dimension r of the constraint box [-c, c]^r", true, 0);
  parser.add<double>("cons-len", 'c', "Length for the constraint box [-c, c]^r", true, 0);
  parser.add<bool>("error-shape", 'u', "Is error uniform?", true, 0);
  parser.add<double>("error-param", 'b', "error bound (for unif) or stddev (for gaussian)", true, 0);
  parser.add<uint32_t>("repeat", '\0', "# of sampling to estimate E[(1-p_rep)^R]", false, 10000);

  parser.parse_check(argc, argv);
  // auto q = parser.get<double>("coord-length");
  // auto rhf = parser.get<double>("rhf");
  auto d = parser.get<uint32_t>("domain-dim");
  auto l = parser.get<double>("domain-len");
  auto r = parser.get<uint32_t>("cons-dim");
  auto c = parser.get<double>("cons-length");
  auto is_error_unif = parser.get<bool>("error-shape");
  auto b = parser.get<double>("error-param");
  uint32_t repeat = parser.get<uint32_t>("repeat");

  // vector<double> GSnorm(r);
  // vector<double> cons_len(r);
  // GSnorm[r-1] = q;
  // for (int i = r-2; i >= 0; i--)
  // {
  //   GSnorm[i] = GSnorm[i+1] * rhf;
  // }
  // for (size_t i = 0; i < r; i++) 
  // {
  //   2*c = (GSnorm[i] > 2*c) ? 2*c : GSnorm[i];
  // }

  cout << "**** Params ****" << endl;
  cout << "Error distribution: ";
  if (is_error_unif) cout << "U[" << -b << ", " << b << "]" << endl;
  else cout << "Gaussian with std.dev " << b << "" << endl;
  cout << endl;

  double vol_ratio = 1;
  // for (size_t i = 0; i < r; i++) 
  // {
  //   vol_ratio *= 2*c / l;
  // }

  // Experiment 1. Probability of v, v + e
  vector<double> error(d);
  auto unif_sampler = uniform_real_distribution<double>(-l, l);
  auto normal_sampler = normal_distribution<double>(0, b);
  for (size_t i = 0; i < d; i++) 
  {
    error[i] = is_error_unif ? unif_sampler(gen) : normal_sampler(gen);
  }

  for (size_t iter = 0; iter < repeat; iter++)
  {
    vector<double> v1(d), v2(d);
    auto unif_sampler = uniform_real_distribution<double>(-l, l);
    for (size_t i = 0; i < d; i++) 
    {
      v1[i] = unif_sampler(gen);
      v2[i] = v1[i] + error[i];
    }

    auto p_rep1 = vol_ratio, p_rep2 = vol_ratio;

    for (size_t i = 0; i < r; i++)
    {
      p_rep1 *= max(0.0, 1.0 - abs(v1[i]) / (2*c));
      p_rep2 *= max(0.0, 1.0 - abs(v2[i]) / (2*c));
    }
    auto p_sp1 = 1 - pow(1 - p_rep1, 10);
    auto p_sp2 = 1 - pow(1 - p_rep2, 10);
  }

  // Experiment 2. Probability of v1, v2 (indep)
  for (size_t iter = 0; iter < repeat; iter++)
  {
    vector<double> v1(d), v2(d);
    auto unif_sampler = uniform_real_distribution<double>(-l, l);
    for (size_t i = 0; i < d; i++) 
    {
      v1[i] = unif_sampler(gen);
      v2[i] = unif_sampler(gen);
    }

    auto p_rep1 = vol_ratio, p_rep2 = vol_ratio;

    for (size_t i = 0; i < r; i++)
    {
      p_rep1 *= max(0.0, 1.0 - abs(v1[i]) / (2*c));
      p_rep2 *= max(0.0, 1.0 - abs(v2[i]) / (2*c));
    }
    auto p_sp1 = 1 - pow(1 - p_rep1, 10);
    auto p_sp2 = 1 - pow(1 - p_rep2, 10);
  }  

  return 0;
}