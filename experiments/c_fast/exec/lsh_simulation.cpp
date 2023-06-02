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
  parser.add<double>("coord-length", 'q', "Last coordinate length", true, 0);
  parser.add<double>("rhf", '\0', "Coordinate length ratio (GSA)", false, 1.12);
  parser.add<uint32_t>("lsh-dim", 'r', "Dimension r of LSH domain", true, 0);
  parser.add<double>("lsh-box-length", 'l', "Length l of lsh box", true, 0);
  parser.add<bool>("error-shape", 'u', "Is error uniform?", true, 0);
  parser.add<double>("error-param", 'b', "error bound (for unif) or stddev (for gaussian)", true, 0);
  parser.add<int32_t>("multiple", '\0', "# of execution: C * 1/E(p_rep)", false, 1);
  parser.add<uint32_t>("repeat", '\0', "# of sampling to estimate E[(1-p_rep)^R]", false, 10000);

  parser.parse_check(argc, argv);
  auto q = parser.get<double>("coord-length");
  auto rhf = parser.get<double>("rhf");
  auto r = parser.get<uint32_t>("lsh-dim");
  auto ell = parser.get<double>("lsh-box-length");
  auto is_error_unif = parser.get<bool>("error-shape");
  auto b = parser.get<double>("error-param");
  auto C = parser.get<int32_t>("multiple");
  uint32_t repeat = parser.get<uint32_t>("repeat");

  vector<double> GSnorm(r);
  vector<double> lsh_len(r);
  GSnorm[r-1] = q;
  for (int i = r-2; i >= 0; i--)
  {
    GSnorm[i] = GSnorm[i+1] * rhf;
  }
  for (size_t i = 0; i < r; i++) 
  {
    size_t box_num = ceil(GSnorm[i] / ell);
    lsh_len[i] = GSnorm[i] / box_num;
  }

  cout << "**** Params ****" << endl;
  cout << "Coordinate Length: [" << GSnorm[0] << ", ... , " 
      << GSnorm[r-1] << "] (Geometrically drops with ratio " << rhf << ")" << endl;
  // cout << "LSH box lengths: ";
  // for (size_t i = 0; i < r-1; i++) {cout << lsh_len[i] << ", ";}
  // cout << lsh_len[r-1] << endl;
  cout << "LSH dim: " << r << " / box length: " << ell << endl;
  cout << "Error distribution: ";
  if (is_error_unif) cout << "U[" << -b << ", " << b << "]" << endl;
  else cout << "Gaussian with std.dev " << b << "" << endl;
  cout << endl;

  // Compute E[p_good], and set R ~ 1/E[p_good]
  double p_good_mean_theory = 1.0;
  for (size_t i = 0; i < r; i++) 
  {
    if (is_error_unif)
    {
      p_good_mean_theory *= prob_admissible_uniform(b, lsh_len[i]);
    }
    else {
      p_good_mean_theory *= prob_admissible_gaussian(b, lsh_len[i]);
    }
  }
  size_t R = C / p_good_mean_theory;

  cout << "**** Theory ****" << endl;
  cout << "- E[p_good] : " << p_good_mean_theory << " (=2^" << log2(p_good_mean_theory) << ")" 
       << ", R = " << C << " / E[p_good]" << endl;
  cout << "- Expect: 1 - (1 - E[p_good])^R ~ 1 - exp(-" << C << ") = " << 1 - exp((double) -C) << endl;
  cout << endl;
  
  // Estimate E[(1 - p_good)^R] 
  double acc = 0.0;
  double p_good_mean = 0.0, p_good_var = 0.0, p_good_max = 0.0, p_good_min = 1.0;
  auto normal_sampler = normal_distribution<double>(0, b);
  auto unif_sampler = uniform_real_distribution<double>(-b, b);
  vector<double> error(r);
  for (size_t k = 0; k < repeat; k++) 
  {
    double p_good = 1.0;
    for (size_t i = 0; i < r; i++)
    {
      error[i] = is_error_unif ? unif_sampler(gen) : normal_sampler(gen);
      p_good *= max(0.0, (1.0 - abs(error[i]) / lsh_len[i]));
    }
    p_good_max = max(p_good, p_good_max);
    p_good_min = min(p_good, p_good_min);
    p_good_mean += p_good;
    p_good_var += p_good*p_good;
    acc += pow(1-p_good, R);
  }
  p_good_mean /= repeat;
  p_good_var /= repeat;
  p_good_var -= p_good_mean*p_good_mean;
  acc /= repeat;

  cout << "**** Real ****" << endl;
  cout << "- p_good range : [" << p_good_min << ", " << p_good_max << "] = " 
       << "[2^" << log2(p_good_min) << ", 2^" << log2(p_good_max) << "]" << ")" << endl;
  cout << "- E[p_good] : " << p_good_mean << " (=2^" << log2(p_good_mean) << ")" << endl;
  cout << "- s[p_good] : " << sqrt(p_good_var) << " (=2^" << 0.5*log2(p_good_var) << ")" << endl;
  cout << "- E[1 - (1 - p_good)^R] : " << 1 - acc << " (=2^" << log2(1 - acc) << ")" << endl;  

  return 0;
}