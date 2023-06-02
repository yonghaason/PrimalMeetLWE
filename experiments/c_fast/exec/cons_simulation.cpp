#include "functions.h"
#include "cmdline.h"
#include <iostream>
#include <iomanip>
#include <random>

using namespace std;

void experiment(
  bool is_error_unif, double vol_ratio, size_t R, 
  vector<double> cons_len, double b, size_t r, 
  double scale, size_t repeat);

int main(int argc, char* argv[]) 
{
  cmdline::parser parser;
  parser.add<double>("coord-length", 'q', "Last coordinate length", true, 0);
  parser.add<double>("rhf", '\0', "Coordinate length ratio (GSA)", false, 1.12);
  parser.add<uint32_t>("cons-dim", 'r', "Dimension r for the constraint box [-l, l]^r", true, 0);
  parser.add<double>("cons-length", 'l', "Length l for the constraint box [-l, l]^r", true, 0);
  parser.add<bool>("error-shape", 'u', "Is error uniform?", true, 0);
  parser.add<double>("error-param", 'b', "error bound (for unif) or stddev (for gaussian)", true, 0);
  parser.add<int32_t>("multiple", '\0', "# of execution: C * 1/E(p_rep)", false, 1);
  parser.add<uint32_t>("repeat", '\0', "# of sampling to estimate E[(1-p_rep)^R]", false, 10000);
  parser.add<int32_t>("scaler", '\0', "scale for precision", false, 25);

  parser.parse_check(argc, argv);
  auto q = parser.get<double>("coord-length");
  auto rhf = parser.get<double>("rhf");
  auto r = parser.get<uint32_t>("cons-dim");
  auto ell = parser.get<double>("cons-length");
  auto is_error_unif = parser.get<bool>("error-shape");
  auto b = parser.get<double>("error-param");
  auto C = parser.get<int32_t>("multiple");
  auto scaler = parser.get<int32_t>("scaler");
  uint32_t repeat = parser.get<uint32_t>("repeat");

  vector<double> GSnorm(r);
  vector<double> cons_len(r);
  GSnorm[r-1] = q;
  for (int i = r-2; i >= 0; i--)
  {
    GSnorm[i] = GSnorm[i+1] * rhf;
  }
  for (size_t i = 0; i < r; i++) 
  {
    cons_len[i] = (GSnorm[i] > 2*ell) ? 2*ell : GSnorm[i];
  }

  cout << "**** Params ****" << endl;
  cout << "Coordinate Length: [" << GSnorm[0] << ", ... , " 
      << GSnorm[r-1] << "] (Geometrically drops with ratio " << rhf << ")" << endl;
  // cout << "Constraint box coordinate lengths: ";
  // for (size_t i = 0; i < r-1; i++) {cout << cons_len[i] << ", ";}
  // cout << cons_len[r-1] << endl;
  cout << "Constraint box: [" << -ell << ", " << ell << "]^" << r << endl;
  cout << "Error distribution: ";
  if (is_error_unif) cout << "U[" << -b << ", " << b << "]" << endl;
  else cout << "Gaussian with std.dev " << b << "" << endl;
  cout << endl;

  double vol_ratio = 1;
  for (size_t i = 0; i < r; i++) 
  {
    if (GSnorm[i] > 2*ell)
    {
      vol_ratio *= 2*ell / GSnorm[i];
    }
  }

  // Compute E[p_rep], and set R ~ 1/E[p_rep]
  double p_rep_mean_theory = vol_ratio;
  for (size_t i = 0; i < r; i++) 
  {
    if (is_error_unif)
    {
      p_rep_mean_theory *= prob_admissible_uniform(b, cons_len[i]);
    }
    else {
      p_rep_mean_theory *= prob_admissible_gaussian(b, cons_len[i]);
    }
  }

  cout << "**** Theory ****" << endl;
  // cout << "- Volume ratio : " << vol_ratio << " (=2^" << log2(vol_ratio) << ")" << endl;
  // cout << "- E[p_rep] : " << p_rep_mean_theory << " (=2^" << log2(p_rep_mean_theory) << ")"
  //      << ", R = " << C << " / E[p_rep]" << endl;
  cout << "- E[p_rep/vol_ratio] : " << p_rep_mean_theory/vol_ratio << " (=2^" << log2(p_rep_mean_theory/vol_ratio) << ")" <<  endl;
  cout << "- R ~ " << C << " / p_rep" << endl; 
  cout << "- Expect: 1 - (1 - E[p_rep])^R ~ 1 - exp(-" << C << ") = " << 1 - exp((double) -C) << endl;
  cout << endl;

  double scale = pow(2, -scaler) / p_rep_mean_theory;
  vol_ratio *= scale;
  p_rep_mean_theory *= scale;

  size_t R = C / p_rep_mean_theory;

  experiment(is_error_unif, vol_ratio, R, cons_len, b, r, scale, repeat);

  return 0;
}

void experiment(
  bool is_error_unif, double vol_ratio, size_t R, 
  vector<double> cons_len, double b, size_t r, 
  double scale, size_t repeat)
{
  random_device rd;
  mt19937 gen(rd());
  
  // Estimate E[(1 - p_rep)^R] 
  double acc = 0.0;
  double p_rep_mean = 0.0, p_rep_var = 0.0, p_rep_max = 0.0, p_rep_min = 1.0;
  auto normal_sampler = normal_distribution<double>(0, b);
  auto unif_sampler = uniform_real_distribution<double>(-b, b);
  vector<double> error(r);
  for (size_t k = 0; k < repeat; k++) 
  {
    double p_rep = vol_ratio;
    for (size_t i = 0; i < r; i++)
    {
      error[i] = is_error_unif ? unif_sampler(gen) : normal_sampler(gen);
      p_rep *= max(0.0, (1.0 - abs(error[i]) / cons_len[i]));
    }
    p_rep_max = max(p_rep, p_rep_max);
    p_rep_min = min(p_rep, p_rep_min);
    p_rep_mean += p_rep;
    p_rep_var += p_rep*p_rep;
    acc += pow(1-p_rep, R);
  }
  p_rep_mean /= repeat;
  p_rep_var /= repeat;
  p_rep_var -= p_rep_mean*p_rep_mean;
  acc /= repeat;

  cout << "**** Real ****" << endl;
  // cout << "- p_rep range : [" << p_rep_min << ", " << p_rep_max << "] = " 
      //  << "[2^" << log2(p_rep_min) << ", 2^" << log2(p_rep_max) << "]" << endl;
  cout << "- p_rep/vol_ratio range : [" << p_rep_min/vol_ratio << ", " << p_rep_max/vol_ratio << "] = " 
       << "[2^" << log2(p_rep_min/vol_ratio) << ", 2^" << log2(p_rep_max/vol_ratio) << "]" << ")" << endl;
  // cout << "- E[p_rep] : " << p_rep_mean << " = 2^" << log2(p_rep_mean) << endl;
  cout << "- E[p_rep/vol_ratio] : " << p_rep_mean/vol_ratio << " (=2^" << log2(p_rep_mean/vol_ratio) << ")" << endl;
  cout << "- s[p_rep/vol_ratio] : " << sqrt(p_rep_var)/vol_ratio 
       << " (=2^" << 0.5*log2(p_rep_var) - log2(vol_ratio) << ")" << endl;
  cout << "- E[1 - (1 - p_rep)^R] : " << 1 - acc << " (=2^" << log2(1 - acc) << ")" << endl;
}