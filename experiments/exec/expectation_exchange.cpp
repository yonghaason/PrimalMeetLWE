#include "functions.h"
#include "cmdline.h"
#include <iostream>
#include <iomanip>
#include <random>

using namespace std;

void experiment(
  bool is_error_unif, size_t R, 
  double ell, double b, size_t r, size_t repeat, double scaler);

int main(int argc, char* argv[]) 
{
  cmdline::parser parser;
  parser.add<uint32_t>("dim", 'r', "Dimension r for the constraint box [-l, l]^r", true, 0);
  parser.add<double>("length", 'l', "Length l for the constraint box [-l, l]^r", true, 0);
  parser.add<bool>("error-shape", 'u', "Is error uniform?", true, 0);
  parser.add<double>("error-param", 'b', "error bound (for unif) or stddev (for gaussian)", true, 0);
  parser.add<double>("multiple", 'C', "# of execution: C * 1/E(p_rep)", false, 1);
  parser.add<uint64_t>("repetition", '\0', "# of execution: R", false, 0);
  parser.add<double>("scaler", '\0', "scale for precision", false, 0);
  parser.add<uint32_t>("experiments", '\0', "# of sampling to estimate E[1-exp(p_rep*R)]", false, 10000);  

  parser.parse_check(argc, argv);
  auto r = parser.get<uint32_t>("dim");
  auto ell = parser.get<double>("length");
  auto is_error_unif = parser.get<bool>("error-shape");
  auto b = parser.get<double>("error-param");
  auto C = parser.get<double>("multiple");
  auto R = parser.get<uint64_t>("repetition");
  auto scaler = parser.get<double>("scaler");
  auto repeat = parser.get<uint32_t>("experiments");

  cout << "**** Params ****" << endl;
  cout << "Constraint box: [" << -ell << ", " << ell << "]^" << r << endl;
  cout << "Error distribution: ";
  if (is_error_unif) cout << "U[" << -b << ", " << b << "]^" << r << endl;
  else cout << "Gaussian with std.dev " << b << "" << endl;
  cout << "---- OR ---- " << endl;
  cout << "LSH-block length: " << 2*ell << endl;
  cout << "Near-collision of distance: ";
  if (is_error_unif) cout << "U[" << -b << ", " << b << "]^" << r << endl;
  else cout << "Gaussian with std.dev " << b << "" << endl;
  cout << endl;

  // Compute E[p_rep], and set R ~ 1/E[p_rep]
  double p_rep_mean_theory = pow(2, scaler);
  for (size_t i = 0; i < r; i++) 
  {
    if (is_error_unif)
    {
      p_rep_mean_theory *= prob_admissible_uniform(b, 2*ell);
    }
    else {
      p_rep_mean_theory *= prob_admissible_gaussian(b, 2*ell);
    }
  }

  cout << "**** Theory ****" << endl;
  cout << "- E[p_rep] : " << p_rep_mean_theory << " (=2^" << log2(p_rep_mean_theory) << ")" <<  endl;
  cout << "- R ~ ";
  if (R != 0) {cout << R << endl;}
  else {cout << C << " / p_rep" << endl;}
  
  if (R == 0) 
  {
    R = C / p_rep_mean_theory;
  }

  cout << "- Expect: 1 - (1 - E[p_rep])^R ~ ";
  cout << 1.0 - pow(1.0 - p_rep_mean_theory, R) << endl;
  cout << endl;

  experiment(is_error_unif, R, ell, b, r, repeat, scaler);

  return 0;
}

void experiment(
  bool is_error_unif, size_t R, 
  double ell, double b, size_t r, size_t repeat, double scaler)
{
  random_device rd;
  mt19937 gen(rd());
  
  double acc = 0.0, acc2 = 0.0;
  double p_rep_mean = 0.0, p_rep_var = 0.0, p_rep_max = 0.0, p_rep_min = 1.0;
  auto normal_sampler = normal_distribution<double>(0, b);
  auto unif_sampler = uniform_real_distribution<double>(-b, b);
  vector<double> error(r);
  for (size_t k = 0; k < repeat; k++) 
  {
    double p_rep = pow(2, scaler);
    for (size_t i = 0; i < r; i++)
    {
      error[i] = is_error_unif ? unif_sampler(gen) : normal_sampler(gen);
      p_rep *= max(0.0, (1.0 - abs(error[i]) / (2*ell)));
    }
    p_rep_max = max(p_rep, p_rep_max);
    p_rep_min = min(p_rep, p_rep_min);
    p_rep_mean += p_rep;
    p_rep_var += p_rep*p_rep;
    acc += pow(1-p_rep, R);
    acc2 += exp(-p_rep*R);
  }
  p_rep_mean /= repeat;
  p_rep_var /= repeat;
  p_rep_var -= p_rep_mean*p_rep_mean;
  acc /= repeat;
  acc2 /= repeat;

  cout << "**** Real ****" << endl;
  cout << "- p_rep range : [" << p_rep_min << ", " << p_rep_max << "] = " 
       << "[2^" << log2(p_rep_min) << ", 2^" << log2(p_rep_max) << "]" << ")" << endl;
  cout << "- E[p_rep] : " << p_rep_mean << " (=2^" << log2(p_rep_mean) << ")" << endl;
  cout << "- E[1 - (1 - p_rep)^R] : " << 1 - acc << " (=2^" << log2(1 - acc) << ")" << endl;
  cout << "- E[exp(R*p_rep)] : " << 1 - acc2 << " (=2^" << log2(1 - acc2) << ")" << endl;
}