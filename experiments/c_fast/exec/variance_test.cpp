#include "functions.h"
#include "cmdline.h"
#include <iostream>
#include <iomanip>
#include <random>

using namespace std;

int main(int argc, char* argv[]) {

  random_device rd;
  mt19937 gen(rd());

  uint64_t repeat = 10000;

  double mean, variance, avg, unif_error_bnd;
  uint64_t meaninv;

  auto normal_sampler = normal_distribution<double>(1, 1);
  auto unif_sampler = uniform_real_distribution<double>(0, 1);

  mean = 0.01;
  variance = 0.001;
  meaninv = 1 / mean;

  normal_sampler = normal_distribution<double>(mean, variance);
  
  avg = 0;
  for (size_t i = 0; i < repeat; i++) {
    auto x = normal_sampler(gen);
    avg += pow(1-x, meaninv);
  }
  avg /= repeat;

  cout << "E[(1-p)^(1/E(p))]: " << avg 
  << " / (1 - E(p))^(1/E(p)): " << pow(1 - mean, meaninv) << endl;


  mean = 0.01;
  unif_error_bnd = 0.009;
  meaninv = 1 / mean;

  unif_sampler = uniform_real_distribution<double>(mean-unif_error_bnd, mean+unif_error_bnd);
  
  avg = 0;
  for (size_t i = 0; i < repeat; i++) {
    auto x = unif_sampler(gen);
    avg += pow(1-x, meaninv);
  }
  avg /= repeat;

  cout << "E[(1-q)^(1/E(q))]: " << avg 
  << " / (1 - E(q))^(1/E(q)): " << pow(1 - mean, meaninv) << endl;

  mean = 0.001;
  unif_error_bnd = 0.0009;
  meaninv = 1 / mean;

  unif_sampler = uniform_real_distribution<double>(mean-unif_error_bnd, mean+unif_error_bnd);
  
  avg = 0;
  for (size_t i = 0; i < repeat; i++) {
    auto x = unif_sampler(gen);
    avg += pow(1-x, meaninv);
  }
  avg /= repeat;

  cout << "E[(1-q)^(1/E(q))]: " << avg 
  << " / (1 - E(q))^(1/E(q)): " << pow(1 - mean, meaninv) << endl;

  return 0;
}