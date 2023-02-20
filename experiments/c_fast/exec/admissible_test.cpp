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
  parser.add<double>("box_length", 'b', "Box length", true, 0);
  parser.add<double>("error", 'e', "error param", true, 0);
  parser.add<bool>("unif", 'u', "unif error?", true, 0);
  parser.add<uint32_t>("repeat", 'N', "repetition", false, 1000000);

  parser.parse_check(argc, argv);
  double box_length = parser.get<double>("box_length");
  double error_param = parser.get<double>("error");
  bool unif = parser.get<bool>("unif");
  uint32_t repeat = parser.get<uint32_t>("repeat");

  auto box_sampler = uniform_real_distribution<double>(-box_length/2, box_length/2);
  uint32_t cnt = 0;
  if (!unif) {
    auto error_sampler = normal_distribution<double>(0, error_param);
    for (size_t i = 0; i < repeat; i++) {
      auto x = box_sampler(gen);
      auto e = error_sampler(gen);
      if (abs(x + e) < box_length/2) cnt++;
    }
    cout << "Gaussian error test." << endl;
    cout << "Theoretical expectation: " << prob_admissible_gaussian(error_param, box_length) << endl;
    cout << "Experiment: " << (double) cnt / (double) repeat << endl;
  }
  else {
    auto error_sampler = uniform_real_distribution<double>(-error_param, error_param);
    for (size_t i = 0; i < repeat; i++) {
      auto x = box_sampler(gen);
      auto e = error_sampler(gen);
      if (abs(x + e) < box_length/2) cnt++;
    }
    cout << "Uniform error test." << endl;
    cout << "Theoretical expectation: " << prob_admissible_uniform(error_param, box_length) << endl;
    cout << "Experiment: " << (double) cnt / (double) repeat << endl;
  }

  return 0;
}