#include "functions.h"
#include "cmdline.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <map>

// ./cons_test -n 25 -m 30 -h 8 -w 5 -q 2048 -e 0.75 -b 1.0 -u 0

using namespace std;

int main(int argc, char* argv[]) {
  
  chrono::system_clock::time_point start, end;

  uint32_t n, m, h; 
  uint64_t q;
  double error_param; 
  double constraint_bound, lsh_length;
  uint32_t guess_weight;
  bool unif;
  uint32_t repetition;

  uint32_t print_term;

  cmdline::parser parser;

  parser.add<uint32_t>("n", 'n', "LWE dimension", true, 20);
  parser.add<uint32_t>("m", 'm', "# of LWE samples", false, 0);
  parser.add<uint32_t>("h", 'h', "Secret HW", true, 0);
  parser.add<uint64_t>("q", 'q', "LWE modulus", true, 1024);
  parser.add<double>("error_param", 'e', "LWE error e stddev (gaussian) or bound (unif)", true, 3.2);
  parser.add<uint32_t>("guess_weight", 'w', "List-1 HW", true, 0);
  parser.add<double>("constraint_bound", 'b', "Constraint bound", true, 0);
  parser.add<bool>("unif", 'u', "Error unif?", true, true);
  parser.add<uint32_t>("repeat", '\0', "Number of experiments", false, 100);
  parser.add<uint32_t>("print_term", '\0', "Print avgs term", false, 50);
  
  parser.parse_check(argc, argv);

  n = parser.get<uint32_t>("n");
  m = parser.get<uint32_t>("m");
  h = parser.get<uint32_t>("h");
  q = parser.get<uint64_t>("q");
  error_param = parser.get<double>("error_param");
  guess_weight = parser.get<uint32_t>("guess_weight");
  constraint_bound = parser.get<double>("constraint_bound");
  unif = parser.get<bool>("unif");
  repetition = parser.get<uint32_t>("repeat");

  print_term = parser.get<uint32_t>("print_term");

  if (m == 0) m = n;
  
  ConstraintTest cons_test(n, error_param, q, h, m, unif);

  cons_test.set_constraint_dim(guess_weight, constraint_bound);
  // auto partial_secrets = cons_test.enumerate_secrets(n, guess_weight);
  
  cout << "**** Start experiments with random M and s ****" << endl;
  matrix M; secret s; vector<double> e;
  cons_test.gen_noisy_instance(M, s, e);
  // auto sol = cons_test.fast_check_near_collision(M, s, e, partial_secrets, guess_weight, constraint_bound);
  // double avg_sol = sol.size();

  auto cnt = cons_test.new_check_near_collision(M, s, e, guess_weight, constraint_bound);
  double avg_sol = cnt;
  
  //  Actual List Construction & Collision Finding
  map<size_t, size_t> stats_sol;
  for (size_t i = 1; i < repetition; i++) {
    cons_test.gen_noisy_instance(M, s, e);
    // auto sol = cons_test.fast_check_near_collision(M, s, e, partial_secrets, guess_weight, constraint_bound);
    // auto cnt = sol.size();

    auto cnt = cons_test.new_check_near_collision(M, s, e, guess_weight, constraint_bound);
    
    // assert(sol.size() == cnt);

    avg_sol = (avg_sol * i + (double) cnt) / (i+1);
    if ((i+1) % print_term == 0) {
      cout << setw(3) << "- " << i+1 << " executions: E(|L_0|) = " << avg_sol << endl;
    }
    auto filled = stats_sol.count(cnt);
    if (filled == 1) {stats_sol[cnt]++;}
    else {stats_sol.insert({cnt, 1});}
  }

  cout << "\n**** Stats ****" << endl;
  for (auto stat: stats_sol) {
    cout << "- \'|L_0| = " << stat.first << "\' occurs " << stat.second << " times" << endl;
  }

  return 0;
}