#include "functions.h"
#include "cmdline.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <map>

// ./cons_test -d 25 -h 8 -w 5 -e 0.75 -b 1.0 -u 0

using namespace std;

int main(int argc, char* argv[]) {
  
  chrono::system_clock::time_point start, end;

  uint32_t d, h; 
  uint64_t q;
  double error_param; 
  double constraint_bound, lsh_length;
  uint32_t guess_weight;
  bool unif;
  uint32_t target_rep_num;

  uint32_t repetition;
  uint32_t print_term;

  cmdline::parser parser;

  parser.add<uint32_t>("d", 'd', "Dimension (r^(i-1) in paper)", true, 20);
  parser.add<uint32_t>("h", 'h', "Secret weight (w^(i-1) in paper)", true, 0);
  parser.add<uint32_t>("guess_weight", 'w', "Splitting weight (w^(i) in paper)", true, 0);
  parser.add<bool>("unif", 'u', "lower list elts sampled from uniform dist?", true, true);
  parser.add<double>("error_param", 'e', "Error stddev (for gaussian) or bound (unif)", true, 3.2);
  parser.add<double>("constraint_bound", 'l', "Constraint bound (ell^(i) in paper)", true, 0);
  parser.add<uint32_t>("target_rep_num", 't', "# of target representation number", false, 1.0);  
  parser.add<uint64_t>("q", 'q', "Modulus param that determines the domain", false, 2048);

  parser.add<uint32_t>("repeat", '\0', "Number of total experiments", false, 100);
  parser.add<uint32_t>("print_term", '\0', "Print avgs after THIS_VALUE experiments", false, 50);
  
  parser.parse_check(argc, argv);

  d = parser.get<uint32_t>("d");
  h = parser.get<uint32_t>("h");
  q = parser.get<uint64_t>("q");
  error_param = parser.get<double>("error_param");
  guess_weight = parser.get<uint32_t>("guess_weight");
  constraint_bound = parser.get<double>("constraint_bound");
  unif = parser.get<bool>("unif");
  target_rep_num = parser.get<uint32_t>("target_rep_num");

  repetition = parser.get<uint32_t>("repeat");
  print_term = parser.get<uint32_t>("print_term");

  auto m = d;
  
  ConstraintTest cons_test(d, error_param, q, h, m, unif);

  cons_test.set_constraint_dim(guess_weight, constraint_bound, target_rep_num);
  
  cout << "**** Start experiments with random M and s ****" << endl;
  matrix M; secret s; vector<double> e;
  cons_test.gen_noisy_instance(M, s, e);

  auto cnt = cons_test.new_check_near_collision(M, s, e, guess_weight, constraint_bound);
  double avg_sol = cnt;
  
  //  Actual List Construction & Collision Finding
  for (size_t i = 1; i < repetition; i++) {
    cons_test.gen_noisy_instance(M, s, e);

    auto cnt = cons_test.new_check_near_collision(M, s, e, guess_weight, constraint_bound);
    
    avg_sol = (avg_sol * i + (double) cnt) / (i+1);
    if ((i+1) % print_term == 0) {
      cout << setw(3) << "- " << i+1 << " executions: E(# Pairs) = " << avg_sol << endl;
    }
  }

  return 0;
}