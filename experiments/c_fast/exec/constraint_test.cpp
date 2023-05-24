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
  parser.add<uint32_t>("print_term", '\0', "Print avgs after THIS_VALUE experiments", false, 100);
  
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
  
  // cout << "**** Start experiments with random M and s ****" << endl;
  matrix M; secret s; vector<double> e;
  cons_test.gen_noisy_instance(M, s, e);

  auto p_rep = prob_admissible_fixed(e, 2*constraint_bound, cons_test.proj_dim);
  p_rep /= cons_test.vol_ratio;
  double p_rep_max = p_rep, p_rep_min = p_rep, p_rep_mean = p_rep, p_rep_var = p_rep*p_rep;
  
  auto num_pairs = cons_test.new_check_near_collision(M, s, e, guess_weight, constraint_bound);

  double num_pairs_avg = num_pairs;

  cout << fixed;
	cout.precision(4);
  
  //  Actual List Construction & Collision Finding
  cout << "**** Reals ****" << endl;
  map<size_t, size_t> num_pairs_stats;
  size_t failures = 0;
  for (size_t i = 1; i < repetition; i++) {
    cons_test.gen_noisy_instance(M, s, e);

    p_rep = prob_admissible_fixed(e, 2*constraint_bound, cons_test.proj_dim);
    p_rep /= cons_test.vol_ratio;
    p_rep_max = max(p_rep, p_rep_max);
    p_rep_min = min(p_rep, p_rep_min);
    p_rep_mean += p_rep;
    p_rep_var += p_rep*p_rep;
    
    auto num_pairs = cons_test.new_check_near_collision(M, s, e, guess_weight, constraint_bound);
    
    num_pairs_avg += num_pairs;
    if ((i+1) % print_term == 0) {
      cout << "Avg of " << setw(5) << i+1 
      << " execs: Pr[# pairs = 0] = " << (double) failures / (double) (i+1) 
      << " / # Pairs = " << num_pairs_avg / (i+1) << endl;
    }
    auto filled = num_pairs_stats.count(num_pairs);
    if (filled == 1) {num_pairs_stats[num_pairs]++;}
    else {num_pairs_stats.insert({num_pairs, 1});}

    if (num_pairs == 0) failures++;
  }

  cout << "\n**** Some Stats ****" << endl;
  p_rep_mean = p_rep_mean / repetition;
  p_rep_var = p_rep_var / repetition - p_rep_mean * p_rep_mean;

  cout << "p_reps lies in [" << p_rep_min << ", " << p_rep_max << "] with E(p_rep): " << p_rep_mean << ", V(p_rep): " << p_rep_var << endl;

  // for (auto stat: num_pairs_stats) {
  //   cout << "- \'|L_0| = " << stat.first << "\' occurs " << stat.second << " times" << endl;
  // }

  return 0;
}