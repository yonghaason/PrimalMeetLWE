#include "functions.h"
#include "cmdline.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <map>

// Example
// ./nc_test -r 25 -q 4096 -e 3 -b 6 -d 10 -u 1

using namespace std;

int main(int argc, char* argv[]) {

  chrono::system_clock::time_point start, end;

  uint32_t r; 
  uint64_t q;
  double e;
  double lsh_length;
  uint32_t lsh_dim;
  uint32_t iteration_multiple;
  uint64_t N;
  bool unif;

  uint32_t repetition;
  uint32_t print_term;

  cmdline::parser parser;

  parser.add<uint32_t>("r", 'r', "Dimension of list vectors", true, 20);
  parser.add<double>("stddev", 'e', "error stddev (gaussian) or bound (unif)", true, 3.2);
  parser.add<double>("lsh_length", 'b', "LSH box length", true, 0);
  parser.add<uint32_t>("lsh_dim", 'd', "LSH dimension", true, 10);
  parser.add<uint32_t>("lsh_iter_mult", '\0', "Iterate THIS_VALUE / p_good times ", false, 3);
  parser.add<bool>("unif", 'u', "Error unif?", true, true);
  parser.add<uint64_t>("list_size", 'N', "List size", false, 0);
  parser.add<uint64_t>("q", 'q', "Domain param", false, 4096);

  parser.add<uint32_t>("repeat", '\0', "# of total experiments", false, 100);
  parser.add<uint32_t>("print_term", '\0', "Print avgs after THIS_VALUE experiments", false, 50);
  
  parser.parse_check(argc, argv);

  r = parser.get<uint32_t>("r");
  q = parser.get<uint64_t>("q");
  e = parser.get<double>("stddev");
  lsh_length = parser.get<double>("lsh_length");
  lsh_dim = parser.get<uint32_t>("lsh_dim");
  N = parser.get<uint64_t>("list_size");
  iteration_multiple = parser.get<uint32_t>("lsh_iter_mult");
  unif = parser.get<bool>("unif");
  
  repetition = parser.get<uint32_t>("repeat");
  print_term = parser.get<uint32_t>("print_term");
  
  if (unif && (lsh_length < e)) {
    cout << "!!! Some error can exceed LSH length !!!\n"
         << "(= LSH never finds some near-collision pairs.)" << endl;
  }
  if (!unif && (lsh_length < 3*e)) { // some arbitrary choice
    cout << "!!! Some error can exceed LSH length !!!\n"
         << "(= LSH never finds some near-collision pairs.)" << endl;
  }
  
  NearCollisionTest nc_test(r, q, e, unif, lsh_dim, lsh_length, iteration_multiple);

  auto near_collision_num = 1000;
  
  auto L = nc_test.gen_instance(near_collision_num, N);
  auto sol = nc_test.lsh_based_search(L);
  double avg_sol = sol.size() / near_collision_num;
  
  for (size_t i = 1; i < repetition; i++) {
    auto L = nc_test.gen_instance(near_collision_num, N);
  auto sol = nc_test.lsh_based_search(L);

  avg_sol = (avg_sol * i + (double) sol.size()) / (i+1);
    if ((i+1) % print_term == 0) {
      cout << setw(3) << "- " << i+1 << " executions: E(# Found pairs) = " << avg_sol / near_collision_num * 100 << " %" << endl;
    }
  }
  
  return 0;
}
