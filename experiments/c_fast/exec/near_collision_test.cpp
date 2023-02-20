#include "functions.h"
#include "cmdline.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <map>

// Example
// ./nc_test -m 25 -q 4096 -e 3 -b 6 -r 10 -u 1

using namespace std;

int main(int argc, char* argv[]) {

  chrono::system_clock::time_point start, end;

  uint32_t m; 
  uint64_t q;
  double e;
  double lsh_length;
  size_t lsh_dim;
  uint64_t iteration_multiple;
  uint64_t N;
  bool unif;

  cmdline::parser parser;

  parser.add<uint32_t>("m", 'm', "Dimension of list vectors", true, 20);
  parser.add<uint64_t>("q", 'q', "Domain length", true, 1024);
  parser.add<double>("stddev", 'e', "error stddev (gaussian) or bound (unif)", true, 3.2);
  parser.add<double>("lsh_length", 'b', "LSH box length", true, 0);
  parser.add<size_t>("lsh_dim", 'r', "LSH dimension", true, 10);
  parser.add<uint64_t>("iter_mult", '\0', "Iteration multiple", false, 3);
  parser.add<bool>("unif", 'u', "Error unif?", true, true);
  parser.add<uint64_t>("list_size", 'N', "List size", false, 0);
  
  parser.parse_check(argc, argv);

  m = parser.get<uint32_t>("m");
  q = parser.get<uint64_t>("q");
  e = parser.get<double>("stddev");
  lsh_length = parser.get<double>("lsh_length");
  lsh_dim = parser.get<size_t>("lsh_dim");
  N = parser.get<uint64_t>("list_size");
  iteration_multiple = parser.get<uint64_t>("iter_mult");
  unif = parser.get<bool>("unif");
  
  if (unif && (lsh_length < e)) {
    cout << "!!! Some error can exceed LSH length !!!\n"
         << "(= LSH never finds some near-collision pairs.)" << endl;
  }
  if (!unif && (lsh_length < 3*e)) { // some arbitrary choice
    cout << "!!! Some error can exceed LSH length !!!\n"
         << "(= LSH never finds some near-collision pairs.)" << endl;
  }
  
  NearCollisionTest nc_test(m, q, e, unif);

  auto near_collision_num = 1000;
  auto L = nc_test.gen_instance(near_collision_num, N);
  auto sol = nc_test.lsh_based_search(L, lsh_length, lsh_dim, iteration_multiple);
  cout << endl;
  cout << "Found " << sol.size() << " good pairs " 
      << "(Expected = " << near_collision_num << ")" << endl;
  
  return 0;
}
