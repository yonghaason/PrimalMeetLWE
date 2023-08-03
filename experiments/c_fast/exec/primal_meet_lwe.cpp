#include "functions.h"
#include "cmdline.h"
#include <iostream>
#include <iomanip>
#include <random>

using namespace std;

int main(int argc, char* argv[])
{
  uint32_t d, h, w;
  uint64_t q;
  double l;

  cmdline::parser parser;

  parser.add<uint32_t>("d", 'd', "Dimension (r^(i-1) in paper)", true, 20);
  parser.add<uint32_t>("h", 'h', "Secret weight (w^(i-1) in paper)", true, 0);
  parser.add<uint32_t>("w", 'w', "Guess weight (w^i in paper)", true, 0);
  parser.add<double>("l", 'l', "Constraint bound (ell^(i) in paper)", true, 0);
  parser.add<uint64_t>("q", 'q', "Modulus param that determines the domain", false, 2048);

  parser.parse_check(argc, argv);

  d = parser.get<uint32_t>("d");
  h = parser.get<uint32_t>("h");
  w = parser.get<uint32_t>("w");
  q = parser.get<uint64_t>("q");
  l = parser.get<double>("l");

  auto m = d;
  
  UniformTest unif_test(d, 1, q, h, m, true);
  matrix B = unif_test.B;
  auto GSnorm = unif_test.GSnorm;
  matrix M; secret s; vector<double> e;
  unif_test.gen_noisy_instance(M, s, e);
  auto Mtrans = transpose(M);

  /**
   * 1. Top level list construction
   * 2. Near-collision
   * 3. Weight-Check 
  */

  auto top_pair = unif_test.sparse_secret_list(Mtrans, d, h/4);
  map<vector<double>, secret> mapping;
  vector<domain> top_Ms;
  for (auto &pair: top_pair)
  {
    mapping[pair.second] = pair.first;
  }
  NearCollisionTest nc_test(r, q, e, unif, lsh_dim, lsh_length, iteration_multiple);
  auto sol = nc_test.lsh_based_search(top_Ms);

  // 미완

  return 0;
}