#include "functions.h"
#include "cmdline.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>

void prepare_LSH(domain D, set<vector<double>>& list, double block_length,
                 vector<vector<double>>& lsh_partitions)
{
  random_device rd;
  mt19937 gen(rd());

  for (size_t iter = 0; iter < R_lsh; iter++)
  {
    // Pick torus-LSH by starting points
    vector<double> lsh_partition(r);
    for (size_t i = 0; i < r; i++) {
      uniform_real_distribution<double> lsh_random_starts(0, b+xxcemsskk[i]);
      lsh_partition[i] = lsh_random_starts(gen);
    }

    // Fill LSH table
    map<vector<int32_t>, list> hash_table;
    for (auto s_Ms: L) {
      auto Ms = s_Ms.second;
      std::vector<int32_t> address(r);
      for (size_t i = 0; i < r; i++) {
        int32_t a = floor((Ms[i] + dom[i] - lsh_partition[i]) / b[i]);
        address[i] = a % n[i];
      }
      if (hash_table.count(address) != 0) {
        hash_table[address].push_back(s_Ms);
      }
      else {
        list s_Ms_copy{s_Ms};
        hash_table[address] = s_Ms_copy;
      }
    }
  }
}

void enumerate(vector<double> c, int k, double& norm_bound, matrix& B, 
               set<vector<double>>& list)
{
  auto n = B.size();
  // auto coeff_range_min = ceil((-norm_bound-c[k-1])/B[k-1][k-1]);
  // auto coeff_range_max = floor((norm_bound-c[k-1])/B[k-1][k-1]);
  auto coeff_range_min = ceil((-norm_bound-c[k-1])/B[k-1][k-1] - 0.5);
  auto coeff_range_max = floor((norm_bound-c[k-1])/B[k-1][k-1] + 0.5);
  for (int coeff = coeff_range_min; coeff < coeff_range_max+1; coeff++) {
    vector<double> next_c(c);
    for (size_t i = 0; i < k; i++) {
      next_c[i] += coeff * B[i][k];
    }
    if (k > n/2) {
      enumerate(next_c, k-1, norm_bound, B, list);
    }
    else { // if k == n/2
      // Size reduce for 0 ~ n/2 -th coords
      for (int i = n/2; i >= 0; i--) 
      {
        for (int j = 0; j < i; j++)
        {
          next_c[j] -= round(next_c[i] / B[i][i]) * B[j][i];
        }
      }
      list.insert(next_c);
    }
  }
}

int main(int argc, char* argv[]) 
{
  chrono::system_clock::time_point start, end;

  cmdline::parser parser;

  parser.add<uint32_t>("n", 'n', "lattice dimension", true, 0);
  parser.add<double>("lastGS", '\0', "modulus param", false, 100);
  parser.add<double>("rhf", '\0', "root-hermite-factor", false, 1.05);

  parser.add<uint64_t>("repeat", '\0', "# of experiments", false, 1000);
  
  parser.parse_check(argc, argv);

  auto n = parser.get<uint32_t>("n");
  auto lastGS = parser.get<double>("lastGS");
  auto rhf = parser.get<double>("rhf");
  auto repeat = parser.get<uint64_t>("repeat");

  // Simulate the basis B by its Gram-Schmidt norm using GSA.
  vector<double> profile(n);
  double det = 1.0;
  for (int i = 0; i < n; i++) {
    profile[n-1-i] = lastGS * pow(rhf, i) * pow(rhf, i);
    det *= profile[n-1-i];
  }

  cout << "Basis GSnorm: " << profile[0] << " -> " << profile[n-1] << endl;
  cout << "determinant: " << det << endl;

  random_device rd;
  mt19937 gen(rd());
  vector<uniform_real_distribution<double>> coord_sampler(n);

  matrix B;

  B.resize(n);
  for (size_t i = 0; i < n; i++) {
    B[i].resize(n);
  }
  for (size_t j = 0; j < n; j++) {
    B[j][j] = profile[j];
    coord_sampler[j] = uniform_real_distribution<double>(-profile[j]/2, profile[j]/2);
    for (size_t i = 0; i < j; i++) {
      B[i][j] = coord_sampler[j](gen);
    }
  }

  double old_norm_bound = profile[n/2 - 1] * 2/3;
  double norm_bound = pow(det, 1.0/n) * 2/3;
  cout << "norm_bound = det(B)^{1/n} * 2/3: " << norm_bound << endl;
  cout << "(old norm_bound = GSnorm[n/2-1] * 2/3: " << old_norm_bound << ")" << endl;
  set<vector<double>> P;
  // Preprocess P = {x in L(B) : ||proj(x, n/2)||_oo <= norm_bound}
  vector<double> c(n);
  enumerate(c, n, norm_bound, B, P);

  cout << "Preprocessed set size: " << P.size() << endl;

  // prepare_LSH();
  size_t succ = 0;
  for (size_t iter = 0; iter < repeat; iter++) {
    vector<double> test(n);
    for (size_t i = 0; i < n; i++) {
      test[i] = coord_sampler[i](gen);
    }
    bool flag = false;
    // near_collision_search(); 
    for (auto prep: P) {
      auto dist = subvec(test, prep);
      if (inf_norm(dist) <= norm_bound) {
        // cout << "t: ";
        // print(test);
        // cout << "l: ";
        // print(prep);
        // cout << "[t]_B: ";
        // print(dist);
        // cout << "mod B norm: " << inf_norm(dist) << endl;
        flag = true;
        // cout << endl;
        break;
      }
    }
    if (flag) succ++;
  }
  cout << "Success : " << succ << " / " << repeat << endl;
  
}