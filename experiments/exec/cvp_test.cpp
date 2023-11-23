#include "functions.h"
#include "cmdline.h"
#include "torus_lsh.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>

void enumerate(vector<double> c, int k, double& normBound, matrix& B, 
               set<vector<double>>& list)
{
  auto n = B.size();
  // auto coeff_range_min = ceil((-normBound-c[k-1])/B[k-1][k-1]);
  // auto coeff_range_max = floor((normBound-c[k-1])/B[k-1][k-1]);
  auto coeff_range_min = ceil((-normBound-c[k-1])/B[k-1][k-1] - 0.5);
  auto coeff_range_max = floor((normBound-c[k-1])/B[k-1][k-1] + 0.5);
  for (int coeff = coeff_range_min; coeff < coeff_range_max+1; coeff++) {
    vector<double> next_c(c);
    for (size_t i = 0; i < k; i++) {
      next_c[i] += coeff * B[i][k];
    }
    if (k > n/2) {
      enumerate(next_c, k-1, normBound, B, list);
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
  chrono::_V2::system_clock::time_point tik, tok;

  cmdline::parser parser;

  parser.add<uint32_t>("n", 'n', "lattice dimension", true, 0);
  parser.add<double>("lastGS", '\0', "modulus param", false, 100);
  parser.add<double>("rhf", '\0', "root-hermite-factor", false, 1.2);

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

  std::cout << "Basis GSnorm: " << profile[0] << " -> " << profile[n-1] << " (rhf = " << rhf << ")" << std::endl;
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
  
  double normBound = pow(det, 1.0/n) * 2/3;
  std::cout << "norm bound = det(B)^{1/n} * 2/3: " << normBound << std::endl;
  // double oldNormBound = profile[n/2 - 1] * 2/3;
  // std::cout << "(old norm bound = GSnorm[n/2-1] * 2/3: " << oldNormBound << ")" << std::endl;

  tik = chrono::system_clock::now();

  set<vector<double>> P;
  // Preprocess P = {x in L(B) : ||proj(x, n/2)||_oo <= normBound}
  vector<double> c(n);
  enumerate(c, n, normBound, B, P);

  tok = chrono::system_clock::now();

  std::cout << "Preprocessed set size: " << P.size() 
            << " / " << chrono::duration_cast<chrono::milliseconds>(tok - tik).count() 
            << " ms" << std::endl;

  // prepare LSH table;
  vector<double> lshDomain(profile.begin(), profile.begin() + n/2);
  TorusLSH torusLSH(4*normBound, lshDomain);
  size_t repeatLsh = 10;
  torusLSH.setPartitions(repeatLsh);
  auto hashTables = torusLSH.computeLshTable(P);
  
  tik = chrono::system_clock::now();

  std::cout << "LSH Table Done" 
            << " / " << chrono::duration_cast<chrono::milliseconds>(tik - tok).count() 
            << " ms" << std::endl;

  vector<vector<double>> testVecs;
  for (size_t iter = 0; iter < repeat; iter++) {
    vector<double> test(n);
    for (size_t i = 0; i < n; i++) {
      test[i] = coord_sampler[i](gen);
    }
    testVecs.push_back(test);
  }

  tok = chrono::system_clock::now();
  
  size_t lshSuccessCount = 0;
  for (size_t iter = 0; iter < repeat; iter++) {
    // LSH Search    
    vector<double> cvp;
    if (torusLSH.searchFromLshTable(testVecs[iter], hashTables, normBound, cvp)) {
      lshSuccessCount++;
    };
  }

  tik = chrono::system_clock::now();

  std::cout << "LSH Success : " << lshSuccessCount 
            << " / " << chrono::duration_cast<chrono::milliseconds>(tik - tok).count() 
            << " ms" << std::endl;

  size_t exactSuccessCount = 0;
  for (size_t iter = 0; iter < repeat; iter++) {
    // Exact Search
    bool flag = false;
    for (auto prep: P) {
      auto dist = subvec(testVecs[iter], prep);
      if (inf_norm(dist) <= normBound) {
        // std::cout << "t: ";
        // print(test);
        // std::cout << "l: ";
        // print(prep);
        // std::cout << "[t]_B: ";
        // print(dist);
        // std::cout << "mod B norm: " << inf_norm(dist) << std::endl;
        flag = true;
        // std::cout << std::endl;
        break;
      }
    }
    if (flag) exactSuccessCount++;
  }

  tok = chrono::system_clock::now();

  std::cout << "Exact Success : " << exactSuccessCount 
            << " / " << chrono::duration_cast<chrono::milliseconds>(tok - tik).count() 
            << " ms" << std::endl;

  std::cout << "Total exps : " << repeat << std::endl;
  
}