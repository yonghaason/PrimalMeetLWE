#include "functions.h"
#include "cmdline.h"
#include "torus_lsh.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>

void enumerate(vector<double> c, vector<int64_t> coeffs, int k, double& normBound, matrix& B, 
               set<pair<vector<int64_t>, vector<double>>>& list)
{
  auto n = B.size();
  auto coeff_range_min = ceil((-normBound-c[k-1])/B[k-1][k-1] - 0.5);
  auto coeff_range_max = floor((normBound-c[k-1])/B[k-1][k-1] + 0.5);
  
  for (int coeff = coeff_range_min; coeff < coeff_range_max+1; coeff++) {
    vector<double> next_c(c);
    vector<int64_t> next_coeffs(coeffs);
    for (size_t i = 0; i < k; i++) {
      next_c[i] += coeff * B[i][k-1];
    }
    next_coeffs[k-1] = coeff;
    
    if (k > n/2 + 1) {
      enumerate(next_c, next_coeffs, k-1, normBound, B, list);
    }
    else { // if k == n/2
      // Size reduce for 0 ~ n/2-1 -th coords
      for (int i = n/2 - 1; i >= 0; i--) 
      {
        auto coeff = round(next_c[i] / B[i][i]);
        next_coeffs[i] = coeff;
        for (int j = 0; j <= i; j++)
        {
          next_c[j] -= coeff * B[j][i];
        }
      }
      list.insert(std::make_pair(next_coeffs, next_c));
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
    if (i < n/2) {
      profile[n-1-i] = lastGS * pow(rhf, 2*i);
    }
    else {
      profile[n-1-i] = profile[n-i];
    }
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
      B[i][j] = coord_sampler[i](gen);
    }
  }

  // std::cout << "B: " << endl;
  // for (size_t i = 0; i < n; i++) {
  //   print(B[i]);
  // }
  // std::cout << endl;
  
  double normBound = pow(det, 1.0/n) * 0.5;
  // double normBound = profile[n/2] * sqrt(1.0/3);
  std::cout << "norm bound = det(B)^{1/n} * 2/3: " << normBound << std::endl;
  // double oldNormBound = profile[n/2 - 1] * 2/3;
  // std::cout << "(old norm bound = GSnorm[n/2-1] * 2/3: " << oldNormBound << ")" << std::endl;

  double expectedListSize = 1.0;
  for (size_t i = n/2; i < n; i++) {
    expectedListSize *= ceil(1 + 2*normBound / profile[i]);
  }



  tik = chrono::system_clock::now();

  set<pair<vector<int64_t>, vector<double>>> P;
  // Preprocess P = {x in L(B) : ||proj(x, n/2)||_oo <= normBound}
  vector<double> x(n);
  vector<int64_t> coeffs(n);
  enumerate(x, coeffs, n, normBound, B, P);
  
  tok = chrono::system_clock::now();

  std::cout << "Preprocessed set size: " << P.size() 
            << " / " << chrono::duration_cast<chrono::milliseconds>(tok - tik).count() 
            << " ms" << std::endl;

  std::cout << "Expected set size: " << (size_t) expectedListSize << std::endl;

  // for (auto x: P) {
  //   print(x.first);
  //   print(x.second);
  //   cin.ignore();
  // }

  set<vector<double>> PP;
  for (auto pair: P) {
    PP.insert(pair.second);
  }

  // prepare LSH table;
  vector<double> lshDomain(profile.begin(), profile.begin() + n/2);
  double block_length = 4*normBound;
  TorusLSH torusLSH(block_length, lshDomain);
  size_t repeatLsh = 10;
  torusLSH.setPartitions(repeatLsh);
  auto hashTables = torusLSH.computeLshTable(PP);
  
  tik = chrono::system_clock::now();

  std::cout << "LSH Table Done (with block-length " << block_length
            << ") / " << chrono::duration_cast<chrono::milliseconds>(tik - tok).count() 
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
    for (auto prep: PP) {
      auto dist = subvec(testVecs[iter], prep);
      if (inf_norm(dist, 0, n) <= normBound) {
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