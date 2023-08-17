#include "functions.h"
#include "cmdline.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>

using namespace std;

int main(int argc, char* argv[])
{ 
  cmdline::parser parser;

  parser.add<uint32_t>("m", 'm', "matrix dimension", true, 0);
  parser.add<uint32_t>("d", 'd', "secret dimension", true, 0);
  parser.add<uint32_t>("h", 'h', "secret weight", true, 0);
  parser.add<double>("stddev", '\0', "error stddev", false, 1);
  parser.add<double>("q", 'q', "modulus param", false, 2048);

  parser.add<uint32_t>("t", 't', "level t", false, 1);
  parser.add<uint32_t>("C", 'C', "lsh iteration multiple", false, 3);
  parser.add<double>("rhf", '\0', "root-hermite-factor", false, 1.05);
  parser.add<bool>("top-opti", '\0', "top-level optimization toggle", false, false);

  parser.add<uint64_t>("repeat", '\0', "# of experiments", false, 1);

  parser.parse_check(argc, argv);

  // LWE parameters
  auto m = parser.get<uint32_t>("m");
  auto d = parser.get<uint32_t>("d");
  auto h = parser.get<uint32_t>("h");
  auto stddev = parser.get<double>("stddev");
  auto q0 = parser.get<double>("q");

  // Attack parameters  
  auto t = parser.get<uint32_t>("t");
  auto C = parser.get<uint32_t>("C");
  auto rhf = parser.get<double>("rhf");
  auto top_opti = parser.get<bool>("top-opti");

  auto repeat = parser.get<uint64_t>("repeat");

  vector<double> q(m);
  for (int i = 0; i < m; i++) {
    q[m-i-1] = q0 * pow(rhf, i);
  }

  vector<double> ell(t+1);
  vector<double> b(t+1);
  ell[0] = 6*stddev;
  b[0] = 0; // dummy; never used
  for (size_t i = 1; i <= t; i++) 
  {
    ell[i] = 3*ell[i-1];
    b[i] = 1.5*ell[i-1];
    // b[i] = ell[i];
  }

  vector<uint32_t> w(t+1);
  w[0] = h;
  for (size_t i = 1; i <= t; i++)
  {
    std::cout << "Insert w[" << i << "]: "; 
    cin >> w[i];
    auto eps = w[i] - w[i-1]/2;
    if (eps < 0)
      throw invalid_argument("w[i] < w[i-1]/2, fail");
    if ((i != t) && (w[i] % 2 != 0))
      throw invalid_argument("Intermediate level w[i] should be even");
    if (eps > d - w[i-1])
      throw invalid_argument("w[i] > d - w[i-1]/2, fail");
  }

  vector<uint64_t> R(t+1);
  R[0] = 0; // dummy: never used
  for (size_t i = 1; i <= t; i++)
  {
    R[i] = ambiguity(d, w[i-1], w[i]); 
  }

  vector<int> r(t+1);
  vector<double> p_rep(t+1);
  vector<double> vol_ratio(t+1);
  r[0] = m;
  for (size_t i = 1; i <= t; i++)
  {
    r[i] = 0;
    p_rep[i] = 1;
    double cur_prob = 1;
    vol_ratio[i] = 1;
    while (true)
    {
      double cur_len = min(q[m-r[i]-1], 2*ell[i]);
      double cur_ratio = cur_len / q[m-r[i]-1];
      if (i == 1)
        cur_prob = cur_ratio * prob_admissible_gaussian(stddev, cur_len);
      else
        cur_prob = cur_ratio * prob_admissible_uniform(ell[i-1], cur_len);
      if (R[i] * p_rep[i] * cur_prob > 2*C) 
      {
        p_rep[i] *= cur_prob;
        vol_ratio[i] *= cur_ratio;
        r[i] += 1;
      }
      else break;
      if (r[i] == r[i-1]) break;
    }
    if (r[i] == 0) 
      throw invalid_argument("Cannot set projection dimension r[i]");
  }

  vector<uint64_t> R_lsh(t+1);
  vector<double> p_goods(t+1);
  vector<double> p_bads(t+1);
  for (int i = t; i >= 1; i--)
  {  
    domain dom(r[i-1]);
    for (int k = 0; k < r[i-1] - r[i]; k++) 
      {dom[k] = q[m - r[i-1] + k];}
    for (int k = r[i-1] - r[i]; k < r[i-1]; k++) 
      {dom[k] = ell[i];}

    p_goods[i] = 1;
    p_bads[i] = 1;

    vector<int> n_tmp(r[i-1]);
    vector<double> b_tmp(r[i-1]);
    for (size_t k = 0; k < r[i-1]; k++) {
      n_tmp[k] = max(1, (int) floor(dom[k] / b[i]));
      b_tmp[k] = (double) dom[k] / n_tmp[k];
      if (n_tmp[k] > 1) 
      {
        if (i != 1) {p_goods[i] *= prob_admissible_uniform(ell[i-1], b_tmp[k]);}
        else {p_goods[i] *= prob_admissible_gaussian(stddev, b_tmp[k]);}
      }
      p_bads[i] *= 1.0 / n_tmp[k];
    }

    R_lsh[i] = C/p_goods[i];
  }

  vector<double> expected_list_size(t+1);
  for (size_t i = 1; i <= t; i++) 
  {
    uint64_t list_size = binom(d, w[i]) * (1ull << w[i]);
    double vol_ratio = 1;
    for (size_t j = 1; j <= r[i]; j++) 
    {
      double cur_len = min(q[m-j], 2*ell[i]);
      vol_ratio *= q[m-j] / cur_len;
    }
    expected_list_size[i] = (double) list_size / vol_ratio;
  }

  
  size_t top_level_full_size = binom(d, w[t]) * (1ull << w[t]);
  size_t top_level_size;
  if (!top_opti)
  {
    top_level_size = top_level_full_size;
  }
  else
  {
    top_level_size = C * (double) top_level_full_size / sqrt(R[t]);
    expected_list_size[t] = top_level_size;
  }

  std::cout << "---------- Parms ----------" << endl;
  std::cout << "Coordinate length (q) = " << q[0] << " -> " << q[m-1] << " (assume GSA)" << endl;
  std::cout << "   w, r, ell, b, R, Rp/2, R_lsh" << endl;
  std::cout << 0 << ": " << w[0] << ", " << r[0] << ", " << ell[0] << ", -, -, -, -" << endl;
  for (size_t i = 1; i <= t; i++) {
    std::cout << i << ": " << w[i] << ", " << r[i] << ", " << ell[i] << ", " << b[i] << ", " << R[i] << ", " << R[i] * p_rep[i] / 2 << ", " << R_lsh[i] << endl;
  }
  std::cout << endl;
  
  //////////////////////////////////////////////////////

  int success_count = 0;
  int unwanted_count = 0;
  vector<double> list_sizes(t+1);
  vector<double> collision_numbers(t+1);

  std::cout << "------ Attack Procedure Overview ------ " << endl;
  for (int i = t; i >= 1; i--)
  {
    std::cout << "Level " << i << endl;
    if ((i == t) && top_opti)
    {
      std::cout << "- Step 1. Construct L" << i << " = {s, [Ms]_{B, " << r[i-1] << "}: HW(s) = " << w[i] << "} of size " << top_level_size << endl;
    }
    else 
    {
      std::cout << "- Step 1. Construct L" << i << " = {s, [Ms]_{B, " << r[i-1] << "}: HW(s) = " << w[i] << " and [Ms]_{B, " << r[i] << "} in [" << -ell[i] << ", " << ell[i] << "]^" << r[i] << "}" << endl;
    }
    std::cout << "- Step 2. Recover S" << i-1 << " = {s: HW(s) = " << w[i-1] 
    << " and [Ms]_{B, " << r[i-1] << "} in [" << -ell[i-1] << ", " << ell[i-1] 
    << "]^" << r[i-1] << "} (NCF with block-length " << b[i] << ")" << endl;
    std::cout << "  - One LSH succeeds with " << p_goods[i] << " probability -> " << R_lsh[i] << " torus-LSH iterations" << endl;
    std::cout << "  - (False) near-collision prob p_bad = " << p_bads[i] << endl;
  }
  cout << endl;

  auto start = chrono::steady_clock::now();

  for (size_t iter = 0; iter < repeat; iter++) 
  {
    std::cout << "... Runnning " << iter + 1 << "-th, until now " 
      << success_count << " successes" << '\r' << std::flush;
    
    matrix M; matrix B; secret s; vector<double> e;
    gen_noisy_instance(
      m, d, h, stddev, q, 
      M, B, s, e);
    auto Mtrans = transpose(M);

    // Sanity check
    // auto MsB = matmul(Mtrans, s);
    // MsB = babaiNP(MsB, B);
    // print(MsB); print(e);

    set<secret> candidate_S;
    
    if (top_opti)
    {
      while (candidate_S.size() < top_level_size) 
      {
        set<secret> basic_secret = enumerate_secrets(w[t], w[t]);
        vector<secret> basic_secret_vec(basic_secret.begin(), basic_secret.end());
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<int64_t> randindex(0, basic_secret_vec.size() - 1);
        auto index = randindex(gen);
        secret tmp(basic_secret_vec[index]);
        tmp.resize(d);
        random_shuffle(tmp.begin(), tmp.end());
        candidate_S.insert(tmp);
      }
    }
    else 
    {
      candidate_S = enumerate_secrets(d, w[t]);
    }
    
    for (int i = t; i >= 1; i--) 
    {
      list L;
      if ((i == t) && !top_opti)
      {
        for (auto s_cand: candidate_S)
        {
          auto Ms = babaiNP(matmul(Mtrans, s_cand), B, r[t-1]);
          if (inf_norm(Ms, r[t-1] - r[t]) < ell[t]) 
            L.push_back(make_pair(s_cand, Ms));
        }
      }
      else 
      {
        for (auto s_cand: candidate_S)
        {
          auto Ms = babaiNP(matmul(Mtrans, s_cand), B, r[i-1]);
          L.push_back(make_pair(s_cand, Ms));
        }
      }

      list_sizes[i] += L.size();
      
      domain dom(r[i-1]);
      for (int k = 0; k < r[i-1] - r[i]; k++) 
        {dom[k] = q[m - r[i-1] + k];}
      for (int k = r[i-1] - r[i]; k < r[i-1]; k++) 
        {dom[k] = ell[i];}
          
      size_t collision_nums = 0;
      candidate_S = NCF_lsh(dom, L, ell[i-1], w[i-1], b[i], R_lsh[i], collision_nums);
      collision_numbers[i] += collision_nums;
    }

    if (candidate_S.size() != 0)
    {
      if (find(candidate_S.begin(), candidate_S.end(), s) != candidate_S.end())
        success_count++;
      if (candidate_S.size() > 2) 
        unwanted_count++;
    }
  }

  auto end = chrono::steady_clock::now();

  std::cout << "----------- Takes " 
      << chrono::duration_cast<chrono::milliseconds>(end - start).count()
      << " ms for " << repeat << " runs -----------------" << endl;
  std::cout << "Success prob: " << (double) success_count / repeat
        << " (# multiple sol: " << (double) unwanted_count /repeat << ")" << endl;
  std::cout << "    L: real, expected" << endl;
  for (size_t i = 1; i <= t; i++)
    std::cout << i << ": " << (double) list_sizes[i] / repeat << ", " 
        << expected_list_size[i] << endl;
  std::cout << "    #cols" << endl;
  for (size_t i = 1; i <= t; i++) 
    std::cout << i << ": " << (double) collision_numbers[i] / repeat << endl;
  std::cout << endl;

  return 0;
}