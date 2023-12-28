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
  parser.add<bool>("top-opti", '\0', "top-level optimization toggle", false, true);

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
    ell[i] = ell[i-1];
    b[i] = 2*ell[i-1];
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
    if (top_opti && i == t) continue;
    
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
      throw invalid_argument("Cannot set projection dimension r" + to_string(i));
  }

  vector<uint64_t> R_ncf(t+1);
  vector<double> p_nc(t+1);
  vector<double> p_lsh(t+1);
  vector<double> p_col(t+1);
  for (int i = t; i >= 1; i--)
  {  
    domain dom(r[i-1]);
    for (int k = 0; k < r[i-1] - r[i]; k++) 
      {dom[k] = q[m - r[i-1] + k];}
    for (int k = r[i-1] - r[i]; k < r[i-1]; k++) 
      {dom[k] = ell[i];}

    p_nc[i] = 1;
    for (int k = 0; k < r[i-1] - r[i]; k++)
    {
      if (i == 1) p_nc[i] *= prob_admissible_gaussian(stddev, q[m - r[i-1] + k]);
      else p_nc[i] *= prob_admissible_uniform(ell[i-1], q[m - r[i-1] + k]);
    }
    p_lsh[i] = 1;
    p_col[i] = 1;

    vector<int> n_tmp(r[i-1]);
    vector<double> b_tmp(r[i-1]);
    for (size_t k = 0; k < r[i-1]; k++) {
      n_tmp[k] = max(1, (int) floor(dom[k] / b[i]));
      b_tmp[k] = (double) dom[k] / n_tmp[k];
      if (n_tmp[k] > 1) 
      {
        if (i != 1) {p_lsh[i] *= prob_admissible_uniform(ell[i-1], b_tmp[k]);}
        else {p_lsh[i] *= prob_admissible_gaussian(stddev, b_tmp[k]);}
      }
      p_col[i] *= 1.0 / n_tmp[k];
    }

    R_ncf[i] = ceil(C/p_lsh[i]);
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
  size_t top_level_size = top_level_full_size;
  if (!top_opti)
  {
    top_level_size = top_level_full_size;
  }
  else
  {
    top_level_size = sqrt(3) * (double) top_level_full_size / sqrt(R[t]);
    expected_list_size[t] = top_level_size;      
  }

  
  std::cout << "---------- Parms ----------" << endl;
  std::cout << "Gram-Schmidt norm = " << q[0] << " -> " << q[m-1] << " (assume GSA with rhf =" << rhf << ")" << endl;
  std::cout << "      w, r, ell, b_lsh, R_ncf, p_nc, R" << endl;
  std::cout << "lv 0: " << w[0] << ", " << r[0] << ", " << ell[0] << ", -, " << R_ncf[1] << ", " << p_nc[1] << ", " << R[1] << endl;
  for (size_t i = 1; i < t; i++) {
    std::cout << "lv " << i << ": " << w[i] << ", " << r[i] << ", " << ell[i] << ", " << b[i] << ", " << R_ncf[i+1] <<  ", " << p_nc[i+1] << ", " << R[i+1] << endl;
  }
  if (top_opti) {
    std::cout << "lv " << t << ": " << w[t] << ", -, -, " << b[t] << ", -, -, -" << endl;  
  }
  else {
    std::cout << "lv " << t << ": " << w[t] << ", " << r[t] << ", " << ell[t] << ", " << b[t] << ", -, -, -" << endl; 
  }
  std::cout << endl;
  
  /**
   * p_sp[0] = Over the last r[1] dimension, norm bound ell[1], Gaussian_stddev
   * p_sp[i] = Over the last r[i+1] dimension, norm bound ell[i+1],  Unif[-ell[i], ell[i]]
   * p_sp[t-1] = (1 - e^-3)**2^{t-1}
   * 
   * p_ncf[0] = Over the last r[0](=m) dimension, lsh block-size b_lsh[1], Gaussian_stddev
   * p_ncf[i] = Over the last r[i] dimension, norm bound ell[i+1], Unif[-ell[i], ell[i]]
   * 
   * p_nc[0] = From 0 to r[0] - r[1] dimension, admissible for box length q[i], Gaussian error
   * p_nc[i] = From r[i-1] - r[i] to r[i-1] - r[i+1] dimension, admissible for box length q[i], Unif[-ell[i], ell[i]] error
   * 
   * p_suc = Prod_{i = 0 to t-1} (p_sp[i] * p_ncf[i] * p_nc[i])^{2^i} 
  */
  vector<double> p_r_gauss{0.864, 0.857, 0.854, 0.849, 0.845, 0.842, 0.836, 0.833, 0.830, 0.826};
  vector<double> p_r_unif{0.815, 0.774, 0.742, 0.706, 0.675, 0.643, 0.617, 0.592, 0.57, 0.549};

  double p_suc_theory = p_r_gauss[r[1]/5] * p_r_gauss[r[0]/5] * p_nc[1];
  for (int i = 1; i < t; i++)
  {
    p_suc_theory *= pow(p_r_unif[r[i+1]/5] * p_r_unif[r[i+1]/5] * p_nc[i+1], (1 << i));
  }

  cout << "* Success probability lower bound (by theoretic analysis): " << p_suc_theory << endl;
  cout << endl;   
  
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
      std::cout << "- Step 1. Construct L" << i << " ⊂ {s, [Ms]_{B, " << r[i-1] << "}: HW(s) = " << w[i] << "} of size " << top_level_size << " (full size = " << top_level_full_size << ")" << endl;;
    }
    else 
    {
      std::cout << "- Step 1. Construct L" << i << " = {s, [Ms]_{B, " << r[i-1] << "}: HW(s) = " << w[i] << " and [Ms]_{B, " << r[i] << "} in [" << -ell[i] << ", " << ell[i] << "]^" << r[i] << "}" << endl;
    }
    std::cout << "- Step 2. Recover S" << i-1 << " = {s: HW(s) = " << w[i-1] 
    << " and [Ms]_{B, " << r[i-1] << "} in [" << -ell[i-1] << ", " << ell[i-1] 
    << "]^" << r[i-1] << "} (NCF with block-length " << b[i] << ")" << endl;
    // std::cout << "  - One LSH succeeds with " << p_lsh[i] << " probability -> " << R_ncf[i] << " torus-LSH iterations" << endl;
    // std::cout << "  - (False) near-collision prob p_bad = " << p_col[i] << endl;
  }
  cout << endl;

  
  /////////////////////////////////////////////////////////////////////////////////////////

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
      candidate_S = NCF_lsh(dom, L, ell[i-1], w[i-1], b[i], R_ncf[i], collision_nums);
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
  std::cout << "|S_i|: real, theory" << endl;
  for (size_t i = 1; i <= t; i++)
    std::cout << "lv " << i << ": " << (double) list_sizes[i] / repeat << ", " 
        << expected_list_size[i] << endl;

  return 0;
}