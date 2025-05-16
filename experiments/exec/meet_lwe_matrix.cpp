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

  parser.add<uint32_t>("t", 't', "level t", false, 2);
  parser.add<uint32_t>("C", 'C', "lsh iteration multiple", false, 2);
  parser.add<double>("rhf", '\0', "root-hermite-factor", false, 1.05);
  
  parser.add<uint32_t>("w1", '\0', "lv 1 weight", false, -1);
  parser.add<uint32_t>("w2", '\0', "lv 2 weight", false, -1);

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
  
  auto repeat = parser.get<uint64_t>("repeat");

  vector<double> q(m);
  for (int i = 0; i < m; i++) {
    q[m-i-1] = q0 * pow(rhf, 2*i);
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
  if (parser.get<uint32_t>("w1") != -1 && parser.get<uint32_t>("w2") != -1) {
    w[1] = parser.get<uint32_t>("w1");
    w[2] = parser.get<uint32_t>("w2");
  }
  else {
    for (size_t i = 1; i <= t; i++)
    {
      std::cout << "Insert w[" << i << "]: "; 
      cin >> w[i];
    }
  }
  
  for (size_t i = 1; i <= t; i++) {
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
    if (i == t) continue;
    
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

  vector<uint64_t> R_lsh(t+1);
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

    R_lsh[i] = ceil(C/p_lsh[i]);
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
  
  top_level_size = sqrt(3) * (double) top_level_full_size / sqrt(R[t]);
  expected_list_size[t] = top_level_size;      
  
  
  std::cout << "---------- Parms ----------" << endl;
  std::cout << "Gram-Schmidt norm = " << q[0] << " -> " << q[m-1] << " (assume GSA with rhf =" << rhf << ")" << endl;
  std::cout << "      " 
    << setw(6) << "w," << setw(6) << "r," << setw(6) << "ell," << setw(7) << "b_lsh," << setw(8) << "R_lsh,"  
    << setw(11) << "C/p_lsh," << setw(11) << "p_nc," << setw(11) << "C/p_rep," << setw(10) << "R," << setw(9) << "size" << endl;
  std::cout << "lv 0: " 
    << setw(5) << w[0] << "," << setw(5) << r[0] << "," << setw(5) << ell[0] << "," << setw(7) << "-," << setw(7) << R_lsh[1] << ","
    << setw(11) << "," << setw(10) << p_nc[1] << "," << setw(11) << "-," << setw(9) << R[1] << ", " << setw(8) << "-" << endl;
  for (size_t i = 1; i < t; i++) {
    std::cout << "lv " << i << ": " 
    << setw(5) << w[i] << "," << setw(5) <<  r[i] << "," << setw(5) << ell[i] << "," << setw(6) <<  b[i] << "," << setw(7) << R_lsh[i+1] << "," 
    << setw(10) << C/p_lsh[i] << "," << setw(10) << p_nc[i+1] << "," << setw(10) << C / p_rep[i] << "," << setw(9) << R[i+1] << "," << setw(9) << expected_list_size[i] << endl;
  }
  
  std::cout << "lv " << t << ": " 
  << setw(5) << w[t] << "," << setw(6) << "-," << setw(6) <<  "-," << setw(6) << b[t] << ",      -,"
  << setw(10) << C/p_lsh[t] << "," << setw(11) << "-," << setw(11) << "-," << "        -," << setw(9) << expected_list_size[t] << endl;  
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
  vector<double> p_r_gauss{0.864, 0.857, 0.854, 0.849, 0.845, 0.842, 0.836, 0.833, 0.830, 0.826, 
  0.822, 0.819, 0.814, 0.812, 0.809, 0.806, 0.802, 0.799, 0.796, 0.792, 0.790};
  vector<double> p_r_unif{0.815, 0.774, 0.742, 0.706, 0.675, 0.643, 0.617, 0.592, 0.57, 0.549,
  0.529, 0.509, 0.490, 0.472, 0.455, 0.437, 0.424, 0.409, 0.396, 0.382, 0.369};

  double p_suc_theory = p_r_gauss[r[1]/5] * p_r_gauss[r[0]/5] * p_nc[1];

  // cout << "p_suc_theory: " << p_r_gauss[r[1]/5] << " * " << p_r_gauss[r[0]/5] << " * " << p_nc[1] << endl;
  for (int i = 1; i < t; i++)
  {
    // cout << " * (" << pow(1 - exp(-3), 2) << " * " << p_r_unif[r[i+1]/5] << " * " << p_nc[i+1] << ")^" << (1 << i) << endl;
    p_suc_theory *= pow((1 - exp(-3)) * p_r_unif[r[i+1]/5] * p_nc[i+1], (1 << i));
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
    if (i == t)
    {
      std::cout << "- Step 1. Construct L" << i << " ⊂ {s, [Ms]_{B, " << r[i-1] << "}: HW(s) = " << w[i] << "} of size " << top_level_size << " (full size = " << top_level_full_size << ")" << endl;;
    }
    else 
    {
      std::cout << "- Step 2. Recover S" << i-1 << " = {s: HW(s) = " << w[i-1] 
      << " and [Ms]_{B, " << r[i-1] << "} in [" << -ell[i-1] << ", " << ell[i-1] 
      << "]^" << r[i-1] << "} (NCF with block-length " << b[i] << ")" << endl;  
    }
  }
  cout << endl;

  string prefix = to_string(d) + "_" + to_string(h) + "_" + to_string(w[1]) + "_" + to_string(w[2]);
  
  auto start = chrono::steady_clock::now();

  ///////////////////////////////////////////////////////////////////////////////////////// Original (RAM)

  for (size_t iter = 0; iter < repeat; iter++) 
  {
    std::cout << "... Runnning " << iter + 1 << "-th, until now " 
      << success_count << " successes" << '\r' << std::flush;
    
    matrix M; matrix B; secret s; vector<double> e;
    gen_noisy_instance(
      m, d, h, stddev, q, 
      M, B, s, e);
    auto Mtrans = transpose(M);

    set<secret> candidate_S;

    set<secret> basic_secret = enumerate_secrets(w[t], w[t]);
      vector<secret> basic_secret_vec(basic_secret.begin(), basic_secret.end());
      random_device rd;
      mt19937 gen(rd());
      uniform_int_distribution<int64_t> randindex(0, basic_secret_vec.size() - 1);

      while (candidate_S.size() < top_level_size) 
      {
        auto index = randindex(gen);
        secret tmp(basic_secret_vec[index]);
        tmp.resize(d);
        shuffle(tmp.begin(), tmp.end(), gen);
        candidate_S.insert(tmp);
      }
  
    for (int i = t; i >= 1; i--) 
    {
      list L;

      for (auto s_cand: candidate_S)
      {
      auto Ms = babaiNP(matmul(Mtrans, s_cand), B, r[i-1]);
      L.push_back(make_pair(s_cand, Ms));
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
  ///////////////////////////////////////////////////////////////////////////////////////// 

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