#include "torus_lsh.h"
#include "utils.h"

inline vector<int32_t> TorusLSH::compute_address(
  vector<double>& value, vector<double>& partition)
{
  vector<int32_t> address(lshDim_);
  for (size_t i = 0; i < lshDim_; i++) {
    int32_t a = floor((value[i] + domain_[i] - partition[i]) / lshLen_);
    address[i] = a % numBox_[i];
  }
  return address;
}

void TorusLSH::setPartitions(size_t repeats) 
{
  random_device rd;
  mt19937 gen(rd());

  repeats_ = repeats;

  partitions_.resize(repeats);
  for (size_t i = 0; i < repeats; i++) {
    vector<double> lsh_partition(lshDim_);
    for (size_t i = 0; i < lshDim_; i++) {
      uniform_real_distribution<double> lsh_random_starts(0, lshLen_);
      lsh_partition[i] = lsh_random_starts(gen);
    }
    partitions_[i] = lsh_partition;
  }
}

vector<HashTable> TorusLSH::computeLshTable(
  set<vector<double>>& list)
{
  vector<HashTable> hashTables;
  for (size_t i = 0; i < repeats_; i++)
  {
    HashTable hashTable;
    for (auto v: list) {
      auto addr = compute_address(v, partitions_[i]);
      hashTable[addr].insert(v);
    }
    hashTables.push_back(hashTable);
  }
  return hashTables;
}

int TorusLSH::searchFromLshTable(
  vector<double>& value,
  vector<HashTable>& hashTables, double dist,
  vector<double>& out)
{
  for (size_t i = 0; i < repeats_; i++) 
  {
    auto addr = compute_address(value, partitions_[i]);
    auto bin = hashTables[i][addr];
    for (auto v: bin) {
      auto norm_check = subvec(value, v);
      if (inf_norm(norm_check) <= dist) {
        out = v;
        return 1;
      }
    }
  }
  return 0;
}