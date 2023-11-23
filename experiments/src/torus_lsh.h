#pragma once

#include "utils.h"
#include <set>
#include <map>
#include <random>

using HashTable = std::map<std::vector<int32_t>, std::set<std::vector<double>>>;

class TorusLSH {
  private:
    double lshLen_;
    double lshDim_;
    vector<double> domain_;
    
    vector<size_t> numBox_; 

    vector<vector<double>> partitions_;
    size_t repeats_ = 0;
    
  public:
    TorusLSH(double lshLen, vector<double> domain) 
    : lshLen_(lshLen), domain_(domain)
    {
      lshDim_ = domain_.size();
      numBox_.resize(lshDim_);
      for (size_t i = 0; i < lshDim_; i++) {
        numBox_[i] = ceil(2*domain_[i] / lshLen_);
      }
    }
    ~TorusLSH() = default;

    inline std::vector<int32_t> compute_address(
      std::vector<double>& value, 
      vector<double>& partition);

    void setPartitions(size_t repeats);

    std::vector<HashTable> computeLshTable(std::set<std::vector<double>>& list);

    int searchFromLshTable(
      std::vector<double>& value,
      std::vector<HashTable>& hashTables, double dist,
      std::vector<double>& out);
};