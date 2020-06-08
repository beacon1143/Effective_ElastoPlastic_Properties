#ifndef __UTILS_HPP__INCLUDED__
#define __UTILS_HPP__INCLUDED__

#include <vector>

namespace EFF_PROPS {

  template <typename T>
  void ResizeXY(std::vector<std::vector<T>>& matr, const int n1, const int n2) {
    matr.resize(n1);
    for (auto& vec: matr) {
      vec.resize(n2);
    }
  }

  template <typename T>
  T GetMaxElement(std::vector<std::vector<T>>& matr) {
    T retValue = std::numeric_limits<T>::min();

    for (const auto& vec : matr) {
      for (const auto& el : vec) {
        if (el > retValue) {
          retValue = el;
        }
      }
    }

    return retValue;
  }

  template <typename T>
  void SetMatrix(std::vector<std::vector<T>>& matr, const T value = static_cast<T>(0)) {
    for (auto& vec : matr) {
      for (auto& el : vec) {
        el = value;
      }
    }
  }

  template <typename T>
  void AverageOverFourPoints(const std::vector<std::vector<T>>& init, std::vector<std::vector<T>>& aver, 
                             const size_t startX = 0, const size_t startY = 0) {
    const size_t n1 = init.size();
    if (n1 < 2) {
      throw std::runtime_error("Error! First dimension of init is too small!\n");
    }
    const size_t n2 = init[0].size();
    for (size_t i = 1; i < n1; i++) {
      if (init[i].size() != n2) {
        throw std::runtime_error("Error! Dimensions of init columns are not consistent!\n");
      }
    }

    if (aver.size() != n1 - 1 + 2 * startX) {
      throw std::runtime_error("Error! First dimension of aver is bad!\n");
    }
    for (size_t i = 2; i < n1 - 1; i++) {
      if (aver[i].size() != n2 - 1 + 2 * startY) {
        throw std::runtime_error("Error! Second dimension of aver is bad!\n");
      }
    }

#pragma omp parallel for
    for (int i = 0; i < n1 - 1; i++) {
      for (int j = 0; j < n2 - 1; j++) {
        aver[startX + i][startY + j] = 0.25 * (init[i][j] + init[i+1][j] + init[i][j+1] + init[i+1][j+1]);
      }
    }
  }

} // namespace

#endif    // __UTILS_HPP__INCLUDED__