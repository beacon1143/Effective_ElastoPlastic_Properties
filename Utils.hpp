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

} // namespace

#endif    // __UTILS_HPP__INCLUDED__