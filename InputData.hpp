#ifndef __INPUT_DATA_HPP__INCLUDED__
#define __INPUT_DATA_HPP__INCLUDED__

#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <array>
#include "Utils.hpp"

namespace EFF_PROPS {

  class InputData {
  public:
    double sizeX, sizeY;    // physical size
    int nX, nY;             // numbers of space steps
    long int nTimeSteps;
    double courant;         // Courant–Friedrichs–Lewy

    double dX, dY;          // space steps

    InputData() = delete;
    explicit InputData(const double size_x, const double size_y,
                       const int nX, const int nY,
                       const long int nTimeSteps, const double courant);
  private:
  };

} // namespace

#endif    // __INPUT_DATA_HPP__INCLUDED__