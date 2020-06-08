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
    EP_FLOAT sizeX, sizeY;    // physical size
    EP_INT nX, nY;             // numbers of space steps
    EP_INT nXp, nYp;
    EP_INT nXm, nYm;
    size_t nTimeSteps;
    size_t nIterations;
    EP_FLOAT courant;         // Courant–Friedrichs–Lewy

    EP_FLOAT dX, dY;          // space steps

    InputData() = delete;
    explicit InputData(const EP_FLOAT size_x, const EP_FLOAT size_y,
                       const EP_INT nX, const EP_INT nY,
                       const size_t nTimeSteps, const size_t nIterations, const EP_FLOAT courant);
  private:
  };

} // namespace

#endif    // __INPUT_DATA_HPP__INCLUDED__