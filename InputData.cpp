#include "InputData.hpp"

namespace EFF_PROPS {

  InputData::InputData(const double sizeX_, const double sizeY_,
                       const int nX_, const int nY_,
                       const size_t nTimeSteps_, const size_t nIterations_, const double courant_) {
    if (sizeX_ <= 0.0) {
      throw std::runtime_error("Error! sizeX is not positive!\n");
    }
    sizeX = sizeX_;

    if (sizeY_ <= 0.0) {
      throw std::runtime_error("Error! sizeY is not positive!\n");
    }
    sizeY = sizeY_;

    if (nX_ <= 0) {
      throw std::runtime_error("Error! nX is not positive!\n");
    }
    nX = nX_;
    nXp = nX + 1;
    nXm = nX - 1;

    if (nY_ <= 0) {
      throw std::runtime_error("Error! nY is not positive!\n");
    }
    nY = nY_;
    nYp = nY + 1;
    nYm = nY - 1;

    if (nTimeSteps_ <= 0) {
      throw std::runtime_error("Error! nTimeSteps is not positive!\n");
    }
    nTimeSteps = nTimeSteps_;
    if (nIterations_ <= 0) {
      throw std::runtime_error("Error! nIterations is not positive!\n");
    }
    nIterations = nIterations_;

    if (courant_ <= 0.0) {
      throw std::runtime_error("Error! courant is not positive!\n");
    }
    courant = courant_;

    dX = sizeX / (nX - 1);
    dY = sizeY / (nY - 1);
  }

} // namespace