#include "InputData.hpp"

namespace EFF_PROPS {

  void InputData::SetMaterials() {
    for (auto& vec : E) {
      for (auto& el : vec) {
        el = 1.0;
      }
    }
    for (auto& vec : nu) {
      for (auto& el : vec) {
        el = 0.25;
      }
    }
    for (auto& vec : rho) {
      for (auto& el : vec) {
        el = 1.0;
      }
    }
  }

  InputData::InputData(const double sizeX_, const double sizeY_,
                       const int nX_, const int nY_,
                       const long int nTimeSteps_, const double courant_) {
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
    if (nY_ <= 0) {
      throw std::runtime_error("Error! nY is not positive!\n");
    }
    nY = nY_;
    if (nTimeSteps_ <= 0) {
      throw std::runtime_error("Error! nTimeSteps is not positive!\n");
    }
    nTimeSteps = nTimeSteps_;
    if (courant_ <= 0.0) {
      throw std::runtime_error("Error! courant is not positive!\n");
    }
    courant = courant_;

    dX = sizeX / (nX - 1);
    dY = sizeY / (nY - 1);

    /* MATERIAL PROPERTIES */
    ResizeXY(E, nX, nY);
    ResizeXY(nu, nX, nY);
    ResizeXY(rho, nX, nY);

    SetMaterials();

    K.resize(nX);
    G.resize(nX);
    for (int i = 0; i < nX; i++) {
      K[i].resize(nY);
      G[i].resize(nY);
      for (int j = 0; j < nY; j++) {
        K[i][j] = E[i][j] / (3.0 - 6.0 * nu[i][j]);
        G[i][j] = E[i][j] / (2.0 + 2.0 * nu[i][j]);
      }
    }

    ResizeXY(Gav, nX - 1, nY - 1);
    AverageOverFourPoints(G, Gav);

    Kmax = GetMaxElement(K);
    Gmax = GetMaxElement(G);
    rho_max = GetMaxElement(rho);

    /*std::cout << "Gmax = " << Gmax << std::endl;
    std::cout << "Kmax = " << Kmax << std::endl;*/

    dT = courant * std::min(dX, dY) / sqrt( (Kmax + 4.0 * Gmax / 3.0) / rho_max );
    /*std::cout << "dT = " << dT << std::endl;*/
    damp = 4.0 / dT / nX;
  }

} // namespace