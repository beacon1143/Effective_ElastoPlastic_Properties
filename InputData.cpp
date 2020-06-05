#include "InputData.hpp"

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

  Kmax = std::numeric_limits<double>::min();
  Gmax = std::numeric_limits<double>::min();
  rho_max = std::numeric_limits<double>::min();
  for (int i = 0; i < nX; i++) {
    for (int j = 0; j < nY; j++) {
      if (K[i][j] > Kmax) {
        Kmax = K[i][j];
      }
      if (G[i][j] > Gmax) {
        Gmax = G[i][j];
      }
      if (rho[i][j] > rho_max) {
        rho_max = rho[i][j];
      }
    }
  }

  dT = courant * std::min(dX, dY) / sqrt( (Kmax + 4.0 * Gmax / 3.0) / rho_max );
  damp = 4.0 / dT / nX;
}

template<typename T>
static void InputData::ResizeXY(std::vector<std::vector<T>>& matr, const int n1, const int n2) {
  matr.resize(n1);
  for (auto& vec: matr) {
    vec.resize(n2);
  }
}