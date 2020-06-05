#include "SigmaCalc.hpp"

namespace EFF_PROPS {

  SigmaCalc::SigmaCalc(std::unique_ptr<InputData> inp_) {
    if (inp_ == nullptr) {
      throw std::runtime_error("Error! inp is nullptr!\n");
    }
    inp = std::move(inp_);

    /* VARIABLES */
    x.resize(inp->nX);
    y.resize(inp->nX);
    for (int i = 0; i < inp->nX; i++) {
      x[i].resize(inp->nY);
      y[i].resize(inp->nY);
      for (int j = 0; j < inp->nY; j++) {
        x[i][j] = -0.5 * inp->sizeX + inp->dX * i;
        y[i][j] = -0.5 * inp->sizeY + inp->dY * j;
      }
    }

    xUx.resize(inp->nX + 1);
    yUx.resize(inp->nX + 1);
    for (int i = 0; i < inp->nX + 1; i++) {
      xUx[i].resize(inp->nY);
      yUx[i].resize(inp->nY);
      for (int j = 0; j < inp->nY; j++) {
        xUx[i][j] = -0.5 * (inp->sizeX + inp->dX) + inp->dX * i;
        yUx[i][j] = -0.5 * inp->sizeY + inp->dY * j;
      }
    }

    xUy.resize(inp->nX);
    yUy.resize(inp->nX);
    for (int i = 0; i < inp->nX; i++) {
      xUy[i].resize(inp->nY + 1);
      yUy[i].resize(inp->nY + 1);
      for (int j = 0; j < inp->nY + 1; j++) {
        xUy[i][j] = -0.5 * inp->sizeX + inp->dX * i;
        yUy[i][j] = -0.5 * (inp->sizeY + inp->dY) + inp->dY * j;
      }
    }

    /* VARIABLES */
    // displacement
    Ux.resize(inp->nX + 1);
    for (auto& vec : Ux) {
      vec.resize(inp->nY, 0.0);
    }
    Uy.resize(inp->nX);
    for (auto& vec : Uy) {
      vec.resize(inp->nY + 1, 0.0);
    }

    // displacement divergence
    divU.resize(inp->nX);
    for (auto& vec : divU) {
      vec.resize(inp->nY, 0.0);
    }

    // velocity
    Vx.resize(inp->nX + 1);
    for (auto& vec : Vx) {
      vec.resize(inp->nY, 0.0);
    }
    Vy.resize(inp->nX);
    for (auto& vec : Vy) {
      vec.resize(inp->nY + 1, 0.0);
    }

    // hydrostatic stress (ball part of tensor)
    Pinit.resize(inp->nX);
    for (auto& vec : Pinit) {
      vec.resize(inp->nY, 0.0);
    }
    P.resize(inp->nX);
    for (auto& vec : P) {
      vec.resize(inp->nY, 0.0);
    }

    // deviatoric stress
    tauXX.resize(inp->nX);
    for (auto& vec : tauXX) {
      vec.resize(inp->nY, 0.0);
    }
    tauYY.resize(inp->nX);
    for (auto& vec : tauYY) {
      vec.resize(inp->nY, 0.0);
    }
    tauXY.resize(inp->nX - 1);
    for (auto& vec : tauXY) {
      vec.resize(inp->nY - 1, 0.0);
    }

  }

  void SigmaCalc::ComputeSigma(const double loadValue, const std::array<int, 3>& loadType) {
    // boundary conditions
    const double dUxdx = loadValue * loadType[0];
    const double dUydy = loadValue * loadType[1];
    const double dUxdy = loadValue * loadType[2];

    // initial conditions
    for (int i = 0; i < inp->nX + 1; i++) {
      for (int j = 0; j < inp->nY; j++) {
        Ux[i][j] = dUxdx * xUx[i][j] + dUxdy * yUx[i][j];
      }
    }
    for (int i = 0; i < inp->nX; i++) {
      for (int j = 0; j < inp->nY + 1; j++) {
        Uy[i][j] = dUydy * yUy[i][j];
      }
    }

    SetMatrix(Vx, 0.0);
    SetMatrix(Vy, 0.0);
    SetMatrix(P, 0.0);
    SetMatrix(tauXX, 0.0);
    SetMatrix(tauYY, 0.0);
    SetMatrix(tauXY, 0.0);

    for (size_t it = 0; it < inp->nTimeSteps; it++) {
      // displacement divergence
      ComputeDivergence(Ux, Uy, divU);

      // constitutive equation - Hooke's law
      for (int i = 0; i < inp->nX; i++) {
        for (int j = 0; j < inp->nY; j++) {
          P[i][j] = Pinit[i][j] - inp->K[i][j] * divU[i][j];
          tauXX[i][j] = 2.0 * inp->G[i][j] * ( (Ux[i+1][j] - Ux[i][j]) / inp->dX - divU[i][j]/3.0);
          tauYY[i][j] = 2.0 * inp->G[i][j] * ( (Uy[i][j+1] - Uy[i][j]) / inp->dY - divU[i][j]/3.0);
          if (i < inp->nX - 1 && j < inp->nY - 1) {
            tauXY[i][j] = inp->Gav[i][j] * ( (Ux[i+1][j+1] - Ux[i+1][j]) / inp->dY + (Uy[i+1][j+1] - Uy[i][j+1]) / inp->dX );
          }
        }
      }

      // motion equation
      for (int i = 1; i < inp->nX; i++) {
        for (int j = 1; j < inp->nY - 1; j++) {
          Vx[i][j] = Vx[i][j] * (1.0 - inp->dT * inp->damp) + (
                     (-P[i][j] + P[i-1][j] + tauXX[i][j] - tauXX[i-1][j]) / inp->dX / inp->rho_max +
                     (tauXY[i-1][j] - tauXY[i-1][j-1]) / inp->dY
                     ) * inp->dT;
        }
      }

      for (int i = 1; i < inp->nX - 1; i++) {
        for (int j = 1; j < inp->nY; j++) {
          Vy[i][j] = Vy[i][j] * (1.0 - inp->dT * inp->damp) + (
                     (-P[i][j] + P[i][j-1] + tauYY[i][j] - tauYY[i][j-1]) / inp->dY / inp->rho_max +
                     (tauXY[i][j-1] - tauXY[i-1][j-1]) / inp->dX
                     ) * inp->dT;
        }
      }

      /*const double Vxmax = GetMaxElement(Vx);
      const double Vymax = GetMaxElement(Vy);
      std::cout << "Vxmax = " << Vxmax << "\nVymax = " << Vymax << "\n";*/

      // displacement
      for (int i = 0; i < inp->nX + 1; i++) {
        for (int j = 0; j < inp->nY; j++) {
          Ux[i][j] = Ux[i][j] + Vx[i][j] * inp->dT;
        }
      }
      for (int i = 0; i < inp->nX; i++) {
        for (int j = 0; j < inp->nY + 1; j++) {
          Uy[i][j] = Uy[i][j] + Vy[i][j] * inp->dT;
        }
      }
    } // time loop

    /*std::cout << "tauYY\n";
    for (int i = 0; i < inp->nX; i++) {
      for (int j = 0; j < inp->nY; j++) {
        std::cout << tauYY[i][j] << ' ';
      }
      std::cout << '\n';
    }*/

    // averaging
    Sigma[0] = 0.0;
    for (int i = 0; i < inp->nX; i++) {
      for (int j = 0; j < inp->nY; j++) {
        Sigma[0] += tauXX[i][j] - P[i][j];
      }
    }
    Sigma[0] /= static_cast<double>(inp->nX * inp->nY);

    Sigma[1] = 0.0;
    for (int i = 0; i < inp->nX; i++) {
      for (int j = 0; j < inp->nY; j++) {
        Sigma[1] += tauYY[i][j] - P[i][j];
      }
    }
    Sigma[1] /= static_cast<double>(inp->nX * inp->nY);

    Sigma[2] = 0.0;
    for (int i = 0; i < inp->nX - 1; i++) {
      for (int j = 0; j < inp->nY - 1; j++) {
        Sigma[2] += tauXY[i][j];
      }
    }
    Sigma[2] /= static_cast<double>((inp->nX - 1) * (inp->nY - 1));

    /*std::cout << "Sigma\n" << Sigma[0] << ' ' << Sigma[1] << ' ' << Sigma[2] << '\n';*/
  }


  void SigmaCalc::ComputeDivergence(const std::vector<std::vector<double>>& Ax,
                                    const std::vector<std::vector<double>>& Ay,
                                    std::vector<std::vector<double>>& divA) {
    size_t length = Ax.size() - 1;
    for (const auto& vec : Ax) {
      if (vec.size() != length) {
        throw std::runtime_error("Error! Wrong dimensions in ComputeDivergence()!\n");
      }
    }
    if (Ay.size() != length) {
      throw std::runtime_error("Error! Wrong dimensions in ComputeDivergence()!\n");
    }
    for (const auto& vec : Ay) {
      if (vec.size() != length + 1) {
        throw std::runtime_error("Error! Wrong dimensions in ComputeDivergence()!\n");
      }
    }
    if (divA.size() != length) {
      throw std::runtime_error("Error! Wrong dimensions in ComputeDivergence()!\n");
    }
    for (const auto& vec : divA) {
      if (vec.size() != length) {
        throw std::runtime_error("Error! Wrong dimensions in ComputeDivergence()!\n");
      }
    }

    for (int i = 0; i < length; i++) {
      for (int j = 0; j < length; j++) {
        divA[i][j] = (Ax[i+1][j] - Ax[i][j]) / inp->dX + (Ay[i][j+1] - Ay[i][j]) / inp->dY;
      }
    }
  }

} // namespace