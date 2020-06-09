#include "SigmaCalc.hpp"

namespace EFF_PROPS {

  void SigmaCalc::SetMaterials() {
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

    /*for (EP_INT i = 0; i < inp->nX; i++) {
      for (EP_INT j = 0; j < inp->nY; j++) {
        if (sqrt(x[i][j] * x[i][j] + y[i][j] * y[i][j]) < 2.85459861019) {
          E[i][j] = 2.0;
          nu[i][j] = 0.2;
        }
        else {
          E[i][j] = 0.002;
          nu[i][j] = 0.3;
        }
      }
    }*/

    for (auto& vec : rho) {
      for (auto& el : vec) {
        el = 2.0;
      }
    }
  }

  SigmaCalc::SigmaCalc(std::unique_ptr<InputData> inp_) {
    if (inp_ == nullptr) {
      throw std::runtime_error("Error! inp is nullptr!\n");
    }
    inp = std::move(inp_);

    /* VARIABLES */
    x.resize(inp->nX);
    y.resize(inp->nX);
    for (EP_INT i = 0; i < inp->nX; i++) {
      x[i].resize(inp->nY);
      y[i].resize(inp->nY);
      for (EP_INT j = 0; j < inp->nY; j++) {
        x[i][j] = -0.5 * inp->sizeX + inp->dX * i;
        y[i][j] = -0.5 * inp->sizeY + inp->dY * j;
      }
    }

    xUx.resize(inp->nX + 1);
    yUx.resize(inp->nX + 1);
    for (EP_INT i = 0; i < inp->nX + 1; i++) {
      xUx[i].resize(inp->nY);
      yUx[i].resize(inp->nY);
      for (EP_INT j = 0; j < inp->nY; j++) {
        xUx[i][j] = -0.5 * (inp->sizeX + inp->dX) + inp->dX * i;
        yUx[i][j] = -0.5 * inp->sizeY + inp->dY * j;
      }
    }

    xUy.resize(inp->nX);
    yUy.resize(inp->nX);
    for (EP_INT i = 0; i < inp->nX; i++) {
      xUy[i].resize(inp->nY + 1);
      yUy[i].resize(inp->nY + 1);
      for (EP_INT j = 0; j < inp->nY + 1; j++) {
        xUy[i][j] = -0.5 * inp->sizeX + inp->dX * i;
        yUy[i][j] = -0.5 * (inp->sizeY + inp->dY) + inp->dY * j;
      }
    }

    /* MATERIAL PROPERTIES */
    ResizeXY(E, inp->nX, inp->nY);
    ResizeXY(nu, inp->nX, inp->nY);
    ResizeXY(rho, inp->nX, inp->nY);

    SetMaterials();

    K.resize(inp->nX);
    G.resize(inp->nX);
    for (EP_INT i = 0; i < inp->nX; i++) {
      K[i].resize(inp->nY);
      G[i].resize(inp->nY);
      for (EP_INT j = 0; j < inp->nY; j++) {
        K[i][j] = E[i][j] / (3.0 - 6.0 * nu[i][j]);
        G[i][j] = E[i][j] / (2.0 + 2.0 * nu[i][j]);
      }
    }

    ResizeXY(Gav, inp->nX - 1, inp->nY - 1);
    AverageOverFourPoints(G, Gav);

    Kmax = GetMaxElement(K);
    Gmax = GetMaxElement(G);
    rho_max = GetMaxElement(rho);
    cohesion = 0.00075;

    /*std::cout << "Gmax = " << Gmax << std::endl;
    std::cout << "Kmax = " << Kmax << std::endl;*/

    /* NUMERIC */
    dT = inp->courant * std::min(inp->dX, inp->dY) / sqrt( (Kmax + 4.0 * Gmax / 3.0) / rho_max );
    /*std::cout << "dT = " << dT << std::endl;*/
    damp = 4.0 / dT / inp->nX;

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

    // deviatoric stress for plasticity criteria
    tauXXav.resize(inp->nX - 1);
    for (auto& vec : tauXXav) {
      vec.resize(inp->nY - 1, 0.0);
    }
    tauYYav.resize(inp->nX - 1);
    for (auto& vec : tauYYav) {
      vec.resize(inp->nY - 1, 0.0);
    }
    tauXYav.resize(inp->nX);
    for (auto& vec : tauXYav) {
      vec.resize(inp->nY, 0.0);
    }

    Sigma.resize(inp->nTimeSteps);
  }

  void SigmaCalc::ComputeSigma(const EP_FLOAT loadValue, const std::array<EP_INT, 3>& loadType) {
    // boundary conditions
    const EP_FLOAT dUxdx = loadValue * loadType[0];
    const EP_FLOAT dUydy = loadValue * loadType[1];
    const EP_FLOAT dUxdy = loadValue * loadType[2];

    for (EP_INT i = 0; i < inp->nX + 1; i++) {
      for (EP_INT j = 0; j < inp->nY; j++) {
        Ux[i][j] = 0.0; //dUxdx * xUx[i][j] + dUxdy * yUx[i][j];
      }
    }
    for (EP_INT i = 0; i < inp->nX; i++) {
      for (EP_INT j = 0; j < inp->nY + 1; j++) {
        Uy[i][j] = 0.0; //dUydy * yUy[i][j];
      }
    }

    SetMatrix(Vx, 0.0);
    SetMatrix(Vy, 0.0);
    SetMatrix(P, 0.0);
    SetMatrix(tauXX, 0.0);
    SetMatrix(tauYY, 0.0);
    SetMatrix(tauXY, 0.0);

    for (size_t tim = 0; tim < inp->nTimeSteps; tim++) {
      // initial conditions
      for (EP_INT i = 0; i < inp->nX + 1; i++) {
        for (EP_INT j = 0; j < inp->nY; j++) {
          Ux[i][j] += (dUxdx * xUx[i][j] + dUxdy * yUx[i][j]) / inp->nTimeSteps;
        }
      }
      for (EP_INT i = 0; i < inp->nX; i++) {
        for (EP_INT j = 0; j < inp->nY + 1; j++) {
          Uy[i][j] += (dUydy * yUy[i][j]) / inp->nTimeSteps;
        }
      }

      for (size_t it = 0; it < inp->nIterations; it++) {
        // displacement divergence
        ComputeDivergence(Ux, Uy, divU);

        // constitutive equation - Hooke's law
#pragma omp parallel for
        for (EP_INT i = 0; i < inp->nX; i++) {
          for (EP_INT j = 0; j < inp->nY; j++) {
            P[i][j] = Pinit[i][j] - K[i][j] * divU[i][j];
            tauXX[i][j] = 2.0 * G[i][j] * ( (Ux[i+1][j] - Ux[i][j]) / inp->dX - divU[i][j]/3.0);
            tauYY[i][j] = 2.0 * G[i][j] * ( (Uy[i][j+1] - Uy[i][j]) / inp->dY - divU[i][j]/3.0);
          }
        }

#pragma omp parallel for
        for (EP_INT i = 0; i < inp->nXm; i++) {
          for (EP_INT j = 0; j < inp->nYm; j++) {
            tauXY[i][j] = Gav[i][j] * ( (Ux[i+1][j+1] - Ux[i+1][j]) / inp->dY + (Uy[i+1][j+1] - Uy[i][j+1]) / inp->dX );
          }
        }

        // tauXY for plasticity
        AverageOverFourPoints(tauXY, tauXYav, 1, 1);
#pragma omp parallel for
        for (EP_INT i = 1; i < inp->nXm; i++) {
          tauXYav[i][0] = tauXYav[i][1];
          tauXYav[i][inp->nXm] = tauXYav[i][inp->nXm - 1];
        }
#pragma omp parallel for
        for (EP_INT i = 1; i < inp->nYm; i++) {
          tauXYav[0][i] = tauXYav[1][i];
          tauXYav[inp->nXm][i] = tauXYav[inp->nXm - 1][i];
        }
        tauXYav[0][0] = 0.5 * (tauXYav[1][0] + tauXYav[0][1]);
        tauXYav[inp->nXm][0] = 0.5 * (tauXYav[inp->nXm][1] + tauXYav[inp->nXm - 1][0]);
        tauXYav[0][inp->nYm] = 0.5 * (tauXYav[1][inp->nYm] + tauXYav[0][inp->nYm - 1]);
        tauXYav[inp->nXm][inp->nYm] = 0.5 * (tauXYav[inp->nXm][inp->nYm - 1] + tauXYav[inp->nXm - 1][inp->nYm]);

        AverageOverFourPoints(tauXX, tauXXav);
        AverageOverFourPoints(tauYY, tauYYav);

        // plasticity
#pragma omp parallel for
        for (EP_INT i = 0; i < inp->nX; i++) {
          for (EP_INT j = 0; j < inp->nY; j++) {
            // Tresca criteria
            EP_FLOAT j2 = sqrt(tauXX[i][j] * tauXX[i][j] + tauYY[i][j] * tauYY[i][j] + 2.0 * tauXYav[i][j] * tauXYav[i][j]);
            if (j2 > cohesion) {
              tauXX[i][j] *= cohesion / j2;
              tauYY[i][j] *= cohesion / j2;
            }
          }
        }

#pragma omp parallel for
        for (EP_INT i = 0; i < inp->nXm; i++) {
          for (EP_INT j = 0; j < inp->nYm; j++) {
            EP_FLOAT j2xy = sqrt(tauXXav[i][j] * tauXXav[i][j] + tauYYav[i][j] * tauYYav[i][j] + 2.0 * tauXY[i][j] * tauXY[i][j]);
            if (j2xy > cohesion) {
              tauXY[i][j] *= cohesion / j2xy;
            }
          }
        }

        // motion equation
#pragma omp parallel for
        for (EP_INT i = 1; i < inp->nX; i++) {
          for (EP_INT j = 1; j < inp->nYm; j++) {
            Vx[i][j] = Vx[i][j] * (1.0 - dT * damp) + (
                       (-P[i][j] + P[i-1][j] + tauXX[i][j] - tauXX[i-1][j]) / inp->dX / rho_max +
                       (tauXY[i-1][j] - tauXY[i-1][j-1]) / inp->dY
                       ) * dT;
          }
        }

#pragma omp parallel for
        for (EP_INT i = 1; i < inp->nXm; i++) {
          for (EP_INT j = 1; j < inp->nY; j++) {
            Vy[i][j] = Vy[i][j] * (1.0 - dT * damp) + (
                       (-P[i][j] + P[i][j-1] + tauYY[i][j] - tauYY[i][j-1]) / inp->dY / rho_max +
                       (tauXY[i][j-1] - tauXY[i-1][j-1]) / inp->dX
                       ) * dT;
          }
        }

        if ((it+1) % 1'000 == 0) {
          std::cout << "Iteration " << it + 1 << " from " << inp->nIterations << "\n";
          const EP_FLOAT Vxmax = GetMaxElement(Vx);
          const EP_FLOAT Vymax = GetMaxElement(Vy);
          std::cout << "Vxmax = " << Vxmax << "\tVymax = " << Vymax << "\n";
        }

        // displacement
#pragma omp parallel for
        for (EP_INT i = 0; i < inp->nXp; i++) {
          for (EP_INT j = 0; j < inp->nY; j++) {
            Ux[i][j] = Ux[i][j] + Vx[i][j] * dT;
          }
        }
#pragma omp parallel for
        for (EP_INT i = 0; i < inp->nX; i++) {
          for (EP_INT j = 0; j < inp->nYp; j++) {
            Uy[i][j] = Uy[i][j] + Vy[i][j] * dT;
          }
        }
      } // for (it)

      /*std::cout << "tauYY\n";
      for (EP_INT i = 0; i < inp->nX; i++) {
        for (EP_INT j = 0; j < inp->nY; j++) {
          std::cout << tauYY[i][j] << ' ';
        }
        std::cout << '\n';
      }*/

      // averaging
      Sigma[tim][0] = 0.0;
      for (EP_INT i = 0; i < inp->nX; i++) {
        for (EP_INT j = 0; j < inp->nY; j++) {
          Sigma[tim][0] += tauXX[i][j] - P[i][j];
        }
      }
      Sigma[tim][0] /= static_cast<EP_FLOAT>(inp->nX * inp->nY);

      Sigma[tim][1] = 0.0;
      for (EP_INT i = 0; i < inp->nX; i++) {
        for (EP_INT j = 0; j < inp->nY; j++) {
          Sigma[tim][1] += tauYY[i][j] - P[i][j];
        }
      }
      Sigma[tim][1] /= static_cast<EP_FLOAT>(inp->nX * inp->nY);

      Sigma[tim][2] = 0.0;
      for (EP_INT i = 0; i < inp->nX - 1; i++) {
        for (EP_INT j = 0; j < inp->nY - 1; j++) {
          Sigma[tim][2] += tauXY[i][j];
        }
      }
      Sigma[tim][2] /= static_cast<EP_FLOAT>((inp->nX - 1) * (inp->nY - 1));

      /*std::cout << "Sigma\n" << Sigma[0] << ' ' << Sigma[1] << ' ' << Sigma[2] << '\n';*/
    } // for (tim)
  }


  void SigmaCalc::ComputeDivergence(const std::vector<std::vector<EP_FLOAT>>& Ax,
                                    const std::vector<std::vector<EP_FLOAT>>& Ay,
                                    std::vector<std::vector<EP_FLOAT>>& divA) const {
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

    for (size_t i = 0; i < length; i++) {
      for (size_t j = 0; j < length; j++) {
        divA[i][j] = (Ax[i+1][j] - Ax[i][j]) / inp->dX + (Ay[i][j+1] - Ay[i][j]) / inp->dY;
      }
    }
  }

} // namespace
