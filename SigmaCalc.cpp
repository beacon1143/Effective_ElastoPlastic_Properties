#include "SigmaCalc.hpp"

namespace EFF_PROPS {

  void SigmaCalc::SetMaterials() {
    /*for (auto& vec : E) {
      for (auto& el : vec) {
        el = 1.0;
      }
    }
    for (auto& vec : nu) {
      for (auto& el : vec) {
        el = 0.25;
      }
    }*/

    double K0 = 1.0;
    double G0 = 0.25;

    for (EP_INT i = 0; i < inp->nX; i++) {
      for (EP_INT j = 0; j < inp->nY; j++) {
        if (sqrt(x[i][j] * x[i][j] + y[i][j] * y[i][j]) < 1.0) {
          K[i][j] = 0.01 * K0;
          G[i][j] = 0.01 * G0;
        }
        else {
          K[i][j] = K0;
          G[i][j] = G0;
        }
      }
    }

    for (auto& vec : rho) {
      for (auto& el : vec) {
        el = 1.0;
      }
    }
  }

  void SigmaCalc::SetPressure(EP_FLOAT P0) {
    for (EP_INT i = 0; i < inp->nX; i++) {
      for (EP_INT j = 0; j < inp->nY; j++) {
        if (sqrt(x[i][j] * x[i][j] + y[i][j] * y[i][j]) < 1.0) {
          Pinit[i][j] = P0;
        }
        else {
          Pinit[i][j] = 0.0;
        }
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
    /*ResizeXY(E, inp->nX, inp->nY);
    ResizeXY(nu, inp->nX, inp->nY);*/
    ResizeXY(rho, inp->nX, inp->nY);

    ResizeXY(K, inp->nX, inp->nY);
    ResizeXY(G, inp->nX, inp->nY);

    SetMaterials();

    /*K.resize(inp->nX);
    G.resize(inp->nX);
    for (EP_INT i = 0; i < inp->nX; i++) {
      K[i].resize(inp->nY);
      G[i].resize(inp->nY);
      for (EP_INT j = 0; j < inp->nY; j++) {
        K[i][j] = E[i][j] / (3.0 - 6.0 * nu[i][j]);
        G[i][j] = E[i][j] / (2.0 + 2.0 * nu[i][j]);
      }
    }*/

    ResizeXY(Gav, inp->nX - 1, inp->nY - 1);
    AverageOverFourPoints(G, Gav);

    Kmax = GetMaxElement(K);
    Gmax = GetMaxElement(G);
    rho_max = GetMaxElement(rho);
    cohesion = 0.001;
    Ppore = 1.0 * cohesion;

    /*std::cout << "Gmax = " << Gmax << std::endl;
    std::cout << "Kmax = " << Kmax << std::endl;*/

    /* NUMERIC */
    dT = inp->courant * std::min(inp->dX, inp->dY) / sqrt( (Kmax + 4.0 * Gmax / 3.0) / rho_max );
    /*std::cout << "dT = " << dT << std::endl;*/
    dampX = 4.0 / dT / inp->nX;
    dampY = 4.0 / dT / inp->nY;

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
    SetPressure(Ppore);

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
    deltaP.resize(inp->nTimeSteps, 0.0);
    tauInfty.resize(inp->nTimeSteps, 0.0);
    Keff.resize(inp->nTimeSteps, 0.0);
    Geff.resize(inp->nTimeSteps);
  }

  void SigmaCalc::ComputeSigma(const EP_FLOAT loadValue, const std::array<EP_FLOAT, 3>& loadType) {
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

    // plasticity indicators
    std::vector<std::vector<EP_FLOAT>> iPlast(inp->nX);
    for (auto& vec : iPlast) {
      vec.resize(inp->nY, 0.0);
    }

    /* TIME LOOP */
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

      /* ITERATIONS LOOP */
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
              iPlast[i][j] = 1.0;
            }
          }
        }

#pragma omp parallel for
        for (EP_INT i = 0; i < inp->nXm; i++) {
          for (EP_INT j = 0; j < inp->nYm; j++) {
            EP_FLOAT j2xy = sqrt(tauXXav[i][j] * tauXXav[i][j] + tauYYav[i][j] * tauYYav[i][j] + 2.0 * tauXY[i][j] * tauXY[i][j]);
            if (j2xy > cohesion) {
              tauXY[i][j] *= cohesion / j2xy;
              iPlast[i][j] = 1.0;
              iPlast[i + 1][j] = 1.0;
              iPlast[i][j + 1] = 1.0;
              iPlast[i + 1][j + 1] = 1.0;
            }
          }
        }

        // motion equation
#pragma omp parallel for
        for (EP_INT i = 1; i < inp->nX; i++) {
          for (EP_INT j = 1; j < inp->nYm; j++) {
            Vx[i][j] = Vx[i][j] * (1.0 - dT * dampX) + (
                       (-P[i][j] + P[i-1][j] + tauXX[i][j] - tauXX[i-1][j]) / inp->dX / rho_max +
                       (tauXY[i-1][j] - tauXY[i-1][j-1]) / inp->dY
                       ) * dT;
          }
        }

#pragma omp parallel for
        for (EP_INT i = 1; i < inp->nXm; i++) {
          for (EP_INT j = 1; j < inp->nY; j++) {
            Vy[i][j] = Vy[i][j] * (1.0 - dT * dampY) + (
                       (-P[i][j] + P[i][j-1] + tauYY[i][j] - tauYY[i][j-1]) / inp->dY / rho_max +
                       (tauXY[i][j-1] - tauXY[i-1][j-1]) / inp->dX
                       ) * dT;
          }
        }

        if ((it+1) % 10'000 == 0) {
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

      deltaP[tim] = 0.0;
      for (EP_INT i = 0; i < inp->nX; i++) {
        deltaP[tim] += (tauXX[i][0] + tauYY[i][0] - 2.0 * P[i][0]) / inp->nX;
        deltaP[tim] += (tauXX[i][inp->nY - 1] + tauYY[i][inp->nY - 1] - 2.0 * P[i][inp->nY - 1]) / inp->nX;
      }
      for (EP_INT j = 0; j < inp->nY; j++) {
        deltaP[tim] += (tauXX[0][j] + tauYY[0][j] - 2.0 * P[0][j]) / inp->nY;
        deltaP[tim] += (tauXX[inp->nX - 1][j] + tauYY[inp->nX - 1][j] - 2.0 * P[inp->nX - 1][j]) / inp->nY;
      }
      deltaP[tim] *= 0.125 / cohesion / sqrt(2.0);

      tauInfty[tim] = 0.0;
      for (EP_INT i = 0; i < inp->nX; i++) {
        tauInfty[tim] += (tauXX[i][0] - tauYY[i][0]) / inp->nX;
        tauInfty[tim] += (tauXX[i][inp->nY - 1] - tauYY[i][inp->nY - 1]) / inp->nX;
      }
      for (EP_INT j = 0; j < inp->nY; j++) {
        tauInfty[tim] += (tauXX[0][j] - tauYY[0][j]) / inp->nY;
        tauInfty[tim] += (tauXX[inp->nX - 1][j] - tauYY[inp->nX - 1][j]) / inp->nY;
      }
      tauInfty[tim] *= 0.125 / cohesion / sqrt(2.0);

      EP_FLOAT divUeff = loadValue * (loadType[0] + loadType[1]);

      Keff[tim] = 0.0;
      for (EP_INT i = 0; i < inp->nX; i++) {
        for (EP_INT j = 0; j < inp->nY; j++) {
          Keff[tim] -= P[i][j];
        }
      }
      Keff[tim] /= static_cast<EP_FLOAT>(inp->nX * inp->nY) * divUeff * (tim + 1) / inp->nTimeSteps;

      Geff[tim][0] = 0.0;
      for (EP_INT i = 0; i < inp->nX; i++) {
        for (EP_INT j = 0; j < inp->nY; j++) {
          Geff[tim][0] += tauXX[i][j];
        }
      }
      Geff[tim][0] /= 2.0 * static_cast<EP_FLOAT>(inp->nX * inp->nY) * (loadValue * loadType[0] - divUeff / 3.0) * (tim + 1) / inp->nTimeSteps;

      Geff[tim][1] = 0.0;
      for (EP_INT i = 0; i < inp->nX; i++) {
        for (EP_INT j = 0; j < inp->nY; j++) {
          Geff[tim][1] += tauYY[i][j];
        }
      }
      Geff[tim][1] /= 2.0 * static_cast<EP_FLOAT>(inp->nX * inp->nY) * (loadValue * loadType[1] - divUeff / 3.0) * (tim + 1) / inp->nTimeSteps;

      std::cout << "deltaP = " << deltaP[tim] << '\n';
      std::cout << "tauInfty = " << tauInfty[tim] << '\n';
      std::cout << "Keff = " << Keff[tim] << '\n';
      std::cout << "GeffXX = " << Geff[tim][0] << '\n';
      std::cout << "GeffYY = " << Geff[tim][1] << '\n';

      std::ofstream fParams("../src/params.txt", std::ofstream::out);
      fParams << inp->sizeX << " " << inp->sizeY << " " << inp->nX << " " << inp->nY << '\n';
      fParams << Ppore << " " << cohesion << " " << loadValue << std::endl;
      fParams.close();

      std::ofstream fPres("../src/pressure.txt", std::ofstream::out);
      for (EP_INT j = 0; j < inp->nY; j++) {
        for (EP_INT i = 0; i < inp->nX; i++) {
          fPres << P[i][j] << " ";
        }
        fPres << '\n';
      }
      fPres.close();

      std::ofstream fTauXX("../src/tau_xx.txt", std::ofstream::out);
      for (EP_INT j = 0; j < inp->nY; j++) {
        for (EP_INT i = 0; i < inp->nX; i++) {
          fTauXX << tauXX[i][j] << " ";
        }
        fTauXX << '\n';
      }
      fTauXX.close();

      std::ofstream fTauYY("../src/tau_yy.txt", std::ofstream::out);
      for (EP_INT j = 0; j < inp->nY; j++) {
        for (EP_INT i = 0; i < inp->nX; i++) {
          fTauYY << tauYY[i][j] << " ";
        }
        fTauYY << '\n';
      }
      fTauYY.close();

      std::ofstream fPlast("../src/plast.txt", std::ofstream::out);
      for (EP_INT j = 0; j < inp->nY; j++) {
        for (EP_INT i = 0; i < inp->nX; i++) {
          fPlast << iPlast[i][j] << " ";
        }
        fPlast << '\n';
      }
      fPlast.close();

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
