#ifndef __SIGMA_CALC_HPP__INCLUDED__
#define __SIGMA_CALC_HPP__INCLUDED__

#include <array>
#include <memory>
#include "InputData.hpp"

namespace EFF_PROPS {

  class SigmaCalc {
  public:
    std::unique_ptr<InputData> inp;

    // coordinates
    std::vector<std::vector<double>> x, y;
    std::vector<std::vector<double>> xUx, yUx;
    std::vector<std::vector<double>> xUy, yUy;

    // material properties
    double Kmax, Gmax, rho_max;
    std::vector<std::vector<double>> E, nu, K, G, rho;
    std::vector<std::vector<double>> Gav;
    void SetMaterials();

    // numeric
    double dT;              // time step
    [[deprecated("only Nx, not Ny, is used in formula")]] double damp;            // damping factor

    // variables
    std::vector<std::vector<double>> Ux, Uy;                 // displacement
    std::vector<std::vector<double>> divU;                   // displacement divergence
    std::vector<std::vector<double>> Vx, Vy;                 // velocity
    std::vector<std::vector<double>> Pinit, P;               // hydrostatic stress (ball part of tensor)
    std::vector<std::vector<double>> tauXX, tauYY, tauXY;    // deviatoric stress

    // effective stress
    std::array<double, 3> Sigma;

    void ComputeSigma(const double loadValue, const std::array<int, 3>& loadType);

    SigmaCalc() = delete;
    explicit SigmaCalc(std::unique_ptr<InputData> inp_);

  private:
    void ComputeDivergence(const std::vector<std::vector<double>>& Ax,
                           const std::vector<std::vector<double>>& Ay,
                           std::vector<std::vector<double>>& divA) const;
  };

} // namespace

#endif    // __SIGMA_CALC_HPP__INCLUDED__