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
    std::vector<std::vector<EP_FLOAT>> x, y;
    std::vector<std::vector<EP_FLOAT>> xUx, yUx;
    std::vector<std::vector<EP_FLOAT>> xUy, yUy;

    // material properties
    EP_FLOAT Kmax, Gmax, rho_max;
    std::vector<std::vector<EP_FLOAT>> /*E, nu,*/ K, G, rho;
    std::vector<std::vector<EP_FLOAT>> Gav;                    // for tauXY
    EP_FLOAT cohesion;                                         // yield strength

    // numeric
    EP_FLOAT dT;              // time step
    EP_FLOAT dampX, dampY;    // damping factor

    // variables
    std::vector<std::vector<EP_FLOAT>> Ux, Uy;                       // displacement
    std::vector<std::vector<EP_FLOAT>> divU;                         // displacement divergence
    std::vector<std::vector<EP_FLOAT>> Vx, Vy;                       // velocity
    std::vector<std::vector<EP_FLOAT>> Pinit, P;                     // hydrostatic stress (ball part of tensor)
    std::vector<std::vector<EP_FLOAT>> tauXX, tauYY, tauXY;          // deviatoric stress
    std::vector<std::vector<EP_FLOAT>> tauXXav, tauYYav, tauXYav;    // deviatoric stress for plasticity criteria

    // effective stress
    std::vector<std::array<EP_FLOAT, 3>> Sigma;
    std::vector<EP_FLOAT> deltaP;
    std::vector<EP_FLOAT> tauInfty;
    std::vector<EP_FLOAT> Keff;
    std::vector<std::array<EP_FLOAT, 3>> Geff;

    void ComputeSigma(const EP_FLOAT loadValue, const std::array<EP_FLOAT, 3>& loadType);

    SigmaCalc() = delete;
    explicit SigmaCalc(std::unique_ptr<InputData> inp_);

  private:
    void ComputeDivergence(const std::vector<std::vector<EP_FLOAT>>& Ax,
                           const std::vector<std::vector<EP_FLOAT>>& Ay,
                           std::vector<std::vector<EP_FLOAT>>& divA) const;
    void SetMaterials();
    void SetPressure(EP_FLOAT coh);
  };

} // namespace

#endif    // __SIGMA_CALC_HPP__INCLUDED__