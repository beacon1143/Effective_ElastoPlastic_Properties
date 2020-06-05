#ifndef __SIGMA_CALC_HPP__INCLUDED__
#define __SIGMA_CALC_HPP__INCLUDED__

#include <array>
#include <memory>
#include "InputData.hpp"

namespace EFF_PROPS {

  class SigmaCalc {
  public:
    std::unique_ptr<InputData> inp;
    /*double loadValue;
    std::array<int, 3> loadType;*/

    // coordinates
    std::vector<std::vector<double>> x, y;
    std::vector<std::vector<double>> xUx, yUx;
    std::vector<std::vector<double>> xUy, yUy;

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
    explicit SigmaCalc(std::unique_ptr<InputData> inp_/*, const double loadValue_, const std::array<int, 3>& loadType_*/);

  private:
    void ComputeDivergence(const std::vector<std::vector<double>>& Ax,
                           const std::vector<std::vector<double>>& Ay,
                           std::vector<std::vector<double>>& divA);
  };

} // namespace

#endif    // __SIGMA_CALC_HPP__INCLUDED__