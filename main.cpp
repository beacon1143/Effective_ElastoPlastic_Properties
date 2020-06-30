#include <iostream>
#include <chrono>
#include "SigmaCalc.hpp"

using namespace EFF_PROPS;

int main() {
  try {
    const auto start = std::chrono::system_clock::now();

    constexpr size_t nTimeSteps = 1;
    auto input = std::make_unique<InputData>(20.0, 20.0, 200, 200, nTimeSteps, 20'000, 0.125);
    SigmaCalc sigma_calc(std::move(input));
    constexpr EP_FLOAT load_value = -0.002;

    /*std::vector<EP_FLOAT> C_1111(nTimeSteps),C_1122(nTimeSteps), C_1112(nTimeSteps),
                        C_2222(nTimeSteps), C_1222(nTimeSteps), C_1212(nTimeSteps);

    sigma_calc.ComputeSigma(load_value, {1, 0, 0});
    for (size_t tim = 0; tim < nTimeSteps; tim++) {
      C_1111[tim] = sigma_calc.Sigma[tim][0] / load_value / static_cast<EP_FLOAT>(tim + 1) * nTimeSteps;
      C_1122[tim] = sigma_calc.Sigma[tim][1] / load_value / static_cast<EP_FLOAT>(tim + 1) * nTimeSteps;
      C_1112[tim] = sigma_calc.Sigma[tim][2] / load_value / static_cast<EP_FLOAT>(tim + 1) * nTimeSteps;
    }

    sigma_calc.ComputeSigma(load_value, {0, 1, 0});
    for (size_t tim = 0; tim < nTimeSteps; tim++) {
      C_2222[tim] = sigma_calc.Sigma[tim][1] / load_value / static_cast<EP_FLOAT>(tim + 1) * nTimeSteps;
      C_1222[tim] = sigma_calc.Sigma[tim][2] / load_value / static_cast<EP_FLOAT>(tim + 1) * nTimeSteps;
    }

    sigma_calc.ComputeSigma(load_value, {0, 0, 1});
    for (size_t tim = 0; tim < nTimeSteps; tim++) {
      C_1212[tim] = sigma_calc.Sigma[tim][2] / load_value / static_cast<EP_FLOAT>(tim + 1) * nTimeSteps;
    }
    
    std::cout << std::endl;

    for (size_t tim = 0; tim < nTimeSteps; tim++) {
      std::cout << "Time step " << tim + 1 << " from " << nTimeSteps << '\n';
      std::cout << "C_1111 = " << C_1111[tim] << '\n';
      std::cout << "C_1122 = " << C_1122[tim] << '\n';
      std::cout << "C_1112 = " << C_1112[tim] << '\n';
      std::cout << "C_2222 = " << C_2222[tim] << '\n';
      std::cout << "C_1222 = " << C_1222[tim] << '\n';
      std::cout << "C_1212 = " << C_1212[tim] << '\n';
      std::cout << std::endl;
    }*/

    sigma_calc.ComputeSigma(load_value, {1.0, 1.0, 0});
    //std::cout << sigma_calc.Sigma[0][0] / load_value << '\t' << sigma_calc.Sigma[0][1] / load_value << '\t' << sigma_calc.Sigma[0][2] / load_value << std::endl;

    const auto end = std::chrono::system_clock::now();

    const int elapsed_sec = static_cast<int>( std::chrono::duration_cast<std::chrono::seconds>(end - start).count() );
    if (elapsed_sec < 60) {
      std::cout << "Calculation time is " << elapsed_sec << " sec\n";
    }
    else {
      const int elapsed_min = elapsed_sec / 60;
      if (elapsed_min < 60) {
        std::cout << "Calculation time is " << elapsed_min << " min " << elapsed_sec % 60 << " sec\n";
      }
      else {
        std::cout << "Calculation time is " << elapsed_min / 60 << " hours " << elapsed_min % 60 << " min " << elapsed_sec % 60 << " sec\n";
      }
    }

    return 0;
  }
  catch (const std::exception& e) {
    std::cerr << e.what();
    return 1;
  }
}