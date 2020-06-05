#include <iostream>
#include "SigmaCalc.hpp"

using namespace EFF_PROPS;

int main() {
  try {
    auto input = std::make_unique<InputData>(10.0, 10.0, 50, 50, 50000, 0.5);
    SigmaCalc sigma_calc(std::move(input));
    const double load_value = 0.002;
    sigma_calc.ComputeSigma(load_value, {1, 0, 0});

    const double C_1111 = sigma_calc.Sigma[0] / load_value;
    const double C_1122 = sigma_calc.Sigma[1] / load_value;
    const double C_1112 = sigma_calc.Sigma[2] / load_value;
    std::cout << "C_1111 = " << C_1111 << std::endl;
    std::cout << "C_1122 = " << C_1122 << std::endl;
    std::cout << "C_1112 = " << C_1112 << std::endl;
    return 0;
  }
  catch (const std::exception& e) {
    std::cerr << e.what();
    return 1;
  }
}