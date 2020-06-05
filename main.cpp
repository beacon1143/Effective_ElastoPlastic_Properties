#include <iostream>
#include "SigmaCalc.hpp"

int main() {
  try {
    auto input = std::make_unique<InputData>(10.0, 10.0, 20, 20, 10000, 0.125);
    SigmaCalc sigma_calc(std::move(input));
    const double load_value = 0.002;
    sigma_calc.ComputeSigma(load_value, {1, 0, 0});
    const double C_1111 = sigma_calc.Sigma[0] / load_value;
    std::cout << "C_1111 = " << C_1111 << std::endl;
    return 0;
  }
  catch (const std::exception& e) {
    std::cerr << e.what();
    return 1;
  }
}