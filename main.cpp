#include <iostream>
#include <chrono>
#include "SigmaCalc.hpp"

using namespace EFF_PROPS;

int main() {
  try {
    auto input = std::make_unique<InputData>(16.0, 16.0, 100, 100, 4'000'000, 0.125);
    SigmaCalc sigma_calc(std::move(input));
    const double load_value = 0.002;

    const auto start = std::chrono::system_clock::now();

    sigma_calc.ComputeSigma(load_value, {1, 0, 0});
    const double C_1111 = sigma_calc.Sigma[0] / load_value;
    const double C_1122 = sigma_calc.Sigma[1] / load_value;
    const double C_1112 = sigma_calc.Sigma[2] / load_value;

    sigma_calc.ComputeSigma(load_value, {0, 1, 0});
    const double C_2222 = sigma_calc.Sigma[1] / load_value;
    const double C_1222 = sigma_calc.Sigma[2] / load_value;

    sigma_calc.ComputeSigma(load_value, {0, 0, 1});
    const double C_1212 = sigma_calc.Sigma[2] / load_value;

    const auto end = std::chrono::system_clock::now();

    std::cout << "C_1111 = " << C_1111 << std::endl;
    std::cout << "C_1122 = " << C_1122 << std::endl;
    std::cout << "C_1112 = " << C_1112 << std::endl;
    std::cout << "C_2222 = " << C_2222 << std::endl;
    std::cout << "C_1222 = " << C_1222 << std::endl;
    std::cout << "C_1212 = " << C_1212 << std::endl;

    const int elapsed_sec = static_cast<int>( std::chrono::duration_cast<std::chrono::seconds>(end - start).count() );
    if (elapsed_sec < 60) {
      std::cout << "Calculation time is " << elapsed_sec << " sec\n";
    }
    else {
      std::cout << "Calculation time is " << elapsed_sec / 60 << " min " << elapsed_sec % 60 << " sec\n";
    }

    return 0;
  }
  catch (const std::exception& e) {
    std::cerr << e.what();
    return 1;
  }
}