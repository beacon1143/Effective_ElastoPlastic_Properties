#ifndef __INPUT_DATA_HPP__INCLUDED__
#define __INPUT_DATA_HPP__INCLUDED__

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <array>

class InputData {
public:
  double sizeX, sizeY;
  int nX, nY;
  long int nTimeSteps;
  double courant;

  double dX, dY;
  double dT;
  double damp;

  // material properties
  double Kmax, Gmax, rho_max;
  std::vector<std::vector<double>> E, nu, K, G, rho;
  void SetMaterials();

  // boundary conditions
  double loadValue;
  std::array<int, 3> loadType;

  InputData() = delete;
  explicit InputData(const double size_x, const double size_y,
                     const int nX, const int nY,
                     const long int nTimeSteps, const double courant);
private:
  template <typename T>
  static void ResizeXY(std::vector<std::vector<T>>& matr, const int n1, const int n2);
};

#endif    // __INPUT_DATA_HPP__INCLUDED__