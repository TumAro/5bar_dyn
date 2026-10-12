#include "dynamics.hpp"
#include <iostream>

int main() {
  // parameters: m1, m2, I1, I2, l, r, b
  Eigen::Matrix<double, 7, 1> p;
  p << 1.0, 0.5, 0.01, 0.005, 0.3, 0.3, 0.2;

  // state: th1, th2, dth1, dth2
  Eigen::Vector4d x;
  x << 1.9, 1.2, 0.5, -0.3;

  auto [phi1, phi2] = Phi(x(0), x(1), p(6), p(4), p(5));

  Eigen::Matrix2d M = Mass(p(0), p(1), p(2), p(3), p(4), p(5), x(0), x(1), phi1, phi2);
  Eigen::Matrix2d C = Coriolis(p(0), p(1), p(2), p(3), p(4), p(5), x(0), x(1), x(2), x(3), phi1, phi2);
  Eigen::Vector2d G = Gravity(p(0), p(1), p(4), p(5), x(0), x(1), phi1, phi2);

  std::cout << "phi1 " << phi1 << " phi2 " << phi2 << "\n";
  std::cout << "M =\n" << M << "\n";
  std::cout << "C =\n" << C << "\n";
  std::cout << "G =\n" << G << "\n";
}
