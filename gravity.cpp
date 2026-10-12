#include "dynamics.hpp"

Eigen::Vector2d Gravity(double m1, double m2, double l, double r, double th1, double th2, double phi1, double phi2) {
  using std::cos;
  using std::sin;
  using std::pow;

  double sphi12 = sin(phi1 - phi2);

  Eigen::Vector2d G;
  G << (g*l*m1*cos(th1))/2. + g*l*m2*cos(th1) - (g*l*m2*cos(phi2)*pow(sphi12,-1)*sin(-phi1 + th1))/2. -
          l*((g*l*m1*cos(th2))/2. + g*l*m2*cos(th2))*pow(r,-1)*pow(sphi12,-1)*sin(-phi2 + th1),
       (g*m2*r*cos(phi1))/2. + (g*l*m2*cos(phi2)*pow(sphi12,-1)*sin(-phi1 + th2))/2. +
          l*((g*l*m1*cos(th2))/2. + g*l*m2*cos(th2))*pow(r,-1)*pow(sphi12,-1)*sin(-phi2 + th2);
  return G;
}
