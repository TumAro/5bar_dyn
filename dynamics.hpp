#pragma once
#include <Eigen/Dense>
#include <tuple>
#include <cmath>

inline constexpr double g = 9.81;

inline std::tuple<double, double> Phi(double th1, double th2, double b, double l, double r) {
  using std::cos;
  using std::sin;
  using std::atan;
  using std::atan2; 
  using std::pow;
  using std::sqrt;

  double phi2 = -2*atan(pow(-2*b*r - 2*l*(b - r)*cos(th1) + 2*b*l*cos(th2) - 2*l*r*cos(th2) + pow(b,2) + 2*pow(l,2) - 2*cos(th1 - th2)*pow(l,2),-1) * (sqrt(-((-2*b*l*cos(th1) + 2*b*l*cos(th2) + pow(b,2) + 2*pow(l,2) - 2*cos(th1 - th2)*pow(l,2)) * (-2*b*l*cos(th1) + 2*b*l*cos(th2) + pow(b,2) + 2*pow(l,2) - 2*cos(th1 - th2)*pow(l,2) - 4*pow(r,2)))) -2*l*r*sin(th1) + 2*l*r*sin(th2)));
  double phi1 = atan2(pow(r,-1)*(r*sin(phi2) - l*sin(th1) + l*sin(th2)),(b + r*cos(phi2) - l*cos(th1) + l*cos(th2))/r);  

  return std::make_tuple(phi1, phi2);
}

inline std::tuple<double, double> dPhi(double th1, double th2, double dth1, double dth2,
                                       double phi1, double phi2, double l, double r) {
  using std::sin;

  double sphi12 = sin(phi1 - phi2);
  double dphi1 = l*(-(dth1*sin(-phi2 + th1))/(r*sphi12) + dth2*sin(-phi2 + th2));
  double dphi2 = l*(-(dth1*sin(-phi1 + th1))/(r*sphi12) + dth2*sin(-phi1 + th2));

  return std::make_tuple(dphi1, dphi2);
}

Eigen::Matrix2d Mass(double m1, double m2, double I1, double I2, double l, double r, double th1, double th2, double phi1, double phi2);
Eigen::Matrix2d Coriolis(double m1, double m2, double I1, double I2, double l, double r, double th1, double th2, double dth1, double dth2, double phi1, double phi2);
Eigen::Vector2d Gravity(double m1, double m2, double l, double r, double th1, double th2, double phi1, double phi2);
