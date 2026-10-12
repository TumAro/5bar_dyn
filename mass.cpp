#include "dynamics.hpp"

Eigen::Matrix2d Mass(double m1, double m2, double I1, double I2, double l, double r, double th1, double th2, double phi1, double phi2) {
  using std::cos;
  using std::sin;
  using std::pow;

  double sphi12 = sin(phi1 - phi2);

  double m11 = (4*I1 + (m1 + 4*m2)*pow(l,2) + pow(l,2)*pow(r,-2)*pow(sphi12,-2)*((4*I2 + m2*pow(r,2))*pow(sin(-phi1 + th1),2) +
      (4*I1 + (m1 + 4*m2)*pow(l,2))*pow(sin(-phi2 + th1),2) + 4*l*m2*r*cos(-phi2 + th2)*sin(-phi1 + th1)*sin(-phi2 + th1)))/4.;

  double m12 = -0.25*(l*pow(r,-2)*(-2*m2*cos(-phi1 + th1)*pow(r,3) + l*pow(sphi12,-2)*
      (sin(-phi2 + th1)*(2*l*m2*r*cos(-phi2 + th2)*sin(-phi1 + th2) + (4*I1 + (m1 + 4*m2)*pow(l,2))*sin(-phi2 + th2)) +
        sin(-phi1 + th1)*((4*I2 + m2*pow(r,2))*sin(-phi1 + th2) + l*m2*r*sin(2*(-phi2 + th2))))));

  double m22 = I2 + (m2*pow(r,2))/4. + (pow(l,2)*pow(r,-2)*pow(sphi12,-2)*((4*I2 + m2*pow(r,2))*pow(sin(-phi1 + th2),2) + (4*I1 + (m1 + 4*m2)*pow(l,2))*pow(sin(-phi2 + th2),2) +
      2*l*m2*r*sin(-phi1 + th2)*sin(2*(-phi2 + th2))))/4.;

  Eigen::Matrix2d M;
  M << m11, m12,
       m12, m22;
  return M;
}
