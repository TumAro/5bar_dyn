#include "dynamics.hpp"

Eigen::Matrix2d Coriolis(double m1, double m2, double I1, double I2, double l, double r, double th1, double th2, double dth1, double dth2, double phi1, double phi2) {
  using std::cos;
  using std::sin;
  using std::pow;

  auto [dphi1, dphi2] = dPhi(th1, th2, dth1, dth2, phi1, phi2, l, r);
  double sphi12 = sin(phi1 - phi2);
  double cphi12 = cos(phi1 - phi2);

  double c11 = -0.5*pow(l,2)*pow(r,-2)*pow(sphi12,-2)*(-1.*sin(-phi2 + th1)*(1.*l*m2*r*cos(-phi2 + th2)*
       ((-dphi1 + dth1)*cos(-phi1 + th1) - cphi12*(dphi1 - dphi2)*pow(sphi12,-1)*sin(-phi1 + th1)) +
      0.5*(4*I1 + (m1 + 4*m2)*pow(l,2))*((-dphi2 + dth1)*cos(-phi2 + th1) - cphi12*(dphi1 - dphi2)*pow(sphi12,-1)*sin(-phi2 + th1)) +
      1.*dphi2*l*m2*r*sin(-phi1 + th1)*sin(-phi2 + th2)) + 1.*sin(-phi1 + th1)*
    (-0.5*(4*I2 + m2*pow(r,2))*((-dphi1 + dth1)*cos(-phi1 + th1) - cphi12*(dphi1 - dphi2)*pow(sphi12,-1)*sin(-phi1 + th1)) -
      1.*l*m2*r*cos(-phi2 + th2)*((-dphi2 + dth1)*cos(-phi2 + th1) - cphi12*(dphi1 - dphi2)*pow(sphi12,-1)*sin(-phi2 + th1)) + 1.*dth2*l*m2*r*sin(-phi2 + th1)*sin(-phi2 + th2)));

  double c12 = 0. + 0.5*dphi1*l*m2*r*sin(-phi1 + th1) - l*pow(r,-1)*pow(sphi12,-1)*sin(-phi1 + th1)*
  (0. - 0.5*dth2*m2*pow(l,2)*pow(sphi12,-1)*pow(sin(-phi2 + th2),2) -
    (l*pow(r,-1)*(4*I2 + m2*pow(r,2))*((dphi1 - dth2)*cos(-phi1 + th2)*pow(sphi12,-1) + cphi12*(dphi1 - dphi2)*pow(sphi12,-2)*sin(-phi1 + th2)))/4. -
    (pow(l,2)*m2*r*cos(-phi2 + th2)*pow(r,-1)*((dphi2 - dth2)*cos(-phi2 + th2)*pow(sphi12,-1) + cphi12*(dphi1 - dphi2)*pow(sphi12,-2)*sin(-phi2 + th2)))/2.) -
 l*pow(r,-1)*pow(sphi12,-1)*sin(-phi2 + th1)*(0. - (pow(l,2)*m2*r*cos(-phi2 + th2)*pow(r,-1)*
       ((dphi1 - dth2)*cos(-phi1 + th2)*pow(sphi12,-1) + cphi12*(dphi1 - dphi2)*pow(sphi12,-2)*sin(-phi1 + th2)))/2. +
    0.5*dphi2*m2*pow(l,2)*pow(sphi12,-1)*sin(-phi1 + th2)*sin(-phi2 + th2) -
    (l*(4*I1 + (m1 + 4*m2)*pow(l,2))*pow(r,-1)*((dphi2 - dth2)*cos(-phi2 + th2)*pow(sphi12,-1) + cphi12*(dphi1 - dphi2)*pow(sphi12,-2)*sin(-phi2 + th2)))/4.);

  double c21 = 0. - 0.5*dth1*l*m2*r*sin(-phi1 + th1) + l*pow(r,-1)*pow(sphi12,-1)*sin(-phi2 + th2)*
  (0. + (pow(l,2)*m2*r*cos(-phi2 + th2)*pow(r,-1)*((dphi1 - dth1)*cos(-phi1 + th1)*pow(sphi12,-1) + cphi12*(dphi1 - dphi2)*pow(sphi12,-2)*sin(-phi1 + th1)))/2. +
    (l*(4*I1 + (m1 + 4*m2)*pow(l,2))*pow(r,-1)*((dphi2 - dth1)*cos(-phi2 + th1)*pow(sphi12,-1) + cphi12*(dphi1 - dphi2)*pow(sphi12,-2)*sin(-phi2 + th1)))/4. -
    0.5*dphi2*m2*pow(l,2)*pow(sphi12,-1)*sin(-phi1 + th1)*sin(-phi2 + th2)) +
 l*pow(r,-1)*pow(sphi12,-1)*sin(-phi1 + th2)*(0. + (l*pow(r,-1)*(4*I2 + m2*pow(r,2))*
       ((dphi1 - dth1)*cos(-phi1 + th1)*pow(sphi12,-1) + cphi12*(dphi1 - dphi2)*pow(sphi12,-2)*sin(-phi1 + th1)))/4. +
    (pow(l,2)*m2*r*cos(-phi2 + th2)*pow(r,-1)*((dphi2 - dth1)*cos(-phi2 + th1)*pow(sphi12,-1) + cphi12*(dphi1 - dphi2)*pow(sphi12,-2)*sin(-phi2 + th1)))/2. +
    0.5*dth2*m2*pow(l,2)*pow(sphi12,-1)*sin(-phi2 + th1)*sin(-phi2 + th2));

  double c22 = 0. + l*pow(r,-1)*pow(sphi12,-1)*sin(-phi1 + th2)*(0. - 0.5*dth2*m2*pow(l,2)*pow(sphi12,-1)*pow(sin(-phi2 + th2),2) -
    (l*pow(r,-1)*(4*I2 + m2*pow(r,2))*((dphi1 - dth2)*cos(-phi1 + th2)*pow(sphi12,-1) + cphi12*(dphi1 - dphi2)*pow(sphi12,-2)*sin(-phi1 + th2)))/4. -
    (pow(l,2)*m2*r*cos(-phi2 + th2)*pow(r,-1)*((dphi2 - dth2)*cos(-phi2 + th2)*pow(sphi12,-1) + cphi12*(dphi1 - dphi2)*pow(sphi12,-2)*sin(-phi2 + th2)))/2.) +
 l*pow(r,-1)*pow(sphi12,-1)*sin(-phi2 + th2)*(0. - (pow(l,2)*m2*r*cos(-phi2 + th2)*pow(r,-1)*
       ((dphi1 - dth2)*cos(-phi1 + th2)*pow(sphi12,-1) + cphi12*(dphi1 - dphi2)*pow(sphi12,-2)*sin(-phi1 + th2)))/2. +
    0.5*dphi2*m2*pow(l,2)*pow(sphi12,-1)*sin(-phi1 + th2)*sin(-phi2 + th2) -
    (l*(4*I1 + (m1 + 4*m2)*pow(l,2))*pow(r,-1)*((dphi2 - dth2)*cos(-phi2 + th2)*pow(sphi12,-1) + cphi12*(dphi1 - dphi2)*pow(sphi12,-2)*sin(-phi2 + th2)))/4.);

  Eigen::Matrix2d C;
  C << c11, c12,
       c21, c22;
  return C;
}
