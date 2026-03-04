#ifndef VAKYUME_POTENTIALENERGYANDCONSERVATIONOFENERGY_HPP
#define VAKYUME_POTENTIALENERGYANDCONSERVATIONOFENERGY_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double>
PotentialEnergyAndConservationOfEnergy_eqn_8_6__K_A(double K_B, double W_net) {
  std::vector<double> result;
  double K_A = (K_B - W_net);
  result.push_back(K_A);
  return result;
}

std::vector<double>
PotentialEnergyAndConservationOfEnergy_eqn_8_6__K_B(double K_A, double W_net) {
  std::vector<double> result;
  double K_B = (K_A + W_net);
  result.push_back(K_B);
  return result;
}

std::vector<double>
PotentialEnergyAndConservationOfEnergy_eqn_8_6__W_net(double K_A, double K_B) {
  std::vector<double> result;
  double W_net = ((-K_A) + K_B);
  result.push_back(W_net);
  return result;
}

std::vector<double>
PotentialEnergyAndConservationOfEnergy_eqn_8_9__A(double B, double U,
                                                  double W) {
  throw std::runtime_error("PotentialEnergyAndConservationOfEnergy_eqn_8_9__A: "
                           "requires numerical solver (not transpilable)");
}

std::vector<double>
PotentialEnergyAndConservationOfEnergy_eqn_8_9__B(double A, double U,
                                                  double W) {
  throw std::runtime_error("PotentialEnergyAndConservationOfEnergy_eqn_8_9__B: "
                           "requires numerical solver (not transpilable)");
}

std::vector<double>
PotentialEnergyAndConservationOfEnergy_eqn_8_9__K(double L) {
  std::vector<double> result;
  return {(std::sqrt((L + 1.0)) - 1.0)};
}

std::vector<double>
PotentialEnergyAndConservationOfEnergy_eqn_8_9__L(double g, double m,
                                                  double v_x, double x) {
  std::vector<double> result;
  double L = ((0.5 * m) * (((2.0 * g) * x) + std::pow(v_x, 2.0)));
  result.push_back(L);
  return result;
}

std::vector<double>
PotentialEnergyAndConservationOfEnergy_eqn_8_9__U(double A, double B,
                                                  double W) {
  throw std::runtime_error("PotentialEnergyAndConservationOfEnergy_eqn_8_9__U: "
                           "requires numerical solver (not transpilable)");
}

std::vector<double>
PotentialEnergyAndConservationOfEnergy_eqn_8_9__W(double A, double B,
                                                  double U) {
  throw std::runtime_error("PotentialEnergyAndConservationOfEnergy_eqn_8_9__W: "
                           "requires numerical solver (not transpilable)");
}

std::vector<double>
PotentialEnergyAndConservationOfEnergy_eqn_8_9__W_C(double W_NC, double W_net) {
  std::vector<double> result;
  double W_C = ((-W_NC) + W_net);
  result.push_back(W_C);
  return result;
}

std::vector<double>
PotentialEnergyAndConservationOfEnergy_eqn_8_9__W_NC(double W_C, double W_net) {
  std::vector<double> result;
  double W_NC = ((-W_C) + W_net);
  result.push_back(W_NC);
  return result;
}

std::vector<double>
PotentialEnergyAndConservationOfEnergy_eqn_8_9__W_net(double W_C, double W_NC) {
  std::vector<double> result;
  double W_net = (W_C + W_NC);
  result.push_back(W_net);
  return result;
}

std::vector<double>
PotentialEnergyAndConservationOfEnergy_eqn_8_9__g(double L, double m,
                                                  double v_x, double x) {
  std::vector<double> result;
  double g = ((L - ((0.5 * m) * std::pow(v_x, 2.0))) / (m * x));
  result.push_back(g);
  return result;
}

std::vector<double>
PotentialEnergyAndConservationOfEnergy_eqn_8_9__m(double L, double g,
                                                  double v_x, double x) {
  std::vector<double> result;
  double m = ((2.0 * L) / (((2.0 * g) * x) + std::pow(v_x, 2.0)));
  result.push_back(m);
  return result;
}

std::vector<double>
PotentialEnergyAndConservationOfEnergy_eqn_8_9__v_x(double L, double g,
                                                    double m, double x) {
  std::vector<double> result;
  double v_x = (-std::sqrt((((2.0 * L) / m) - ((2.0 * g) * x))));
  result.push_back(v_x);
  v_x = std::sqrt((((2.0 * L) / m) - ((2.0 * g) * x)));
  result.push_back(v_x);
  return result;
}

std::vector<double>
PotentialEnergyAndConservationOfEnergy_eqn_8_9__x(double L, double g, double m,
                                                  double v_x) {
  std::vector<double> result;
  double x = ((L - ((0.5 * m) * std::pow(v_x, 2.0))) / (g * m));
  result.push_back(x);
  return result;
}

#endif // VAKYUME_POTENTIALENERGYANDCONSERVATIONOFENERGY_HPP
