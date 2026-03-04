#ifndef VAKYUME_ELECTRICCIRCUITS_HPP
#define VAKYUME_ELECTRICCIRCUITS_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double> ElectricCircuits_eqn_20_3__I1(double I2, double I3) {
  std::vector<double> result;
  double I1 = (I2 + I3);
  result.push_back(I1);
  return result;
}

std::vector<double> ElectricCircuits_eqn_20_3__I2(double I1, double I3) {
  std::vector<double> result;
  double I2 = (I1 - I3);
  result.push_back(I2);
  return result;
}

std::vector<double> ElectricCircuits_eqn_20_3__I3(double I1, double I2) {
  std::vector<double> result;
  double I3 = (I1 - I2);
  result.push_back(I3);
  return result;
}

std::vector<std::complex<double>>
ElectricCircuits_eqn_20_4__I(double Reff, double Vvoltmeter) {
  std::vector<std::complex<double>> result;
  std::complex<double> I = (Vvoltmeter / Reff);
  result.push_back(std::complex<double>(0.0, 1.0));
  return result;
}

std::vector<std::complex<double>>
ElectricCircuits_eqn_20_4__Reff(double I, double Vvoltmeter) {
  std::vector<std::complex<double>> result;
  std::complex<double> Reff = (Vvoltmeter / std::complex<double>(0.0, 1.0));
  result.push_back(Reff);
  return result;
}

std::vector<std::complex<double>>
ElectricCircuits_eqn_20_4__Vvoltmeter(double I, double Reff) {
  std::vector<std::complex<double>> result;
  std::complex<double> Vvoltmeter = (std::complex<double>(0.0, 1.0) * Reff);
  result.push_back(Vvoltmeter);
  return result;
}

std::vector<double> ElectricCircuits_eqn_20_5__C(double IR, double Q,
                                                 double V) {
  std::vector<double> result;
  double C = ((-Q) / (IR - V));
  result.push_back(C);
  return result;
}

std::vector<double> ElectricCircuits_eqn_20_5__IR(double C, double Q,
                                                  double V) {
  std::vector<double> result;
  double IR = (V - (Q / C));
  result.push_back(IR);
  return result;
}

std::vector<double> ElectricCircuits_eqn_20_5__Q(double C, double IR,
                                                 double V) {
  std::vector<double> result;
  double Q = (C * ((-IR) + V));
  result.push_back(Q);
  return result;
}

std::vector<double> ElectricCircuits_eqn_20_5__V(double C, double IR,
                                                 double Q) {
  std::vector<double> result;
  double V = (IR + (Q / C));
  result.push_back(V);
  return result;
}

#endif // VAKYUME_ELECTRICCIRCUITS_HPP
