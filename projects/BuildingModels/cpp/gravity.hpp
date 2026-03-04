#ifndef VAKYUME_GRAVITY_HPP
#define VAKYUME_GRAVITY_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double> Gravity_eqn_9_2__G(double M, double R, double T,
                                       double pi) {
  std::vector<double> result;
  double G = (((4.0 * std::pow(T, 2.0)) * std::pow(M_PI, 2.0)) /
              (M * std::pow(R, 3.0)));
  result.push_back(G);
  return result;
}

std::vector<double> Gravity_eqn_9_2__M(double G, double R, double T,
                                       double pi) {
  std::vector<double> result;
  double M = (((4.0 * std::pow(T, 2.0)) * std::pow(M_PI, 2.0)) /
              (G * std::pow(R, 3.0)));
  result.push_back(M);
  return result;
}

std::vector<std::complex<double>> Gravity_eqn_9_2__R(double G, double M,
                                                     double T, double pi) {
  std::vector<std::complex<double>> result;
  std::complex<double> R =
      (std::pow(2.0, (2.0 / 3.0)) *
       std::pow(((std::pow(T, 2.0) * std::pow(M_PI, 2.0)) / (G * M)),
                (1.0 / 3.0)));
  result.push_back(R);
  R = ((((-std::pow(2.0, (2.0 / 3.0))) *
         std::pow(((std::pow(T, 2.0) * std::pow(M_PI, 2.0)) / (G * M)),
                  (1.0 / 3.0))) /
        2.0) -
       ((((std::pow(2.0, (2.0 / 3.0)) * std::sqrt(3.0)) *
          std::complex<double>(0.0, 1.0)) *
         std::pow(((std::pow(T, 2.0) * std::pow(M_PI, 2.0)) / (G * M)),
                  (1.0 / 3.0))) /
        2.0));
  result.push_back(R);
  R = ((((-std::pow(2.0, (2.0 / 3.0))) *
         std::pow(((std::pow(T, 2.0) * std::pow(M_PI, 2.0)) / (G * M)),
                  (1.0 / 3.0))) /
        2.0) +
       ((((std::pow(2.0, (2.0 / 3.0)) * std::sqrt(3.0)) *
          std::complex<double>(0.0, 1.0)) *
         std::pow(((std::pow(T, 2.0) * std::pow(M_PI, 2.0)) / (G * M)),
                  (1.0 / 3.0))) /
        2.0));
  result.push_back(R);
  return result;
}

std::vector<double> Gravity_eqn_9_2__T(double G, double M, double R,
                                       double pi) {
  std::vector<double> result;
  double T = ((-std::sqrt(((G * M) * std::pow(R, 3.0)))) / (2.0 * M_PI));
  result.push_back(T);
  T = (std::sqrt(((G * M) * std::pow(R, 3.0))) / (2.0 * M_PI));
  result.push_back(T);
  return result;
}

std::vector<double> Gravity_eqn_9_2__pi(double G, double M, double R,
                                        double T) {
  std::vector<double> result;
  double pi = ((-std::sqrt(((G * M) * std::pow(R, 3.0)))) / (2.0 * T));
  result.push_back(M_PI);
  pi = (std::sqrt(((G * M) * std::pow(R, 3.0))) / (2.0 * T));
  result.push_back(M_PI);
  return result;
}

std::vector<double> Gravity_eqn_9_3__G(double M, double R, double a, double m) {
  std::vector<double> result;
  double G = ((std::pow(R, 2.0) * a) / (M * m));
  result.push_back(G);
  return result;
}

std::vector<double> Gravity_eqn_9_3__M(double G, double R, double a, double m) {
  std::vector<double> result;
  double M = ((std::pow(R, 2.0) * a) / (G * m));
  result.push_back(M);
  return result;
}

std::vector<double> Gravity_eqn_9_3__R(double G, double M, double a, double m) {
  std::vector<double> result;
  double R = (-std::sqrt((((G * M) * m) / a)));
  result.push_back(R);
  R = std::sqrt((((G * M) * m) / a));
  result.push_back(R);
  return result;
}

std::vector<double> Gravity_eqn_9_3__a(double G, double M, double R, double m) {
  std::vector<double> result;
  double a = (((G * M) * m) / std::pow(R, 2.0));
  result.push_back(a);
  return result;
}

std::vector<double> Gravity_eqn_9_3__m(double G, double M, double R, double a) {
  std::vector<double> result;
  double m = ((std::pow(R, 2.0) * a) / (G * M));
  result.push_back(m);
  return result;
}

std::vector<double> Gravity_eqn_9_4__G(double M, double U, double m, double r) {
  std::vector<double> result;
  double G = (((-U) * r) / (M * m));
  result.push_back(G);
  return result;
}

std::vector<double> Gravity_eqn_9_4__M(double G, double U, double m, double r) {
  std::vector<double> result;
  double M = (((-U) * r) / (G * m));
  result.push_back(M);
  return result;
}

std::vector<double> Gravity_eqn_9_4__U(double G, double M, double m, double r) {
  std::vector<double> result;
  double U = ((((-G) * M) * m) / r);
  result.push_back(U);
  return result;
}

std::vector<double> Gravity_eqn_9_4__m(double G, double M, double U, double r) {
  std::vector<double> result;
  double m = (((-U) * r) / (G * M));
  result.push_back(m);
  return result;
}

std::vector<double> Gravity_eqn_9_4__r(double G, double M, double U, double m) {
  std::vector<double> result;
  double r = ((((-G) * M) * m) / U);
  result.push_back(r);
  return result;
}

#endif // VAKYUME_GRAVITY_HPP
