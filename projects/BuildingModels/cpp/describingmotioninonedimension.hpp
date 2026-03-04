#ifndef VAKYUME_DESCRIBINGMOTIONINONEDIMENSION_HPP
#define VAKYUME_DESCRIBINGMOTIONINONEDIMENSION_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double>
DescribingMotionInOneDimension_eqn_3_1__t(double v_x, double x, double x_0) {
  std::vector<double> result;
  double t = ((x - x_0) / v_x);
  result.push_back(t);
  return result;
}

std::vector<double>
DescribingMotionInOneDimension_eqn_3_1__v_x(double t, double x, double x_0) {
  std::vector<double> result;
  double v_x = ((x - x_0) / t);
  result.push_back(v_x);
  return result;
}

std::vector<double>
DescribingMotionInOneDimension_eqn_3_1__x(double t, double v_x, double x_0) {
  std::vector<double> result;
  double x = ((t * v_x) + x_0);
  result.push_back(x);
  return result;
}

std::vector<double>
DescribingMotionInOneDimension_eqn_3_1__x_0(double t, double v_x, double x) {
  std::vector<double> result;
  double x_0 = (((-t) * v_x) + x);
  result.push_back(x_0);
  return result;
}

std::vector<double> DescribingMotionInOneDimension_eqn_3_14__a(double a_A) {
  std::vector<double> result;
  double a = a_A;
  result.push_back(a);
  return result;
}

std::vector<double> DescribingMotionInOneDimension_eqn_3_14__a_A(double a) {
  std::vector<double> result;
  double a_A = a;
  result.push_back(a_A);
  return result;
}

std::vector<double>
DescribingMotionInOneDimension_eqn_3_15__t(double v, double v_A, double v_B) {
  std::vector<double> result;
  double t = (v / (v_A + v_B));
  result.push_back(t);
  return result;
}

std::vector<double>
DescribingMotionInOneDimension_eqn_3_15__v(double t, double v_A, double v_B) {
  std::vector<double> result;
  double v = (t * (v_A + v_B));
  result.push_back(v);
  return result;
}

std::vector<double>
DescribingMotionInOneDimension_eqn_3_15__v_A(double t, double v, double v_B) {
  std::vector<double> result;
  double v_A = ((-v_B) + (v / t));
  result.push_back(v_A);
  return result;
}

std::vector<double>
DescribingMotionInOneDimension_eqn_3_15__v_B(double t, double v, double v_A) {
  std::vector<double> result;
  double v_B = ((-v_A) + (v / t));
  result.push_back(v_B);
  return result;
}

std::vector<double>
DescribingMotionInOneDimension_eqn_3_2__ax(double t, double v, double v_0x) {
  std::vector<double> result;
  double ax = ((v - v_0x) / t);
  result.push_back(ax);
  return result;
}

std::vector<double>
DescribingMotionInOneDimension_eqn_3_2__t(double ax, double v, double v_0x) {
  std::vector<double> result;
  double t = ((v - v_0x) / ax);
  result.push_back(t);
  return result;
}

std::vector<double>
DescribingMotionInOneDimension_eqn_3_2__v(double ax, double t, double v_0x) {
  std::vector<double> result;
  double v = ((ax * t) + v_0x);
  result.push_back(v);
  return result;
}

std::vector<double>
DescribingMotionInOneDimension_eqn_3_2__v_0x(double ax, double t, double v) {
  std::vector<double> result;
  double v_0x = (((-ax) * t) + v);
  result.push_back(v_0x);
  return result;
}

std::vector<double> DescribingMotionInOneDimension_eqn_3_3__ax(double t,
                                                               double v_0x,
                                                               double x,
                                                               double x_0) {
  std::vector<double> result;
  double ax = ((2.0 * ((((-t) * v_0x) + x) - x_0)) / std::pow(t, 2.0));
  result.push_back(ax);
  return result;
}

std::vector<double> DescribingMotionInOneDimension_eqn_3_3__t(double ax,
                                                              double v_0x,
                                                              double x,
                                                              double x_0) {
  std::vector<double> result;
  double t = (((-v_0x) - std::sqrt(((((2.0 * ax) * x) - ((2.0 * ax) * x_0)) +
                                    std::pow(v_0x, 2.0)))) /
              ax);
  result.push_back(t);
  t = (((-v_0x) + std::sqrt(((((2.0 * ax) * x) - ((2.0 * ax) * x_0)) +
                             std::pow(v_0x, 2.0)))) /
       ax);
  result.push_back(t);
  return result;
}

std::vector<double> DescribingMotionInOneDimension_eqn_3_3__v_0x(double ax,
                                                                 double t,
                                                                 double x,
                                                                 double x_0) {
  std::vector<double> result;
  double v_0x = ((((((-0.5) * ax) * std::pow(t, 2.0)) + x) - x_0) / t);
  result.push_back(v_0x);
  return result;
}

std::vector<double> DescribingMotionInOneDimension_eqn_3_3__x(double ax,
                                                              double t,
                                                              double v_0x,
                                                              double x_0) {
  std::vector<double> result;
  double x = ((((0.5 * ax) * std::pow(t, 2.0)) + (t * v_0x)) + x_0);
  result.push_back(x);
  return result;
}

std::vector<double> DescribingMotionInOneDimension_eqn_3_3__x_0(double ax,
                                                                double t,
                                                                double v_0x,
                                                                double x) {
  std::vector<double> result;
  double x_0 = (((((-0.5) * ax) * std::pow(t, 2.0)) - (t * v_0x)) + x);
  result.push_back(x_0);
  return result;
}

std::vector<double>
DescribingMotionInOneDimension_eqn_3_4__a(double t, double v, double v_0) {
  std::vector<double> result;
  double a = ((v - v_0) / t);
  result.push_back(a);
  return result;
}

std::vector<double>
DescribingMotionInOneDimension_eqn_3_4__t(double a, double v, double v_0) {
  std::vector<double> result;
  double t = ((v - v_0) / a);
  result.push_back(t);
  return result;
}

std::vector<double>
DescribingMotionInOneDimension_eqn_3_4__v(double a, double t, double v_0) {
  std::vector<double> result;
  double v = ((a * t) + v_0);
  result.push_back(v);
  return result;
}

std::vector<double>
DescribingMotionInOneDimension_eqn_3_4__v_0(double a, double t, double v) {
  std::vector<double> result;
  double v_0 = (((-a) * t) + v);
  result.push_back(v_0);
  return result;
}

std::vector<double>
DescribingMotionInOneDimension_eqn_3_6__t(double v_0, double x, double x_0) {
  std::vector<double> result;
  double t = ((x - x_0) / v_0);
  result.push_back(t);
  return result;
}

std::vector<double>
DescribingMotionInOneDimension_eqn_3_6__v_0(double t, double x, double x_0) {
  std::vector<double> result;
  double v_0 = ((x - x_0) / t);
  result.push_back(v_0);
  return result;
}

std::vector<double>
DescribingMotionInOneDimension_eqn_3_6__x(double t, double v_0, double x_0) {
  std::vector<double> result;
  double x = ((t * v_0) + x_0);
  result.push_back(x);
  return result;
}

std::vector<double>
DescribingMotionInOneDimension_eqn_3_6__x_0(double t, double v_0, double x) {
  std::vector<double> result;
  double x_0 = (((-t) * v_0) + x);
  result.push_back(x_0);
  return result;
}

#endif // VAKYUME_DESCRIBINGMOTIONINONEDIMENSION_HPP
