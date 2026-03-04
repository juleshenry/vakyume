#ifndef VAKYUME_DESCRIBINGMOTIONINMULTIPLEDIMENSIONS_HPP
#define VAKYUME_DESCRIBINGMOTIONINMULTIPLEDIMENSIONS_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double> DescribingMotionInMultipleDimensions_eqn_4_1__t(
    double v, double x_1, double x_2, double y_1, double y_2) {
  std::vector<double> result;
  double t = (((((-x_1) + x_2) - y_1) + y_2) / v);
  result.push_back(t);
  return result;
}

std::vector<double> DescribingMotionInMultipleDimensions_eqn_4_1__v(
    double t, double x_1, double x_2, double y_1, double y_2) {
  std::vector<double> result;
  double v = (((((-x_1) + x_2) - y_1) + y_2) / t);
  result.push_back(v);
  return result;
}

std::vector<double> DescribingMotionInMultipleDimensions_eqn_4_1__x_1(
    double t, double v, double x_2, double y_1, double y_2) {
  std::vector<double> result;
  double x_1 = (((((-t) * v) + x_2) - y_1) + y_2);
  result.push_back(x_1);
  return result;
}

std::vector<double> DescribingMotionInMultipleDimensions_eqn_4_1__x_2(
    double t, double v, double x_1, double y_1, double y_2) {
  std::vector<double> result;
  double x_2 = ((((t * v) + x_1) + y_1) - y_2);
  result.push_back(x_2);
  return result;
}

std::vector<double> DescribingMotionInMultipleDimensions_eqn_4_1__y_1(
    double t, double v, double x_1, double x_2, double y_2) {
  std::vector<double> result;
  double y_1 = (((((-t) * v) - x_1) + x_2) + y_2);
  result.push_back(y_1);
  return result;
}

std::vector<double> DescribingMotionInMultipleDimensions_eqn_4_1__y_2(
    double t, double v, double x_1, double x_2, double y_1) {
  std::vector<double> result;
  double y_2 = ((((t * v) + x_1) - x_2) + y_1);
  result.push_back(y_2);
  return result;
}

std::vector<double> DescribingMotionInMultipleDimensions_eqn_4_3__v(double vx,
                                                                    double vy,
                                                                    double x,
                                                                    double y) {
  std::vector<double> result;
  double v = ((vx * x) + (vy * y));
  result.push_back(v);
  return result;
}

std::vector<double> DescribingMotionInMultipleDimensions_eqn_4_3__vx(double v,
                                                                     double vy,
                                                                     double x,
                                                                     double y) {
  std::vector<double> result;
  double vx = ((v - (vy * y)) / x);
  result.push_back(vx);
  return result;
}

std::vector<double> DescribingMotionInMultipleDimensions_eqn_4_3__vy(double v,
                                                                     double vx,
                                                                     double x,
                                                                     double y) {
  std::vector<double> result;
  double vy = ((v - (vx * x)) / y);
  result.push_back(vy);
  return result;
}

std::vector<double> DescribingMotionInMultipleDimensions_eqn_4_3__x(double v,
                                                                    double vx,
                                                                    double vy,
                                                                    double y) {
  std::vector<double> result;
  double x = ((v - (vy * y)) / vx);
  result.push_back(x);
  return result;
}

std::vector<double> DescribingMotionInMultipleDimensions_eqn_4_3__y(double v,
                                                                    double vx,
                                                                    double vy,
                                                                    double x) {
  std::vector<double> result;
  double y = ((v - (vx * x)) / vy);
  result.push_back(y);
  return result;
}

std::vector<double> DescribingMotionInMultipleDimensions_eqn_4_4__a(double ax,
                                                                    double ay,
                                                                    double x,
                                                                    double y) {
  std::vector<double> result;
  double a = ((ax * x) + (ay * y));
  result.push_back(a);
  return result;
}

std::vector<double> DescribingMotionInMultipleDimensions_eqn_4_4__ax(double a,
                                                                     double ay,
                                                                     double x,
                                                                     double y) {
  std::vector<double> result;
  double ax = ((a - (ay * y)) / x);
  result.push_back(ax);
  return result;
}

std::vector<double> DescribingMotionInMultipleDimensions_eqn_4_4__ay(double a,
                                                                     double ax,
                                                                     double x,
                                                                     double y) {
  std::vector<double> result;
  double ay = ((a - (ax * x)) / y);
  result.push_back(ay);
  return result;
}

std::vector<double> DescribingMotionInMultipleDimensions_eqn_4_4__x(double a,
                                                                    double ax,
                                                                    double ay,
                                                                    double y) {
  std::vector<double> result;
  double x = ((a - (ay * y)) / ax);
  result.push_back(x);
  return result;
}

std::vector<double> DescribingMotionInMultipleDimensions_eqn_4_4__y(double a,
                                                                    double ax,
                                                                    double ay,
                                                                    double x) {
  std::vector<double> result;
  double y = ((a - (ax * x)) / ay);
  result.push_back(y);
  return result;
}

std::vector<double>
DescribingMotionInMultipleDimensions_eqn_4_5__a(double t, double v,
                                                double v_0) {
  std::vector<double> result;
  double a = ((v - v_0) / t);
  result.push_back(a);
  return result;
}

std::vector<double>
DescribingMotionInMultipleDimensions_eqn_4_5__t(double a, double v,
                                                double v_0) {
  std::vector<double> result;
  double t = ((v - v_0) / a);
  result.push_back(t);
  return result;
}

std::vector<double>
DescribingMotionInMultipleDimensions_eqn_4_5__v(double a, double t,
                                                double v_0) {
  std::vector<double> result;
  double v = ((a * t) + v_0);
  result.push_back(v);
  return result;
}

std::vector<double>
DescribingMotionInMultipleDimensions_eqn_4_5__v_0(double a, double t,
                                                  double v) {
  std::vector<double> result;
  double v_0 = (((-a) * t) + v);
  result.push_back(v_0);
  return result;
}

std::vector<double>
DescribingMotionInMultipleDimensions_eqn_4_6__a(double r, double r_0, double t,
                                                double v_0) {
  std::vector<double> result;
  double a = ((2.0 * ((r - r_0) - (t * v_0))) / std::pow(t, 2.0));
  result.push_back(a);
  return result;
}

std::vector<double>
DescribingMotionInMultipleDimensions_eqn_4_6__r(double a, double r_0, double t,
                                                double v_0) {
  std::vector<double> result;
  double r = ((((0.5 * a) * std::pow(t, 2.0)) + r_0) + (t * v_0));
  result.push_back(r);
  return result;
}

std::vector<double>
DescribingMotionInMultipleDimensions_eqn_4_6__r_0(double a, double r, double t,
                                                  double v_0) {
  std::vector<double> result;
  double r_0 = (((((-0.5) * a) * std::pow(t, 2.0)) + r) - (t * v_0));
  result.push_back(r_0);
  return result;
}

std::vector<double>
DescribingMotionInMultipleDimensions_eqn_4_6__t(double a, double r, double r_0,
                                                double v_0) {
  std::vector<double> result;
  double t = (((-v_0) - std::sqrt(((((2.0 * a) * r) - ((2.0 * a) * r_0)) +
                                   std::pow(v_0, 2.0)))) /
              a);
  result.push_back(t);
  t = (((-v_0) + std::sqrt(((((2.0 * a) * r) - ((2.0 * a) * r_0)) +
                            std::pow(v_0, 2.0)))) /
       a);
  result.push_back(t);
  return result;
}

std::vector<double>
DescribingMotionInMultipleDimensions_eqn_4_6__v_0(double a, double r,
                                                  double r_0, double t) {
  std::vector<double> result;
  double v_0 = ((((((-0.5) * a) * std::pow(t, 2.0)) + r) - r_0) / t);
  result.push_back(v_0);
  return result;
}

#endif // VAKYUME_DESCRIBINGMOTIONINMULTIPLEDIMENSIONS_HPP
