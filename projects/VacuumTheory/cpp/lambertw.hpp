#ifndef VAKYUME_LAMBERTW_HPP
#define VAKYUME_LAMBERTW_HPP

#include <cmath>

// Lambert W₀ (principal branch) via Halley iteration
inline double _lambertW(double x) {
  if (x == 0.0)
    return 0.0;
  double w = (x > 1.0) ? std::log(x) - std::log(std::log(x)) : 0.5;
  for (int i = 0; i < 50; ++i) {
    double ew = std::exp(w);
    double wew = w * ew;
    double f = wew - x;
    double fp = ew * (w + 1.0);
    double d = f / (fp - (w + 2.0) * f / (2.0 * (w + 1.0)));
    w -= d;
    if (std::abs(d) < 1e-15)
      break;
  }
  return w;
}

#endif // VAKYUME_LAMBERTW_HPP
