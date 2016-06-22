//--- Aurora/AuroraLib/Math.cpp ------------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Math.hpp"

#include <limits>

namespace Aurora
{

EMath::EMath(const std::string& message, const std::string& function,
             const std::string& filename, size_t line) :
  Exception(message, function, filename, line, "EMath") { }

Real besselI0(Real x)
{
  // this function is even
  const Real absX = std::fabs(x);

  if (absX <= R(3.75))
  {
    // use 10.25.2 [DLMF]
    // the constants are 1 / (i!)^2
    const Real c[7] = {R(1.0), R(1.0), R(0.25), R(0.02777777778),
                       R(0.001736111111), R(6.944444444e-5), R(1.929012346e-6)};

    const Real t = R(0.25) * x * x;
    Real sum = c[6];
    for (int i = 5; i >= 0; --i)
      sum = sum * t + c[i];

    return sum;
  }
  else
  {
    // use 10.40.1 [DLMF]
    const Real c[9] = {R( 0.39894228), R( 0.01328592), R( 0.00225319),
                       R(-0.00157565), R( 0.00916281), R(-0.02057706),
                       R( 0.02635537), R(-0.01647633), R( 0.00392377)};

    const Real t = R(3.75) / absX;
    Real sum = c[8];
    for (int i = 7; i >= 0; --i)
      sum = sum * t + c[i];

    return std::exp(absX) * sum / std::sqrt(absX);
  }
}

Real besselI1(Real x)
{
  // this function is odd
  const Real absX = std::fabs(x);

  if (absX <= R(3.75))
  {
    // use 10.25.2 [DLMF]
    // the constants are 1 / (2 * i! * (i + 1)!)
    const Real c[7] = {R(0.5), R(0.25), R(0.0416666666667), R(3.47222222222e-3),
                       R(1.73611111111e-4), R(5.78703703704e-6),
                       R(1.3778659612e-7)};

    const Real t = R(0.25) * x * x;
    Real sum = c[6];
    for (int i = 5; i >= 0; --i)
      sum = sum * t + c[i];

    return x * sum;
  }
  else
  {
    // use 10.40.1 [DLMF]
    const Real c[7] = {R( 0.398942280401 ), R(-0.149603355151),
                       R(-0.0467510484845), R(-0.040907167424),
                       R(-0.05752570419  ), R(-0.110736980566),
                       R(-0.269921390129)};

    const Real t = 1 / absX;
    Real sum = c[6];
    for (int i = 5; i >= 0; --i)
      sum = sum * t + c[i];

    return std::copysign(std::exp(absX) * sum / std::sqrt(absX), x);
  }
}

Real besselK0(Real x)
{
  if (x < 0)
    return std::numeric_limits<Real>::quiet_NaN();

  if (0 == x)
    return std::numeric_limits<Real>::infinity();

  if (x <= 2)
  {
    // use 10.31.1 [DLMF]
    // the constants are digamma[i + 1] / (i!)^2
    const Real c[7] = {R(-0.5772156649  ), R(0.4227843351   ),
                       R( 0.2306960838  ), R(0.03489215746  ),
                       R( 0.002614787619), R(0.0001184803936),
                       R( 3.612624103e-6)};

    const Real t = R(0.25) * x * x;
    Real sum = c[6];
    for(int i = 5; i >= 0; --i)
      sum = sum * t + c[i];

    return -std::log(R(0.5) * x) * besselI0(x) + sum;
  }
  else
  {
    // use 10.40.2 [DLMF]
    const Real c[7] = {R( 1.25331414), R(-0.07832358), R( 0.02189568),
                       R(-0.01062446), R( 0.00587872), R(-0.00251540),
                       R( 0.00053208)};

    const Real t = 2 / x;
    Real sum = c[6];
    for (int i = 5; i >= 0; --i)
      sum = sum * t + c[i];

    return std::exp(-x) * sum / std::sqrt(x);
  }
}

Real besselK1(Real x)
{
  if (x < 0)
    return std::numeric_limits<Real>::quiet_NaN();

  if (0 == x)
    return std::numeric_limits<Real>::infinity();

  if (x <= 2)
  {
    // use 10.31.1 [DLMF]
    // the constants are (digamma(i + 1) + digamma(i + 2)) / (4 * i! (i + 1)!)
    const Real c[7] = {R(-3.86078324508e-2), R(0.168196083775  ),
                       R( 4.53937917402e-2), R(4.79554745983e-3),
                       R( 2.78839872992e-4), R(1.03556192232e-5),
                       R( 2.67886478523e-7)};

    const Real t = R(0.25) * x * x;
    Real sum = c[6];
    for (int i = 5; i >= 0; --i)
      sum = sum * t + c[i];

    return 1 / x + std::log(R(0.5) * x) * besselI1(x) - x * sum;
  }
  else
  {
    // use 10.40.2 [DLMF]
    const Real c[7] = {R( 1.25331413732 ), R(0.469992801493),
                       R(-0.146872750467), R(0.128513656658),
                       R(-0.180722329676), R(0.347890484626),
                       R(-0.847983056276)};

    const Real t = 1 / x;
    Real sum = c[6];
    for (int i = 5; i >= 0; --i)
      sum = sum * t + c[i];

    return std::exp(-x) * sum / std::sqrt(x);
  }
}

} // namespace Aurora
