//--- Aurora/Clients/Matrix/PassBandlimit.cpp ----------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "PassBandlimit.hpp"

#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/Lapack.hpp"
#include "AuroraLib/PhaseShift.hpp"
#include "AuroraLib/MultisliceParameters.hpp"

#include <chrono>
#include <iostream>

namespace Aurora
{

void PassBandlimit::run()
{
  using namespace std::chrono;
  auto startTime = high_resolution_clock::now();

  MatrixNC bandlimit;
  calcBandlimitMatrix(mParams, bandlimit);
  saveMatrix(bandlimit, matrixDir() + "/bandlimit.mat");

  if (verbosity >= 1)
  {
    std::cout << "Bandlimit matrix calculated: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }
}

void calcBandlimitMatrix(const MultisliceParameters& params,
                         MatrixNC& matrix)
{
  using namespace std::chrono;
  auto startTime = high_resolution_clock::now();

  if (verbosity >= 1)
    std::cout << "- calcBandlimitMatrix started." << std::endl;

  const int64_t cols = params.cols();
  const int64_t rows = params.rows();
  const int64_t numPixels = params.numPixels();

  const Real deltaKX = params.deltaKX();
  const Real deltaKY = params.deltaKY();

  // inverse Fourier trafo
  MatrixNC invFourierTrafoMatrix;
  calcFourierTrafoMatrix(params, invFourierTrafoMatrix, -1.0);
  if (verbosity >= 1)
    std::cout << "   * inverse Fourier trafo matrix calculated." << std::endl;

  DiagonalMatrixNC diagMatrix(numPixels);

  const Real kXNyquistSqr = powerOf<2>(params.kXNyquist());
  const Real kYNyquistSqr = powerOf<2>(params.kYNyquist());

  if (Bandlimit::Smooth == params.bandlimit())
  {
    // We are using a fermi function with mu = 2/3 - alpha and k_B T = 0.015
    const Real alpha = R(0.05);
    const Real mu = R(2.0) / R(3.0) - alpha;
    // Inverse of k_B T
    const Real factor = R(1.0) / R(0.015);

    #pragma omp parallel for
    for (int64_t j = 0; j < rows; ++j)
    {
      // coordinates in reciprocal space
      const Real kY = deltaKY * (j > rows / 2 ? j - rows : j);
      const Real kYSqr = kY * kY;

      for (int64_t i = 0; i < cols; ++i)
      {
        const Real kX = deltaKX * (i > cols / 2 ? i - cols : i);
        const Real kXSqr = kX * kX;

        // normalized spatial frequency
        const Real r = std::sqrt(kXSqr / kXNyquistSqr + kYSqr / kYNyquistSqr);

        diagMatrix.diagonal()[j * cols + i] =
          R(1.0) / ((R(1.0) + std::exp((r - mu) * factor)) * numPixels);
      }
    }
  }
  else if (Bandlimit::Sharp == params.bandlimit())
  {
    #pragma omp parallel for
    for (int64_t j = 0; j < rows; ++j)
    {
      // coordinates in reciprocal space
      const Real kY = deltaKY * (j > rows / 2 ? j - rows : j);
      const Real kYSqr = kY * kY;

      for (int64_t i = 0; i < cols; ++i)
      {
        const Real kX = deltaKX * (i > cols / 2 ? i - cols : i);
        const Real kXSqr = kX * kX;

        // normalized spatial frequency
        const Real radiusSqr = kXSqr / kXNyquistSqr + kYSqr / kYNyquistSqr;

        if (radiusSqr > R(4.0) / R(9.0))
          diagMatrix.diagonal()[j * cols + i] = 0;
        else
          diagMatrix.diagonal()[j * cols + i] = R(1.0) / numPixels;
      }
    }
  }
  else
  {
    AURORA_THROW(EInvalidParameter, "Bandlimit disbaled");
  }

  invFourierTrafoMatrix *= diagMatrix;

  if (verbosity >= 1)
    std::cout << "   * Diagonal matrix applied." << std::endl;

  MatrixNC fourierTrafoMatrix;
  calcFourierTrafoMatrix(params, fourierTrafoMatrix, +1.0);
  if (verbosity >= 1)
    std::cout << "   * Fourier trafo matrix calculated." << std::endl;
  matrix.noalias() = invFourierTrafoMatrix * fourierTrafoMatrix;

  if (verbosity >= 1)
  {
    std::cout << "   * Fourier trafo matrix applied.\n"
                 " - Time for calcBandlimitMatrix: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }
}

} // namespace Aurora
