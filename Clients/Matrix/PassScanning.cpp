//--- Aurora/Clients/Matrix/PassScanning.cpp -----------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "PassScanning.hpp"

#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Ctf.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/StemDetector.hpp"
#include "AuroraLib/StemProbe.hpp"

#include <chrono>
#include <iostream>

namespace Aurora
{

PassScanning::PassScanning(const MultisliceParameters& params,
    const std::string& matrixDir, size_t power, const std::string& outTemplate,
    const std::shared_ptr<StemDetector>& stemDetector) :
  PassFinal(params, matrixDir, power, outTemplate), mStemDetector(stemDetector)
{
  assert(mPower >= 1 && "Invalid power specified.");
}

/// create an real space image for the electron probe
void PassScanning::run()
{
  const int64_t cols = mParams.cols();
  const int64_t rows = mParams.rows();
  const Real deltaX = mParams.deltaX();
  const Real deltaY = mParams.deltaY();

  // load the necessary matrices
  MatrixNC productMatrix;
  loadMatrix(matrixName( "ProductMatrix-", mPower), productMatrix);

  MatrixNC t11;
  loadMatrix(matrixName("T11-", mPower), t11);

  MatrixNC t12;
  loadMatrix(matrixName("T12-", mPower), t12);

  MatrixNC u;
  VectorNR eigenValuesSqr;
  loadEVDecomposition(matrixName("A-", 0), eigenValuesSqr, u);

  // We need the eigenvalues for the calculation of the current
  VectorNC eigenValues;
  sqrtEigenValues(eigenValuesSqr, eigenValues);

  // the mask descibing the STEM detector in Fourier space
  CBuffer2D stemDetectorMask = mStemDetector->detectorMask(R(1.0) / (cols * rows));

  CBuffer2D tempBuffer(cols, rows);

  // the incommming wave
  CBuffer2D waveIn(cols, rows);
  auto ctf = CtfCoherent::fromCommandLine(mParams.energy());
  StemProbe probe(ctf, mParams, waveIn);

  // the SEM and STEM signal
  CBuffer2D semSignal(cols, rows);
  CBuffer2D stemSignal(cols, rows);

  for (int64_t y = 0; y < rows; ++y)
  {
    if (verbosity >= 2)
      std::cout << "Line " << y << "... ";

    for (int64_t x = 0; x < cols; ++x)
    {
      // calculate the probe
      probe.generate(x * deltaX, y * deltaY);

      // This formula and the even conditions in the next line are used to
      // undersample the data four times. This makes the produced data easier to
      // handle.
      size_t index = (y / 2) * (cols / 2) + (x / 2);
      if (signal(Signal::probe) && isEven(x) && isEven(y))
        saveWave(waveIn, mParams, "PsiIn", index);

      VectorNC psiIn;
      bufferToVector(waveIn, psiIn);

      // the SEM image
      VectorNC f0 = u.adjoint() * psiIn;
      VectorNC b0 = -productMatrix * f0;

      if (signal(Signal::psiBack) && isEven(x) && isEven(y))
      {
        // go to real space
        vectorToBuffer(u * b0, tempBuffer);
        saveWave(tempBuffer, mParams, "psiBack", index);
      }

      if (signal(Signal::sem))
      {
        // calculate the current in backward direction
        // we stay in the eigen base
        Real semResult = 0.0;
        #pragma omp parallel for reduction (+:semResult)
        for (int64_t i = 0; i < b0.size(); ++i)
          semResult += absSqr(b0[i]) * eigenValues[i].real();
        semSignal.pixel(x, y) = semResult;
      }

      // skip expensive calculations if they are not requested
      if (!signal(Signal::stem) || !signal(Signal::psiEx))
        continue;

      // the exit-wave in real space
      VectorNC psiEx = u * (t11 * f0 + t12 * b0);
      vectorToBuffer(psiEx, tempBuffer);

      if (signal(Signal::psiEx) && isEven(x) && isEven(y))
        saveWave(tempBuffer, mParams, "psiEx", index);

      if (!signal(Signal::stem))
        continue;

      // go to Fourier space
      CBuffer2D bufferFft = tempBuffer.fft(Direction::forward);

      // apply the detector geometry
      bufferFft.assign(stemDetectorMask * bufferFft);

      // go back to the eigen base
      tempBuffer = bufferFft.fft(Direction::backward);
      bufferToVector(tempBuffer, psiEx);
      VectorNC fN = u.adjoint() * psiEx;

      // calculate the current in forward direction
      Real stemResult = 0.0;
      #pragma omp parallel for reduction (+:stemResult)
      for (int64_t i = 0; i < fN.size(); ++i)
        stemResult += absSqr(fN[i]) * eigenValues[i].real();

      stemSignal.pixel(x, y) = stemResult;
    }

    if (verbosity >= 2)
      std::cout << " finished." << std::endl;
  }

  if (signal(Signal::sem))
    saveWave(semSignal, mParams, "SEM", mPower);

  if (signal(Signal::stem))
    saveWave(stemSignal, mParams, "STEM", mPower);
}

PassScanningIterative::PassScanningIterative(const MultisliceParameters& params,
    const std::string& matrixDir, size_t power, const std::string& outTemplate,
    const std::shared_ptr<StemDetector>& stemDetector, int64_t scanCols,
    int64_t scanRows, size_t minIterations, size_t maxIterations) :
  PassScanning(params, matrixDir, power, outTemplate, stemDetector),
  mScanCols(scanCols), mScanRows(scanRows), mMinIterations(minIterations),
  mMaxIterations(maxIterations)
{ }

void PassScanningIterative::run()
{
  using namespace std::chrono;
  auto startTime = high_resolution_clock::now();

  const int64_t rows = mParams.rows();
  const int64_t cols = mParams.cols();
  const int64_t numPixels = mParams.numPixels();
  const int64_t scanPixels = mScanCols * mScanRows;

  MatrixNC probeMatrix(numPixels, scanPixels);

  Buffer2DInfo scanInfo(mScanCols, mScanRows,
    mParams.width() * mScanCols / cols, mParams.height() * mScanRows / rows,
    mParams.energy());

  // the incommming wave
  CBuffer2D waveIn(cols, rows);
  auto ctf = CtfCoherent::fromCommandLine(mParams.energy());
  StemProbe probe(ctf, mParams, waveIn);

  // store the incomming wave as columns of the matrix probeMatrix
  for (int64_t y = 0; y < mScanRows; ++y)
  {
    if (verbosity >= 2)
      std::cout << "Line " << y << "... ";

    for (int64_t x = 0; x < mScanCols; ++x)
    {
      // calculate the probe
      probe.generate(x * mParams.deltaX(), y * mParams.deltaY());

      // This formula and the even conditions in the next line are used to
      // undersample the data four times. This makes the produced data easier to
      // handle.
      size_t index = (y / 2) * (mScanCols / 2) + (x / 2);
      if (signal(Signal::probe) && isEven(x) && isEven(y))
        saveWave(waveIn, mParams, "psiIn", index);

      // copy the probe into the matrix
      // for every pixel
      #pragma omp parallel for
      for (int64_t j = 0; j < rows; ++j)
        for (int64_t i = 0; i < cols; ++i)
          probeMatrix(j * cols + i, y * mScanCols + x) = waveIn.pixel(i, j);
    }
  }

  if (verbosity >= 1)
  {
    std::cout << "Probe matrix generated. Elapsed time: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }

  // calculate the probes in the eigenbase
  VectorNC eigenValues;

  {
    VectorNR eigenValuesSqr;
    MatrixNC u;
    loadEVDecomposition(matrixName("A-", 0), eigenValuesSqr, u);
    probeMatrix = u.adjoint() * probeMatrix;

    sqrtEigenValues(eigenValuesSqr, eigenValues);
  }

  // transfer the wave through the sample
  {
    MatrixNC t11;
    loadMatrix(matrixDir() + "/T-power" + toString(mPower, 4) + "-11.mat", t11);
    probeMatrix = t11 * probeMatrix;
  }

  if (verbosity >= 1)
  {
    std::cout << "Probe matrix transferred. Elapsed time: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }

  // initial values: no backscattered wave
  MatrixNC b0 = MatrixNC::Zero(numPixels, scanPixels);

  MatrixNC fN = probeMatrix;

  // start a block here so matrices will be automatically unloaded when we are
  // finished
  {
    // the matrices are needed during iterative solving
    MatrixNC t12;
    loadMatrix(matrixDir() + "/T-power" + toString(mPower, 4) + "-12.mat", t12);

    MatrixNC s21;
    loadMatrix(matrixDir() + "/S-power" + toString(mPower, 4) + "-21.mat", s21);

    // the matrices of the previous iteration
    MatrixNC b0Prev;
    MatrixNC fNPrev;

    size_t i = 0;
    while (true)
    {
      ++i;

      b0Prev = b0;
      fNPrev = fN;

      b0.noalias() = s21 * fNPrev;
      fN.noalias() = t12 * b0Prev + probeMatrix;

      const Real distanceB = (b0 - b0Prev).squaredNorm();
      const Real distanceF = (fN - fNPrev).squaredNorm();

      if (verbosity >= 1)
      {
        std::cout << "Iteration " << i << ":\n"
                     "- DistanceB = " << distanceB << "\n"
                     "- DistanceF = " << distanceF << "\n"
                     "Elapsed time: "
                  << duration(high_resolution_clock::now() - startTime)
                  << std::endl;
      }

      if ((i > mMinIterations) &&
          (distanceF < R(1e-12)) && (distanceB < R(1e-12)))
        break;

      if ((i > mMaxIterations) ||
          std::isnan(distanceF) || std::isnan(distanceB))
      {
        std::cout << "Divergent!\n";
        break;
      }
    }
  }


  if (verbosity >= 1)
  {
    std::cout << "Iteration convergent. Elapsed time: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }

  if (signal(Signal::sem))
  {
    CBuffer2D semSignal(mScanCols, mScanRows);

    #pragma omp parallel for
    for (int64_t y = 0; y < mScanRows; ++y)
    {
      for (int64_t x = 0; x < mScanCols; ++x)
      {
        Real semResult = 0.0;
        for (int64_t i = 0; i < numPixels; ++i)
          semResult += absSqr(b0(i, y * mScanCols + x)) * eigenValues[i].real();
        semSignal.pixel(x, y) = semResult;
      }
    }

    if (verbosity >= 1)
    {
      std::cout << "SEM signal calculated. Elapsed time: "
                << duration(high_resolution_clock::now() - startTime)
                << std::endl;
    }

    saveWave(semSignal, scanInfo, "SEM", mPower);
  }

  if (signal(Signal::psiBack))
  {
    MatrixNC u;
    if (signal(Signal::psiBack))
      loadEigenvectors(matrixName("A-", 0), u);

    MatrixNC matrixBack = u * b0;

    #pragma omp parallel for
    for (int64_t y = 0; y < mScanRows; y += 2)
    {
      for (int64_t x = 0; x < mScanCols; x += 2)
      {
        // extract the wave from the columns
        CBuffer2D psiBack(cols, rows);
        for (int64_t j = 0; j < rows; ++j)
          for (int64_t i = 0; i < cols; ++i)
            psiBack.pixel(i, j) = matrixBack(j * cols + i, y * mScanCols + x);

        size_t index = (y / 2) * (mScanCols / 2) + (x / 2);
        saveWave(psiBack, mParams, "psiBack", index);
      }
    }
  }

  // skip expensive calculations if they are not requested
  if (!signal(Signal::stem) && !signal(Signal::psiEx))
    return;

  // the mask descibing the STEM detector in Fourier space
  CBuffer2D stemDetectorMask = mStemDetector->detectorMask(R(1.0) / numPixels);

  MatrixNC u;
  loadEigenvectors(matrixName("A-", 0), u);

  // the exit-waves in real space
  MatrixNC matrixEx = u * fN;

  CBuffer2D stemSignal(mScanCols, mScanRows);

  #pragma omp parallel for
  for (int64_t y = 0; y < mScanRows; ++y)
  {
    for (int64_t x = 0; x < mScanCols; ++x)
    {
      // extract the exit wave from the columns
      CBuffer2D psiEx(cols, rows);
      for (int64_t j = 0; j < rows; ++j)
        for (int64_t i = 0; i < cols; ++i)
          psiEx.pixel(i, j) = matrixEx(j * cols + i, y * mScanCols + x);

      if (signal(Signal::psiEx) && isEven(x) && isEven(y))
      {
        size_t index = (y / 2) * (mScanCols / 2) + (x / 2);
        saveWave(psiEx, mParams, "psiEx", index);
      }

      // go to Fourier space
      CBuffer2D psiExFft = psiEx.fft(Direction::forward);

      // apply the detector geometry
      CBuffer2D psiDetectorFft(cols, rows);
      psiDetectorFft.assign(stemDetectorMask * psiExFft);

      // go back to the eigen base
      CBuffer2D psiDetector = psiDetectorFft.fft(Direction::backward);
      VectorNC vecTemp;
      bufferToVector(psiDetector, vecTemp);
      VectorNC vecDetector = u.adjoint() * vecTemp;

      // calculate the current in forward direction
      Real stemResult = 0.0;
      for (int64_t i = 0; i < numPixels; ++i)
        stemResult += absSqr(vecDetector[i]) * eigenValues[i].real();

      stemSignal.pixel(x, y) = stemResult;
    }
  }

  if (signal(Signal::stem))
    saveWave(stemSignal, scanInfo, "STEM", mPower);

  if (verbosity >= 1)
  {
    std::cout << "STEM signal calculated. Elapsed time: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }
}

} // namespace Aurora
