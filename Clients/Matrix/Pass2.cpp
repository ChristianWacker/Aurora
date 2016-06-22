//--- Aurora/Clients/Matrix/Pass2.cpp ------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "Pass2.hpp"

#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/Histogram.hpp"
#include "AuroraLib/MultisliceParameters.hpp"
#include "PassBandlimit.hpp"

#include <chrono>
#include <iostream>

namespace Aurora
{

using namespace std::chrono;

void Pass2F::run()
{
  auto startTime = high_resolution_clock::now();

  const int64_t numPixels = mParams.numPixels();
  MatrixNC transferMatrix = MatrixNC::Identity(numPixels, numPixels);

  VectorNR eigenvaluesSqr;
  VectorNC prevEigenValues;
  VectorNC curEigenvalues;
  MatrixNC prevEigenvectors;
  MatrixNC curEigenvectors;
  DiagonalMatrixNC expMatrix;

  if (mMatrixA2)
  {
    mMatrixA2->calc(mPart.startSlice(), eigenvaluesSqr, curEigenvectors);
  }
  else
  {
    loadEVDecomposition(matrixName("A-", mPart.startSlice()),
                        eigenvaluesSqr, curEigenvectors);
  }

  // calculate the square root of the eigen values
  sqrtEigenValues(eigenvaluesSqr, curEigenvalues);

  MatrixNC bandlimitMatrix;
  if (Bandlimit::Disabled != mParams.bandlimit())
    calcBandlimitMatrix(mParams, bandlimitMatrix);

  for (size_t n = mPart.startSlice() + 1; n < mPart.endSlice() + 1; ++n)
  {
    prevEigenValues  = curEigenvalues;
    prevEigenvectors = curEigenvectors;

    if (mMatrixA2)
      mMatrixA2->calc(n, eigenvaluesSqr, curEigenvectors);
    else
      loadEVDecomposition(matrixName("A-", n), eigenvaluesSqr, curEigenvectors);

    // calculate the square root of the eigen values
    sqrtEigenValues(eigenvaluesSqr, curEigenvalues);

    expMatrix = (Complex(0.0, mParams.deltaZ()) *
                 prevEigenValues.array()).exp().matrix().asDiagonal();

    if (Bandlimit::Disabled != mParams.bandlimit())
    {
      transferMatrix = curEigenvectors.adjoint() * bandlimitMatrix *
        prevEigenvectors * expMatrix * transferMatrix;
    }
    else
    {
      transferMatrix = curEigenvectors.adjoint() *
        prevEigenvectors * expMatrix * transferMatrix;
    }

    if (verbosity >= 1)
    {
      std::cout << "Current runtime (Slice " << n << "): "
                << duration(high_resolution_clock::now() - startTime)
                << std::endl;
    }
  }

  saveMatrix(transferMatrix, matrixName("F-part", mPart));
}

void Pass2FLarge::run()
{
  auto startTime = high_resolution_clock::now();

  if (mMatrixA2)
  {
    VectorNR eigenvaluesSqr;
    MatrixNC eigenvectors;
    mMatrixA2->calc(mPart.startSlice(), eigenvaluesSqr, eigenvectors);
    saveEVDecomposition(eigenvaluesSqr, eigenvectors,
      tempName("A-", mPart.startSlice()));
  }

  const int64_t numPixels = mParams.numPixels();
  MatrixNC transferMatrix = MatrixNC::Identity(numPixels, numPixels);

  for (size_t n = mPart.startSlice() + 1; n < mPart.endSlice() + 1; ++n)
  {
    {
      VectorNR eigenvaluesSqr;
      if (mMatrixA2)
        loadEigenvalues(tempName("A-", n - 1), eigenvaluesSqr);
      else
        loadEigenvalues(matrixName("A-", n - 1), eigenvaluesSqr);

      // calculate the square root of the eigen values
      VectorNC prevEigenvalues;
      sqrtEigenValues(eigenvaluesSqr, prevEigenvalues);

      DiagonalMatrixNC expMatrix =
        (Complex(0.0, mParams.deltaZ()) *
         prevEigenvalues.array()).exp().matrix().asDiagonal();

      transferMatrix = expMatrix * transferMatrix;
    }

    {
      VectorNR eigenvaluesSqr;
      MatrixNC prevEigenvectors;
      if (mMatrixA2)
      {
        loadEVDecomposition(tempName("A-", n - 1), eigenvaluesSqr,
                            prevEigenvectors);
        remove(tempName("A-", n - 1).c_str());
      }
      else
      {
        loadEVDecomposition(matrixName("A-", n - 1), eigenvaluesSqr,
                            prevEigenvectors);
      }

      transferMatrix = prevEigenvectors * transferMatrix;
    }

    if (Bandlimit::Disabled != mParams.bandlimit())
    {
      MatrixNC bandlimitMatrix;
      loadMatrix(matrixDir() + "/bandlimit.mat", bandlimitMatrix);

      transferMatrix = bandlimitMatrix * transferMatrix;
    }

    {
      VectorNR eigenvaluesSqr;
      MatrixNC curEigenvectors;
      if (mMatrixA2)
      {
        mMatrixA2->calc(n, eigenvaluesSqr, curEigenvectors);
        saveEVDecomposition(eigenvaluesSqr, curEigenvectors, tempName("A-", n));
      }
      else
      {
        loadEVDecomposition(matrixName("A-", n), eigenvaluesSqr,
                            curEigenvectors);
      }

      transferMatrix = curEigenvectors.adjoint() * transferMatrix;
    }

    if (verbosity >= 1)
    {
      std::cout << "Current runtime (Slice " << n << "): "
                << duration(high_resolution_clock::now() - startTime)
                << std::endl;
    }
  }

  saveMatrix(transferMatrix, matrixName("F-part", mPart));
}

void Pass2T::run()
{
  auto startTime = high_resolution_clock::now();

  const int64_t numPix = mParams.numPixels();

  MatrixNC tMatrix;

  VectorNR eigenvaluesSqr;
  VectorNC prevEigenvalues;
  VectorNC curEigenvalues;
  MatrixNC prevEigenvectors;
  MatrixNC curEigenvectors;

  if (mMatrixA2)
  {
    mMatrixA2->calc(mPart.startSlice(), eigenvaluesSqr, curEigenvectors);
  }
  else
  {
    loadEVDecomposition(matrixName("A-", mPart.startSlice()),
                        eigenvaluesSqr, curEigenvectors);
  }

  // calculate the square root of the eigen values
  sqrtEigenValues(eigenvaluesSqr, curEigenvalues);

  for (size_t n = mPart.startSlice() + 1; n < mPart.endSlice() + 1; ++n)
  {
    prevEigenvalues  = curEigenvalues;
    prevEigenvectors = curEigenvectors;

    if (mMatrixA2)
      mMatrixA2->calc(n, eigenvaluesSqr, curEigenvectors);
    else
      loadEVDecomposition(matrixName("A-", n), eigenvaluesSqr, curEigenvectors);

    // calculate the square root of the eigen values
    sqrtEigenValues(eigenvaluesSqr, curEigenvalues);

    MatrixNC uDaggerU = curEigenvectors.adjoint() * prevEigenvectors;

    MatrixNC f = uDaggerU +
                 curEigenvalues.cwiseInverse().asDiagonal() *
                 uDaggerU * prevEigenvalues.asDiagonal();
    f *= 0.5;
    MatrixNC b = uDaggerU -
                 curEigenvalues.cwiseInverse().asDiagonal() *
                 uDaggerU * prevEigenvalues.asDiagonal();
    b *= 0.5;

    if (mParams.debug(Debug::pass2))
    {
      dumpHistogram(f, matrixName("T11-mult", n) + ".txt");
      dumpHistogram(b, matrixName("T21-mult", n) + ".txt");
    }

    // first column of the block matrix
    DiagonalMatrixNC expMatrixPos = (Complex(0.0, mParams.deltaZ()) *
                           prevEigenvalues.array()).exp().matrix().asDiagonal();

    MatrixNC tMultiplier(2 * numPix, 2 * numPix);

    tMultiplier.block(0, 0, numPix, numPix).noalias() = f * expMatrixPos;
    tMultiplier.block(numPix, 0, numPix, numPix).noalias() = b * expMatrixPos;

    // second column of the block matrix
    DiagonalMatrixNC expMatrixNeg = (Complex(0.0, -mParams.deltaZ()) *
                           prevEigenvalues.array()).exp().matrix().asDiagonal();

    tMultiplier.block(0, numPix, numPix, numPix).noalias() = b * expMatrixNeg;
    tMultiplier.block(numPix, numPix, numPix, numPix).noalias() =
      f * expMatrixNeg;

    if (mPart.startSlice() + 1 == n)
      tMatrix = tMultiplier; // first iteration => simply copy
    else
      tMatrix = tMultiplier * tMatrix;

    if (verbosity >= 1)
    {
      std::cout << "Current runtime (Slice " << n << "): "
                << duration(high_resolution_clock::now() - startTime)
                << std::endl;
    }

    if (mParams.debug(Debug::pass2))
    {
      dumpHistogram(tMatrix.block(0, 0, numPix, numPix),
                    matrixName("T11-part", n) + ".txt");
      dumpHistogram(tMatrix.block(0, numPix, numPix, numPix),
                    matrixName("T21-part", n) + ".txt");
    }
  }

  saveMatrix(tMatrix, matrixName("T-part", mPart));
}

void Pass2TLarge::run()
{
  auto startTime = high_resolution_clock::now();

  if (mMatrixA2)
  {
    VectorNR eigenvaluesSqr;
    MatrixNC eigenvectors;
    mMatrixA2->calc(mPart.startSlice(), eigenvaluesSqr, eigenvectors);
    saveEVDecomposition(eigenvaluesSqr, eigenvectors,
      tempName("A-", mPart.startSlice()));
  }

  for (size_t n = mPart.startSlice() + 1; n < mPart.endSlice() + 1; ++n)
  {
    {
      VectorNC prevEigenValues;
      VectorNC curEigenvalues;
      MatrixNC uDaggerU;

      {
        VectorNR eigenvaluesSqr;

        // the eigenvalues and eigenvector for the previous slice
        MatrixNC prevEigenvectors;
        if (mMatrixA2)
        {
          loadEVDecomposition(tempName("A-", n - 1), eigenvaluesSqr,
                              prevEigenvectors);
          remove(tempName("A-", n - 1).c_str());
        }
        else
        {
          loadEVDecomposition(matrixName("A-", n - 1), eigenvaluesSqr,
                              prevEigenvectors);
        }

        // calculate the square root of the eigenvalues
        sqrtEigenValues(eigenvaluesSqr, prevEigenValues);

        MatrixNC curEigenvectors;
        if (mMatrixA2)
        {
          mMatrixA2->calc(n, eigenvaluesSqr, curEigenvectors);
          saveEVDecomposition(eigenvaluesSqr, curEigenvectors,
                              tempName("A-", n));
        }
        else
        {
          loadEVDecomposition(matrixName("A-", n), eigenvaluesSqr,
                              curEigenvectors);
        }

        // calculate the square root of the eigenvalues
        sqrtEigenValues(eigenvaluesSqr, curEigenvalues);

        uDaggerU.noalias() = curEigenvectors.adjoint() * prevEigenvectors;
      }

      // forward-scattering matrix F^n
      MatrixNC f = R(0.5) * uDaggerU +
        Complex(R(0.5)) * curEigenvalues.cwiseInverse().asDiagonal() *
        uDaggerU * prevEigenValues.asDiagonal();

      std::string multiplierTarget;
      if (mPart.startSlice() + 1 == n)
      {
        // first iteration
        if (mPart.endSlice() == n)
        {
          // also last iteration => write directly to the final matrix
          multiplierTarget = matrixDir() + "/T-part" + toString(mPart, 4);
        }
        else
        {
          // write the matrices into temporary buffers
          multiplierTarget = tempDirectory() + "/T-" + toString(n, 6);
        }
      }
      else
      {
        // write the matrices as multipliers into temporary buffers
        multiplierTarget = tempDirectory() + "/TMultiplier-" + toString(n, 6);
      }

      // first column of the block matrix
      DiagonalMatrixNC expMatrixPos = (Complex(0.0, mParams.deltaZ()) *
        prevEigenValues.array()).exp().matrix().asDiagonal();

      {
        MatrixNC tMultiplier11 = f * expMatrixPos;
        saveMatrix(tMultiplier11, multiplierTarget + "-11.mat");
      }

      {
        MatrixNC tMultiplier21 = (uDaggerU - f) * expMatrixPos;
        saveMatrix(tMultiplier21, multiplierTarget + "-21.mat");
      }

      // second column of the block matrix
      DiagonalMatrixNC expMatrixNeg = (Complex(0.0, -mParams.deltaZ()) *
        prevEigenValues.array()).exp().matrix().asDiagonal();

      {
        MatrixNC tMultiplier12 = (uDaggerU - f) * expMatrixNeg;
        saveMatrix(tMultiplier12, multiplierTarget + "-12.mat");
      }

      {
        MatrixNC tMultiplier22 = f * expMatrixNeg;
        saveMatrix(tMultiplier22, multiplierTarget + "-22.mat");
      }
    }

    if (mPart.startSlice() + 1 != n)
    {
      // this is not the first iteration
      // multiply the current result with the previous result
      std::string dstName;

      if (mPart.endSlice() == n)  // last iteration?
        dstName = matrixDir() + "/T-part" + toString(mPart, 4);
      else
        dstName = tempDirectory() + "/T-" + toString(n, 6);

      multiplyLarge(tempDirectory() + "/TMultiplier-" + toString(n, 6),
                    tempDirectory() + "/T-" + toString(n - 1, 6),
                    dstName);

      // we can remove the multipliers
      removeMatricesLarge(tempDirectory() + "/TMultiplier-" + toString(n, 6));

      // and the previous matrices
      removeMatricesLarge(tempDirectory() + "/T-" + toString(n - 1, 6));
    }

    if (verbosity >= 1)
    {
      std::cout << "Current runtime (Slice " << n << "): "
                << duration(high_resolution_clock::now() - startTime)
                << std::endl;
    }
  }
}

void Pass2TSLarge::run()
{
  auto startTime = high_resolution_clock::now();

  {
    // The first matrix A
    VectorNR eigenvaluesSqr;
    MatrixNC eigenvectors;

    if (mMatrixA2)
    {
      mMatrixA2->calc(mPart.startSlice(), eigenvaluesSqr, eigenvectors);
      saveEVDecomposition(eigenvaluesSqr, eigenvectors,
        tempName("A-", mPart.startSlice()));
    }
    else
    {
      loadEVDecomposition(matrixName("A-", mPart.startSlice()),
        eigenvaluesSqr, eigenvectors);
    }
  }

  for (size_t n = mPart.startSlice() + 1; n < mPart.endSlice() + 1; ++n)
  {
    const std::string curUPrevUName  = tempName("curUPrevU-", n);
    const std::string prevUCurUName  = tempName("prevUCurU-", n);

    VectorNC prevEigenvalues;
    VectorNC curEigenvalues;
    {
      VectorNR curEigenvaluesSqr;
      VectorNR prevEigenvaluesSqr;
      MatrixNC curEigenvectors;
      MatrixNC prevEigenvectors;

      // load / calculate the matrix A^2
      if (mMatrixA2)
      {
        mMatrixA2->calc(n, curEigenvaluesSqr, curEigenvectors);
        if (mPart.endSlice() != n)
        {
          saveEVDecomposition(curEigenvaluesSqr, curEigenvectors,
                              tempName("A-", n));
        }

        loadEVDecomposition(tempName("A-", n - 1), prevEigenvaluesSqr,
                            prevEigenvectors);
        remove(tempName("A-", n - 1).c_str());
      }
      else
      {
        loadEVDecomposition(matrixName("A-", n),
                            curEigenvaluesSqr, curEigenvectors);
        loadEVDecomposition(matrixName("A-", n - 1),
                            prevEigenvaluesSqr, prevEigenvectors);
      }

      // calculate the square root of the eigenvalues
      sqrtEigenValues(prevEigenvaluesSqr, prevEigenvalues);
      sqrtEigenValues(curEigenvaluesSqr, curEigenvalues);

      {
        MatrixNC curUPrevU;
        curUPrevU.noalias() = curEigenvectors.adjoint() * prevEigenvectors;
        saveMatrix(curUPrevU, curUPrevUName);
      }
      {
        MatrixNC prevUCurU;
        prevUCurU.noalias() = prevEigenvectors.adjoint() * curEigenvectors;
        saveMatrix(prevUCurU, prevUCurUName);
      }
    }

    // positive exp matrix
    DiagonalMatrixNC expMatrixPos = (Complex(0.0, mParams.deltaZ()) *
                           prevEigenvalues.array()).exp().matrix().asDiagonal();

    // negative exp matrix
    DiagonalMatrixNC expMatrixNeg = (Complex(0.0, -mParams.deltaZ()) *
                           prevEigenvalues.array()).exp().matrix().asDiagonal();

    // generate multiplier target names
    std::string tMultiplierTarget;
    std::string sMultiplierTarget;
    if (mPart.startSlice() + 1 == n) // first iteration?
    {
      if (mPart.endSlice() == n)
      {
        // also last iteration => write directly to the final matrices
        tMultiplierTarget = matrixDir() + "/T-part" + toString(mPart, 4);
        sMultiplierTarget = matrixDir() + "/S-part" + toString(mPart, 4);
      }
      else
      {
        // write the matrices into temporary buffers
        tMultiplierTarget = tempDirectory() + "/T-" + toString(n, 6);
        sMultiplierTarget = tempDirectory() + "/S-" + toString(n, 6);
      }
    }
    else
    {
      // write the matrices as multipliers into temporary buffers
      tMultiplierTarget = tempDirectory() + "/TMultiplier-" + toString(n, 6);
      sMultiplierTarget = tempDirectory() + "/SMultiplier-" + toString(n, 6);
    }

    // T matrix
    {
      MatrixNC f;
      {
        MatrixNC curUPrevU;
        loadMatrix(curUPrevUName, curUPrevU);
        f.noalias() = curEigenvalues.cwiseInverse().asDiagonal() * curUPrevU;
        f *= prevEigenvalues.asDiagonal();
        f += curUPrevU;
      }

      f *= 0.5;

      // diagonal elements
      {
        MatrixNC tMultiplier11;
        tMultiplier11.noalias() = f * expMatrixPos;
        saveMatrix(tMultiplier11, tMultiplierTarget + "-11.mat");

        if (mParams.debug(Debug::pass2))
          dumpHistogram(tMultiplier11, matrixName("T11-mult", n) + ".txt");
      }
      {
        MatrixNC tMultiplier22;
        tMultiplier22.noalias() = f * expMatrixNeg;
        saveMatrix(tMultiplier22, tMultiplierTarget + "-22.mat");
      }

      MatrixNC b;
      {
        MatrixNC curUPrevU;
        loadMatrix(curUPrevUName, curUPrevU);

        b.noalias() = curUPrevU - f;
      }

      // free the memory
      f.resize(0, 0);

      // off-diagonal elements
      {
        MatrixNC tMultiplier21;
        tMultiplier21.noalias() = b * expMatrixPos;
        saveMatrix(tMultiplier21, tMultiplierTarget + "-21.mat");

        if (mParams.debug(Debug::pass2))
          dumpHistogram(tMultiplier21, matrixName("T21-mult", n) + ".txt");
      }
      {
        MatrixNC tMultiplier12;
        tMultiplier12.noalias() = b * expMatrixNeg;
        saveMatrix(tMultiplier12, tMultiplierTarget + "-12.mat");
      }
    }

    if (verbosity >= 1)
    {
      std::cout << "t Multiplier constructed (Slice " << n << "): "
                << duration(high_resolution_clock::now() - startTime)
                << std::endl;
    }

    // S matrix
    {
      MatrixNC fTilde;
      {
        MatrixNC prevUCurU;
        loadMatrix(prevUCurUName, prevUCurU);
        fTilde.noalias() =
          prevEigenvalues.cwiseInverse().asDiagonal() * prevUCurU;
        fTilde *= curEigenvalues.asDiagonal();
        fTilde += prevUCurU;
      }
      fTilde *= 0.5;

      // diagonal elements
      {
        MatrixNC sMultiplier11;
        sMultiplier11.noalias() = expMatrixNeg * fTilde;
        saveMatrix(sMultiplier11, sMultiplierTarget + "-11.mat");

        if (mParams.debug(Debug::pass2))
          dumpHistogram(sMultiplier11, matrixName("S11-mult", n) + ".txt");
      }
      {
        MatrixNC sMultiplier22;
        sMultiplier22.noalias() = expMatrixPos * fTilde;
        saveMatrix(sMultiplier22, sMultiplierTarget + "-22.mat");
      }

      MatrixNC bTilde;
      {
        MatrixNC prevUCurU;
        loadMatrix(prevUCurUName, prevUCurU);

        bTilde.noalias() = prevUCurU - fTilde;
      }

      // free the memory
      fTilde.resize(0, 0);

      // off-diagonal elements
      {
        MatrixNC sMultiplier21;
        sMultiplier21.noalias() = expMatrixPos * bTilde;
        saveMatrix(sMultiplier21, sMultiplierTarget + "-21.mat");

        if (mParams.debug(Debug::pass2))
          dumpHistogram(sMultiplier21, matrixName("S21-mult", n) + ".txt");
      }
      {
        MatrixNC sMultiplier12;
        sMultiplier12.noalias() = expMatrixNeg * bTilde;
        saveMatrix(sMultiplier12, sMultiplierTarget + "-12.mat");
      }
    }

    if (verbosity >= 1)
    {
      std::cout << "s Multiplier constructed (Slice " << n << "): "
                << duration(high_resolution_clock::now() - startTime)
                << std::endl;
    }

    if (mPart.startSlice() + 1 != n)
    {
      // this is not the first iteration
      // multiply the current results with the previous results
      std::string tDstName;
      std::string sDstName;

      if (mPart.endSlice() == n) // last iteration?
      {
        tDstName = matrixDir() + "/T-part" + toString(mPart, 4);
        sDstName = matrixDir() + "/S-part" + toString(mPart, 4);
      }
      else
      {
        tDstName = tempDirectory() + "/T-" + toString(n, 6);
        sDstName = tempDirectory() + "/S-" + toString(n, 6);
      }

      multiplyLarge(tempDirectory() + "/TMultiplier-" + toString(n, 6),
                    tempDirectory() + "/T-" + toString(n - 1, 6),
                    tDstName);

      // we can remove the multipliers and the previous matrices
      removeMatricesLarge(tempDirectory() + "/TMultiplier-" + toString(n, 6));
      removeMatricesLarge(tempDirectory() + "/T-" + toString(n - 1, 6));

      multiplyLarge(tempDirectory() + "/S-"+ toString(n - 1, 6),
                    tempDirectory() + "/SMultiplier-" + toString(n, 6),
                    sDstName);

      removeMatricesLarge(tempDirectory() + "/SMultiplier-" + toString(n, 6));
      removeMatricesLarge(tempDirectory() + "/S-" + toString(n - 1, 6));
    }

    // clean up
    remove(curUPrevUName.c_str());
    remove(prevUCurUName.c_str());

    if (verbosity >= 1)
    {
      std::cout << "Matrices multiplied (Slice " << n << "): "
                << duration(high_resolution_clock::now() - startTime)
                << std::endl;
    }
  }
}

void Pass2TSBandlimitLarge::run()
{
  auto startTime = high_resolution_clock::now();

  {
    // The first matrix A
    VectorNR eigenvaluesSqr;
    MatrixNC eigenvectors;

    if (mMatrixA2)
    {
      mMatrixA2->calc(mPart.startSlice(), eigenvaluesSqr, eigenvectors);
      saveEVDecomposition(eigenvaluesSqr, eigenvectors,
        tempName("A-", mPart.startSlice()));
    }
    else
    {
      loadEVDecomposition(matrixName("A-", mPart.startSlice()),
        eigenvaluesSqr, eigenvectors);
    }

    // calculate curULCurUName for n = mPart.startSlice()
    MatrixNC bandlimit;
    loadMatrix(matrixDir() + "/bandlimit.mat", bandlimit);

    MatrixNC lCurU;
    lCurU.noalias() = bandlimit * eigenvectors;

    // free the bandlimit matrix
    bandlimit.resize(0, 0);

    MatrixNC curULCurU;
    curULCurU.noalias() = eigenvectors.adjoint() * lCurU;

    saveMatrix(curULCurU, tempName("curULCurU-", mPart.startSlice()));
  }

  for (size_t n = mPart.startSlice() + 1; n < mPart.endSlice() + 1; ++n)
  {
    const std::string prevULPrevUName = tempName("curULCurU-" , n - 1);
    const std::string curULCurUName   = tempName("curULCurU-" , n    );
    const std::string curULPrevUName  = tempName("curULPrevU-", n    );
    const std::string prevULCurUName  = tempName("prevULCurU-", n    );

    VectorNC prevEigenvalues;
    VectorNC curEigenvalues;
    VectorNR eigenvaluesSqr;
    {
      // calculate the matrix A
      MatrixNC curEigenvectors;
      if (mMatrixA2)
      {
        mMatrixA2->calc(n, eigenvaluesSqr, curEigenvectors);
        if (mPart.endSlice() != n)
        {
          saveEVDecomposition(eigenvaluesSqr, curEigenvectors,
                              tempName("A-", n));
        }
      }
      else
      {
        loadEVDecomposition(matrixName("A-", n),
                            eigenvaluesSqr, curEigenvectors);
      }

      // calculate the square root of the eigen values
      sqrtEigenValues(eigenvaluesSqr, curEigenvalues);

      MatrixNC bandlimit;
      loadMatrix(matrixDir() + "/bandlimit.mat", bandlimit);

      // U matrices - first step
      {
        MatrixNC curUL = curEigenvectors.adjoint() * bandlimit;
        saveMatrix(curUL, curULPrevUName);
      }
      {
        MatrixNC lCurU;
        lCurU.noalias() = bandlimit * curEigenvectors;
        saveMatrix(lCurU, prevULCurUName);

        // free the bandlimit matrix
        bandlimit.resize(0, 0);

        MatrixNC curULCurU;
        curULCurU.noalias() = curEigenvectors.adjoint() * lCurU;
        saveMatrix(curULCurU, curULCurUName);
      }
    }

    {
      MatrixNC prevEigenvectors;
      if (mMatrixA2)
      {
        loadEVDecomposition(tempName("A-", n - 1), eigenvaluesSqr,
                            prevEigenvectors);
        remove(tempName("A-", n - 1).c_str());
      }
      else
      {
        loadEVDecomposition(matrixName("A-", n - 1),
                            eigenvaluesSqr, prevEigenvectors);
      }

      // calculate the square root of the eigen values
      sqrtEigenValues(eigenvaluesSqr, prevEigenvalues);

      // U matrices - second step
      {
        MatrixNC curUL;
        loadMatrix(curULPrevUName, curUL);
        MatrixNC curULPrevU;
        curULPrevU.noalias() = curUL * prevEigenvectors;
        saveMatrix(curULPrevU, curULPrevUName);
      }
      {
        MatrixNC lCurU;
        loadMatrix(prevULCurUName, lCurU);
        MatrixNC prevULCurU;
        prevULCurU.noalias() = prevEigenvectors.adjoint() * lCurU;
        saveMatrix(prevULCurU, prevULCurUName);
      }
    }

    // positive exp matrix
    DiagonalMatrixNC expMatrixPos = (Complex(0.0, mParams.deltaZ()) *
                           prevEigenvalues.array()).exp().matrix().asDiagonal();

    // negative exp matrix
    DiagonalMatrixNC expMatrixNeg = (Complex(0.0, -mParams.deltaZ()) *
                           prevEigenvalues.array()).exp().matrix().asDiagonal();

    // generate multiplier target names
    std::string tMultiplierTarget;
    std::string sMultiplierTarget;
    if (mPart.startSlice() + 1 == n) // first iteration?
    {
      if (mPart.endSlice() == n)
      {
        // also last iteration => write directly to the final matrices
        tMultiplierTarget = matrixDir() + "/T-part" + toString(mPart, 4);
        sMultiplierTarget = matrixDir() + "/S-part" + toString(mPart, 4);
      }
      else
      {
        // write the matrices into temporary buffers
        tMultiplierTarget = tempDirectory() + "/T-" + toString(n, 6);
        sMultiplierTarget = tempDirectory() + "/S-" + toString(n, 6);
      }
    }
    else
    {
      // write the matrices as multipliers into temporary buffers
      tMultiplierTarget = tempDirectory() + "/TMultiplier-" + toString(n, 6);
      sMultiplierTarget = tempDirectory() + "/SMultiplier-" + toString(n, 6);
    }

    // T matrix
    {
      MatrixNC f;
      {
        MatrixNC curULCurU;
        loadMatrix(curULCurUName, curULCurU);
        curULCurU *= curEigenvalues.cwiseInverse().asDiagonal();

        MatrixNC curULPrevU;
        loadMatrix(curULPrevUName, curULPrevU);
        curULPrevU *= prevEigenvalues.asDiagonal();

        f.noalias() = curULCurU * curULPrevU;
      }

      {
        MatrixNC curULPrevU;
        loadMatrix(curULPrevUName, curULPrevU);
        f += curULPrevU;
      }
      f *= 0.5;

      // diagonal elements
      {
        MatrixNC prevULPrevU;
        loadMatrix(prevULPrevUName, prevULPrevU);
        prevULPrevU *= expMatrixPos;

        MatrixNC tMultiplier11;
        tMultiplier11.noalias() = f * prevULPrevU;
        saveMatrix(tMultiplier11, tMultiplierTarget + "-11.mat");

        if (mParams.debug(Debug::pass2))
          dumpHistogram(tMultiplier11, matrixName("T11-mult", n) + ".txt");
      }
      {
        MatrixNC prevULPrevU;
        loadMatrix(prevULPrevUName, prevULPrevU);
        prevULPrevU *= expMatrixNeg;

        MatrixNC tMultiplier22;
        tMultiplier22.noalias() = f * prevULPrevU;
        saveMatrix(tMultiplier22, tMultiplierTarget + "-22.mat");
      }

      MatrixNC b;
      {
        MatrixNC curULPrevU;
        loadMatrix(curULPrevUName, curULPrevU);

        b.noalias() = curULPrevU - f;
      }

      // free the memory
      f.resize(0, 0);

      // off-diagonal elements
      {
        MatrixNC prevULPrevU;
        loadMatrix(prevULPrevUName, prevULPrevU);
        prevULPrevU *= expMatrixPos;

        MatrixNC tMultiplier21;
        tMultiplier21.noalias() = b * prevULPrevU;
        saveMatrix(tMultiplier21, tMultiplierTarget + "-21.mat");

        if (mParams.debug(Debug::pass2))
          dumpHistogram(tMultiplier21, matrixName("T21-mult", n) + ".txt");
      }
      {
        MatrixNC prevULPrevU;
        loadMatrix(prevULPrevUName, prevULPrevU);
        prevULPrevU *= expMatrixNeg;

        MatrixNC tMultiplier12;
        tMultiplier12.noalias() = b * prevULPrevU;
        saveMatrix(tMultiplier12, tMultiplierTarget + "-12.mat");
      }
    }

    if (verbosity >= 1)
    {
      std::cout << "t Multiplier constructed (Slice " << n << "): "
                << duration(high_resolution_clock::now() - startTime)
                << std::endl;
    }

    // S matrix
    {
      MatrixNC fTilde;
      {
        MatrixNC prevULprevU;
        loadMatrix(prevULPrevUName, prevULprevU);
        prevULprevU *= prevEigenvalues.cwiseInverse().asDiagonal();

        MatrixNC prevULCurU;
        loadMatrix(prevULCurUName, prevULCurU);
        prevULCurU *= curEigenvalues.asDiagonal();

        fTilde.noalias() = prevULprevU * prevULCurU;
      }
      {
        MatrixNC prevULCurU;
        loadMatrix(prevULCurUName, prevULCurU);
        fTilde += prevULCurU;
      }
      fTilde *= 0.5;

      // diagonal elements
      {
        MatrixNC prevULPrevU;
        loadMatrix(prevULPrevUName, prevULPrevU);
        prevULPrevU *= expMatrixNeg;

        MatrixNC sMultiplier11;
        sMultiplier11.noalias() = prevULPrevU  * fTilde;
        saveMatrix(sMultiplier11, sMultiplierTarget + "-11.mat");

        if (mParams.debug(Debug::pass2))
          dumpHistogram(sMultiplier11, matrixName("S11-mult", n) + ".txt");
      }
      {
        MatrixNC prevULPrevU;
        loadMatrix(prevULPrevUName, prevULPrevU);
        prevULPrevU *= expMatrixPos;

        MatrixNC sMultiplier22;
        sMultiplier22.noalias() = prevULPrevU * fTilde;
        saveMatrix(sMultiplier22, sMultiplierTarget + "-22.mat");
      }

      MatrixNC bTilde;
      {
        MatrixNC prevULCurU;
        loadMatrix(prevULCurUName, prevULCurU);

        bTilde.noalias() = prevULCurU - fTilde;
      }

      // free the memory
      fTilde.resize(0, 0);

      // off-diagonal elements
      {
        MatrixNC prevULPrevU;
        loadMatrix(prevULPrevUName, prevULPrevU);
        prevULPrevU *= expMatrixPos;

        MatrixNC sMultiplier21;
        sMultiplier21.noalias() = prevULPrevU * bTilde;
        saveMatrix(sMultiplier21, sMultiplierTarget + "-21.mat");

        if (mParams.debug(Debug::pass2))
          dumpHistogram(sMultiplier21, matrixName("S21-mult", n) + ".txt");
      }
      {
        MatrixNC prevULPrevU;
        loadMatrix(prevULPrevUName, prevULPrevU);
        prevULPrevU *= expMatrixNeg;

        MatrixNC sMultiplier12;
        sMultiplier12.noalias() = prevULPrevU * bTilde;
        saveMatrix(sMultiplier12, sMultiplierTarget + "-12.mat");
      }
    }

    if (verbosity >= 1)
    {
      std::cout << "s Multiplier constructed (Slice " << n << "): "
                << duration(high_resolution_clock::now() - startTime)
                << std::endl;
    }

        // clean up
    remove(curULPrevUName.c_str());
    remove(prevULCurUName.c_str());
    remove(prevULPrevUName.c_str());

    if (mPart.endSlice() == n) // last slice
      remove(curULCurUName.c_str());

    if (mPart.startSlice() + 1 != n)
    {
      // this is not the first iteration
      // multiply the current results with the previous results
      std::string tDstName;
      std::string sDstName;

      if (mPart.endSlice() == n) // last iteration?
      {
        tDstName = matrixDir() + "/T-part" + toString(mPart, 4);
        sDstName = matrixDir() + "/S-part" + toString(mPart, 4);
      }
      else
      {
        tDstName = tempDirectory() + "/T-" + toString(n, 6);
        sDstName = tempDirectory() + "/S-" + toString(n, 6);
      }

      multiplyLarge(tempDirectory() + "/TMultiplier-" + toString(n, 6),
                    tempDirectory() + "/T-" + toString(n - 1, 6),
                    tDstName);

      // we can remove the multipliers and the previous matrices
      removeMatricesLarge(tempDirectory() + "/TMultiplier-" + toString(n, 6));
      removeMatricesLarge(tempDirectory() + "/T-" + toString(n - 1, 6));

      multiplyLarge(tempDirectory() + "/S-"+ toString(n - 1, 6),
                    tempDirectory() + "/SMultiplier-" + toString(n, 6),
                    sDstName);

      removeMatricesLarge(tempDirectory() + "/SMultiplier-" + toString(n, 6));
      removeMatricesLarge(tempDirectory() + "/S-" + toString(n - 1, 6));
    }

    if (verbosity >= 1)
    {
      std::cout << "Matrices multiplied (Slice " << n << "): "
                << duration(high_resolution_clock::now() - startTime)
                << std::endl;
    }
  }
}

} // namespace Aurora
