//--- Aurora/Clients/Matrix/Pass3.cpp ------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "Pass3.hpp"

#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Formatter.hpp"

#include <chrono>
#include <iostream>

namespace Aurora
{

using namespace std::chrono;

void Pass3::run()
{
  auto startTime = high_resolution_clock::now();

  MatrixNC transferMatrix;
  loadMatrix(matrixName(mMatrixPrefix + "-part", 1), transferMatrix);
  for (size_t i = 2; i <= mNumParts; ++i)
  {
    MatrixNC temp;
    loadMatrix(matrixName(mMatrixPrefix + "-part", i), temp);

    // the order of the factors depend on the propagation direction
    if (Direction::forward == mDirection)
      transferMatrix = temp * transferMatrix;
    else
      transferMatrix = transferMatrix * temp;

    if (verbosity >= 1)
    {
      std::cout << "Elapsed time (Part " << i << "): "
                << duration(high_resolution_clock::now() - startTime)
                << std::endl;
    }
  }

  // transferMatrix contains now the transfer matrix for one symmetry unit
  saveMatrix(transferMatrix, matrixName(mMatrixPrefix + "-power", 1));
}

void Pass3Large::run()
{
  auto startTime = high_resolution_clock::now();

  for (size_t part = 2; part <= mNumParts; ++part)
  {
    std::string srcName1;
    if (2 == part) // first multiplication?
    {
      srcName1 = matrixDir() + "/" + mMatrixPrefix + "-part" + toString(1, 4);
    }
    else
    {
      srcName1 =
        tempDirectory() + "/" + mMatrixPrefix + "-cum" + toString(part - 1, 4);
    }

    std::string dstName;
    if (mNumParts == part) // last multiplication?
    {
      // write directly to the target
      dstName = matrixDir() + "/" + mMatrixPrefix + "-power" + toString(1, 4);
    }
    else
    {
      // cumulative result in a temporaty directory
      dstName =
        tempDirectory() + "/" + mMatrixPrefix + "-cum" + toString(part, 4);
    }

    std::string srcName2 =
      matrixDir() + "/" + mMatrixPrefix + "-part" + toString(part, 4);

    // the order of the factors depend on the propagation direction
    if (Direction::forward == mDirection)
      multiplyLarge(srcName2, srcName1, dstName);
    else
      multiplyLarge(srcName1, srcName2, dstName);

    // remove the temporary files we have created
    if (part > 2)
    {
      removeMatricesLarge(tempDirectory() + "/" + mMatrixPrefix + "-cum" +
                          toString(part - 1, 4));
    }

    if (verbosity >= 1)
    {
      std::cout << "Elapsed time (Part " << part << "): "
                << duration(high_resolution_clock::now() - startTime)
                << std::endl;
    }
  }
}

} // namespace Aurora
