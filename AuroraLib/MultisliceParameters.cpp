//--- Aurora/AuroraLib/MultisliceParameters.cpp --------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/MultisliceParameters.hpp"

#include "AuroraLib/CommandLine.hpp"

#include <iostream>

#ifdef _OPENMP
# include <omp.h>
#endif

namespace Aurora
{

void MultisliceParameters::setBox(const Vector3R& minBox,
                                  const Vector3R& maxBox,
                                  Real deltaZ)
{
  mMinBox = minBox;
  mMaxBox = maxBox;

  setSize(maxBox.x() - minBox.x(), maxBox.y() - minBox.y());
  mDepth = maxBox.z() - minBox.z();
  mDeltaZ = deltaZ;

  if (mDeltaZ <= 0)
    AURORA_THROW(EInvalidParameter, "deltaZ must be positive");

  mMinBoxes.clear();
  Real z = minBox.z();
  const Real maxZ = maxBox.z();

  mMinBoxes.push_back(Vector3R(minBox.x(), minBox.y(), z));

  do
  {
    z += mDeltaZ;
    mMinBoxes.push_back(Vector3R(minBox.x(), minBox.y(), z));
  }
  while (z < maxZ);

  if (verbosity >= 2)
  {
    for (size_t i = 0; i < numSlices(); ++i)
    {
      std::cout << " slice " << i << " (z1 = " << zStart(i) << ", z2 = "
                << zEnd(i) << ")\n";
    }
  }

  if (verbosity >= 1)
    std::cout << "Number of Slices: " << numSlices() << std::endl;
}

void MultisliceParameters::setJitNumSlices(int64_t numSlices)
{
  assert(numSlices >= -1 && "numSlices out of range.");

  if (-1 == numSlices)
  {
  #ifdef _OPENMP
    mJitNumSlices = omp_get_max_threads();
  #else
    mJitNumSlices = 1;
  #endif
  }
  else
  {
    mJitNumSlices = numSlices;
  }
}

} // namespace Aurora
