//--- Aurora/AuroraLib/Sample.cpp --------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Sample.hpp"

#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Element.hpp"
#include "AuroraLib/Exception.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/Math.hpp"
#include "AuroraLib/RBuffer2D.hpp"

#include <algorithm>
#include <chrono>
#include <iostream>

namespace Aurora
{

Sample::Sample(const std::string& filename)
{
  load(filename);
}

void Sample::preprocess()
{
  Real minReal = std::numeric_limits<Real>::lowest();
  Real maxReal = std::numeric_limits<Real>::max();

  // determine the greatest and smallest value for x, y and z. That is an
  // axis-aligned bounding box for all atoms
  mMaxBoundingBox = Vector3R(minReal, minReal, minReal);
  mMinBoundingBox = Vector3R(maxReal, maxReal, maxReal);

  for (auto i : mAtoms)
  {
    mMaxBoundingBox = max(mMaxBoundingBox, i.position());
    mMinBoundingBox = min(mMinBoundingBox, i.position());
  }

  // lambda expression
  auto compareZCoords = [] (const Atom& a, const Atom& b)
  {
    return a.position().z() < b.position().z();
  };

  // sort the atoms for their z-Coordinates
  std::sort(mAtoms.begin(), mAtoms.end(), compareZCoords);

  if (verbosity >= 3)
  {
    for (auto i : mAtoms)
      std::cout << i.position() << '\n';
  }
}

/*static*/ void Sample::load(const std::string& filename)
{
  if (verbosity >= 1)
    std::cout << "Load sample...\n";

  using namespace std::chrono;
  auto startTime = high_resolution_clock::now();

  // call the load function of the base class
  Loadable<Sample>::load(filename);
  auto stopTime = high_resolution_clock::now();

  if (verbosity >= 1)
  {
    std::cout << "Loaded " << mAtoms.size() << " atoms in "
              << duration(stopTime - startTime) << ":\n";
  }

  if (verbosity == 0)
    return;

  std::vector<uint32_t> counts(Element::mNumElements, 0);

  for (auto atom : mAtoms)
    ++counts[atom.atomicNumber() - 1];

  // for all atomic numbers
  for (size_t i = 1, numElements = Element::mNumElements; i <= numElements; ++i)
  {
    if (0 != counts[i - 1])
    {
      std::cout << ' ' << Element::symbol(i) << ": "
        << counts[i - 1] << " atoms\n";
    }
  }
}

RBuffer2D Sample::projectedPotential(int64_t cols, int64_t rows,
                                     Vector3R minCoords, Vector3R maxCoords,
                                     Real cutoff,
                                     const std::string& parametrization) const
{
  RBuffer2D projPotential(cols, rows, R(0.0));

  auto periodicTable = FormFactorParametrization::create(parametrization);

  const Real cutoffSqr = cutoff * cutoff;

  for (auto atom : mAtoms)
  {
    Vector3R atomPosition = atom.position();
    const Element& element = periodicTable->element(atom.atomicNumber(),
                                                    atom.charge());

    #pragma omp parallel for
    for (int64_t j = 0; j < rows; ++j)
    {
      const Real fracY = static_cast<Real>(j) / static_cast<Real>(rows);
      const Real y = lerp(minCoords.y(), maxCoords.y(), fracY);
      const Real deltaYSqr = powerOf<2>(y - atomPosition.y());

      if (deltaYSqr > cutoffSqr)
        continue;

      for (int64_t i = 0; i < cols; ++i)
      {
        const Real fracX = static_cast<Real>(i) / static_cast<Real>(cols);
        const Real x = lerp(minCoords.x(), maxCoords.x(), fracX);
        const Real deltaXSqr = powerOf<2>(x - atomPosition.x());
        const Real rhoSqr = deltaXSqr + deltaYSqr;

        if (rhoSqr > cutoffSqr)
          continue;

        projPotential.pixel(i, j) += element.projectedPotential(rhoSqr);
      }
    }
  }

  return projPotential;
}

Sample::Loaders& Sample::loaders()
{
  static Loaders loaders;
  return loaders;
}

Sample::Savers& Sample::savers()
{
  static Savers savers;
  return savers;
}

} // namespace Aurora
