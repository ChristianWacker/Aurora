//--- Aurora/AuroraLib/Sample.hpp --------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_SAMPLE_HPP
#define AURORA_AURORA_LIB_SAMPLE_HPP

#include "AuroraLib/Atom.hpp"
#include "AuroraLib/Loader.hpp"
#include "AuroraLib/Saver.hpp"

#include <map>
#include <vector>

namespace Aurora
{

typedef std::vector<Atom> AtomVector;

class AURORA_API Sample : public Loadable<Sample>, public Savable<Sample>
{
public:
  Sample() = default;

  Sample(const std::string& filename);

  /// Calculate a bound box where the faces parallel to the x-y-plane are
  /// quadratic
  Vector3R minQuadraticBoundingBox(Real scale = 1) const
  {
    // the center of the bounding box
    Vector3R center = R(0.5) * (mMinBoundingBox + mMaxBoundingBox);

    // the diagonal of the bounding box
    Vector3R diagonal = mMaxBoundingBox - mMinBoundingBox;

    // the greatest extension of the bounding box in the x-y-plane, this will
    // be the edge length of the square.
    Real maxExtension = std::max(diagonal.x(), diagonal.y());
    return Vector3R(center.x() - maxExtension * R(0.5) * scale,
                    center.y() - maxExtension * R(0.5) * scale,
                    mMinBoundingBox.z());
  }

  /// Calculate a bound box where the faces parallel to the x-y-plane are
  /// quadratic
  Vector3R maxQuadraticBoundingBox(Real scale = 1) const
  {
    // the center of the bounding box
    Vector3R center = R(0.5) * (mMinBoundingBox + mMaxBoundingBox);

    // the diagonal of the bounding box
    Vector3R diagonal = mMaxBoundingBox - mMinBoundingBox;

    // the greatest extension of the bounding box in the x-y-plane, this will
    // be the edge length of the square.
    Real maxExtension = std::max(diagonal.x(), diagonal.y());
    return Vector3R(center.x() + maxExtension * R(0.5) * scale,
                    center.y() + maxExtension * R(0.5) * scale,
                    mMaxBoundingBox.z());
  }

  Vector3R minBoundingBox() const
  {
    return mMinBoundingBox;
  }

  Vector3R maxBoundingBox() const
  {
    return mMaxBoundingBox;
  }

  size_t numAtoms() const
  {
    return mAtoms.size();
  }

  void preprocess();

  RBuffer2D projectedPotential(int64_t cols, int64_t rows, Vector3R minCoords,
    Vector3R maxCoords, Real cutoff,
    const std::string& parametrization = "Peng") const;

  void load(const std::string& filename);

  AtomVector& atoms()
  {
    return mAtoms;
  }

  const AtomVector& atoms() const
  {
    return mAtoms;
  }

  static Loadable<Sample>::Loaders& loaders();
  static Savable<Sample>::Savers& savers();

private:
  AtomVector mAtoms;
  Vector3R mMinBoundingBox;
  Vector3R mMaxBoundingBox;
};

} // namespace Aurora

#endif
