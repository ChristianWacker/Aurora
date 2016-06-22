//--- Aurora/AuroraLib/Crytal.hpp ----------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_CRYSTAL_HPP
#define AURORA_AURORA_LIB_CRYSTAL_HPP

#include "AuroraLib/Atom.hpp"
#include "AuroraLib/Factory.hpp"
#include "AuroraLib/Math.hpp"

#include <Eigen/Geometry>
#include <vector>

namespace Aurora
{

// todo: Documentation
// todo: Test Case
class AURORA_API Crystal : public Factory<Crystal>
{
public:
  Real volume() const
  {
    return std::fabs(mLatticeVec1.dot(mLatticeVec2.cross(mLatticeVec3)));
  }

  Vector3R latticeVec1() const
  {
    return mLatticeVec1;
  }

  Vector3R latticeVec2() const
  {
    return mLatticeVec2;
  }

  Vector3R latticeVec3() const
  {
    return mLatticeVec3;
  }

  Vector3R reciprocalLatticeVec1() const
  {
    return Math::twoPi * mLatticeVec2.cross(mLatticeVec3) /
           mLatticeVec1.dot(mLatticeVec2.cross(mLatticeVec3));
  }

  Vector3R reciprocalLatticeVec2() const
  {
    return Math::twoPi * mLatticeVec3.cross(mLatticeVec1) /
           mLatticeVec2.dot(mLatticeVec3.cross(mLatticeVec1));
  }

  Vector3R reciprocalLatticeVec3() const
  {
    return Math::twoPi * mLatticeVec1.cross( mLatticeVec2) /
           mLatticeVec3.dot(mLatticeVec1.cross(mLatticeVec2));
  }

  SamplePtr generateSample(const Vector3R& minBox,
                           const Vector3R& maxBox,
                           const Matrix3x3R& rotation,
                           const Vector3R& translation,
                           size_t maxAtoms = 100000000) const;

  void setLatticeVec1(const Vector3R& vec)
  {
    mLatticeVec1 = vec;
  }

  void setLatticeVec2(const Vector3R& vec)
  {
    mLatticeVec2 = vec;
  }

  void setLatticeVec3(const Vector3R& vec)
  {
    mLatticeVec3 = vec;
  }

  typedef std::vector<Atom> BaseVector;
  BaseVector& base()
  {
    return mBase;
  }

  static Creators& creators();

private:
  Vector3R mLatticeVec1;
  Vector3R mLatticeVec2;
  Vector3R mLatticeVec3;

  BaseVector mBase;
};

/// Creates a graphene crystal. The C-C distance is
/// @f$a_{CC} = 1.42\,\text{A}@f$. The distance of two graphene planes is half
/// the unit cell height of A-B stacked graphite @f$3.35\,\text{A}@f$
/// @details Periodic Box:
/// @f$\sqrt{3}\,a_{CC}\times3\,a_{CC}=2.46\,\text{A}\times4.26\,\text{A}@f$
class Graphene : public Crystal
{
public:
  Graphene();
};

/// Creates a silicon crystal. The length of the edge of the cubic super cell
/// is @f$a=5.43\,\text{A}@f$.
/// @details Periodic Box [110]:
/// @f$a\times a/\sqrt{2}=5.430\,\text{A}\times3.840\,\text{A}@f$
class Silicon : public Crystal
{
public:
  Silicon();
};

/// Creates the crystal @f$ \text{SmBa}_2\text{Cu}_3\text{O}_{7-x} @f$
/// [Ming2013]. The size of the super cell is
/// @f$3.878\,\text{A} \times 3.878\,\text{A} \times 11.818\,\text{A}@f$.
class SmBaCuO : public Crystal
{
public:
  SmBaCuO();
};

class Gold : public Crystal
{
public:
  Gold();
};

} // namespace Aurora

#endif
