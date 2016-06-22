//--- Aurora/AuroraLib/Crytal.cpp ----------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Crystal.hpp"

#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/Sample.hpp"

#include <iostream>
#include <set>

namespace Aurora
{

struct LatticePoint
{
  LatticePoint(int64_t _i, int64_t _j, int64_t _k) :
    i(_i), j(_j), k(_k)
  {}

  int64_t i, j, k;

  bool operator<(const LatticePoint& other) const
  {
    return (i <  other.i) ||
           (i == other.i && j <  other.j) ||
           (i == other.i && j == other.j && k < other.k);
  }
};

SamplePtr Crystal::generateSample(const Vector3R& minBox,
                                  const Vector3R& maxBox,
                                  const Matrix3x3R& rotation,
                                  const Vector3R& translation,
                                  size_t maxAtoms) const
{
  // For every given lattice point we look at its neighbour to check whether or
  // not they are inside the box. If they are we will check these points in
  // turn. For the management of the points we use two different: the already
  // checked points and the unchecked points.

  const size_t numAtomsBase = mBase.size();

  // convert the coordinates of the base to the final coordinate system
  std::vector<Vector3R> atomCoordsBase;

  for (auto atom : mBase)
  {
    // coordinates of the atom iwith respect to the unit cell
    const Vector3R baseVector = atom.position();
    // coordinates in an unrotated cartesian coordinate system
    const Vector3R atomVectorCrystal = baseVector.x() * mLatticeVec1 +
                                       baseVector.y() * mLatticeVec2 +
                                       baseVector.z() * mLatticeVec3;
    // coordinates in the fina cartesian coordinate system
    atomCoordsBase.push_back(rotation * atomVectorCrystal);
  }

  auto result = std::make_shared<Sample>();
  AtomVector& atoms = result->atoms();

  std::set<LatticePoint> checkedPoints;
  std::set<LatticePoint> uncheckedPoints;
  uncheckedPoints.insert(LatticePoint(0, 0, 0));

  while (!uncheckedPoints.empty())
  {
    auto currentPointIt = uncheckedPoints.begin();

    // coordinates of the lattice point in an unrotated cartesian coordinate
    // system
    const Vector3R latticePointCrystal =
      currentPointIt->i * mLatticeVec1 +
      currentPointIt->j * mLatticeVec2 +
      currentPointIt->k * mLatticeVec3;
    // coordinates of the lattice point in the final cartesian coordinate system
    const Vector3R latticePoint = rotation * latticePointCrystal + translation;

    // check, if the lattice point is inside our box
    if (minBox.x() <= latticePoint.x() && latticePoint.x() < maxBox.x() &&
        minBox.y() <= latticePoint.y() && latticePoint.y() < maxBox.y() &&
        minBox.z() <= latticePoint.z() && latticePoint.z() < maxBox.z())
    {
      // add atoms of the base to the sample
      for (size_t i = 0; i < numAtomsBase; ++i)
      {
        atoms.push_back(Atom(latticePoint + atomCoordsBase[i],
                        mBase[i].atomicNumber()));
      }

      // add all 26 neighbours
      for (int64_t i = -1; i <= 1; ++i)
      {
        for (int64_t j = -1; j <= 1; ++j)
        {
          for (int64_t k = -1; k <= 1; ++k)
          {
            // skip the center
            if ((0 == i) && (0 == j) && (0 == k))
              continue;

            LatticePoint point(currentPointIt->i + i,
                               currentPointIt->j + j,
                               currentPointIt->k + k);
            if (checkedPoints.find(point) == checkedPoints.end())
              uncheckedPoints.insert(point);
          }
        }
      }

      if (atoms.size() > maxAtoms)
      {
        std::cerr << Aurora::Color::red << Aurora::Terminal::bold
                  << "Maximum number of atoms reached.\n"
                  << Aurora::Terminal::reset;
        return result;
      }
    }

    checkedPoints.insert(*currentPointIt);
    uncheckedPoints.erase(currentPointIt);
  }

  return result;
}

Graphene::Graphene()
{
#if 0
  // C-C distance
  Real aCC = 1.42;
  // The lattice vectors for graphene. It is a triangular lattice.
  setLatticeVec1(Vector3(std::sqrt(3.0) * aCC, 0.0, 0.0));
  setLatticeVec2(Vector3(0.5 * std::sqrt(3.0) * aCC, 1.5 * aCC, 0.0));
  // half the unit cell height of A-B stacked graphite
  setLatticeVec3(Vector3(0.5 * std::sqrt(3.0) * aCC, 0.5 * aCC, 3.35));
#endif

  setLatticeVec1(Vector3R(R(2.46), R(0.0 ), R(0.0 )));
  setLatticeVec2(Vector3R(R(1.23), R(2.13), R(0.0 )));
  setLatticeVec3(Vector3R(R(1.23), R(0.71), R(3.35)));

  // Graphene has a base with two elements
  base().push_back(Atom(Vector3R(0, 0, 0), 6));
  base().push_back(Atom(Vector3R(R(1.0) / R(3.0), R(1.0) / R(3.0), 0), 6));
}

Silicon::Silicon()
{
  // supercube size
  Real a = R(5.43);
  setLatticeVec1(Vector3R(a, 0, 0));
  setLatticeVec2(Vector3R(0, a, 0));
  setLatticeVec3(Vector3R(0, 0, a));

  base().push_back(Atom(Vector3R(R(0.0 ), R(0.0 ), R(0.0 )), 14));
  base().push_back(Atom(Vector3R(R(0.5 ), R(0.5 ), R(0.0 )), 14));
  base().push_back(Atom(Vector3R(R(0.0 ), R(0.5 ), R(0.5 )), 14));
  base().push_back(Atom(Vector3R(R(0.5 ), R(0.0 ), R(0.5 )), 14));
  base().push_back(Atom(Vector3R(R(0.25), R(0.25), R(0.25)), 14));
  base().push_back(Atom(Vector3R(R(0.75), R(0.75), R(0.25)), 14));
  base().push_back(Atom(Vector3R(R(0.75), R(0.25), R(0.75)), 14));
  base().push_back(Atom(Vector3R(R(0.25), R(0.75), R(0.75)), 14));
}

SmBaCuO::SmBaCuO()
{
  setLatticeVec1(Vector3R(R(3.878), R(0.0  ), R( 0.0  )));
  setLatticeVec2(Vector3R(R(  0.0), R(3.878), R( 0.0  )));
  setLatticeVec3(Vector3R(R(  0.0), R(0.0  ), R(11.818)));

  Real z1 = R(0.19249);
  Real z2 = R(0.35607);
  Real z3 = R(0.1531 );
  Real z4 = R(0.3735 );

  // Sm
  base().push_back(Atom(Vector3R(R(0.5), R(0.5), R(0.5)), 62));
  // Ba
  base().push_back(Atom(Vector3R(R(0.5), R(0.5), z1), 56));
  base().push_back(Atom(Vector3R(R(0.5), R(0.5), R(1.0) - z1), 56));
  // Cu
  base().push_back(Atom(Vector3R(R(0.0), R(0.0), R(0.0)), 29));
  base().push_back(Atom(Vector3R(R(0.0), R(0.0), z2), 29));
  base().push_back(Atom(Vector3R(R(0.0), R(0.0), R(1.0) - z2), 29));
  // O
  base().push_back(Atom(Vector3R(R(0.0), R(0.0), z3), 8));
  base().push_back(Atom(Vector3R(R(0.0), R(0.0), R(1.0) - z3), 8));
  base().push_back(Atom(Vector3R(R(0.0), R(0.5), R(0.0)), 8));
  base().push_back(Atom(Vector3R(R(0.5), R(0.0), R(0.0)), 8));
  base().push_back(Atom(Vector3R(R(0.0), R(0.5), z4), 8));
  base().push_back(Atom(Vector3R(R(0.5), R(0.0), z4), 8));
  base().push_back(Atom(Vector3R(R(0.0), R(0.5), R(1.0) - z4), 8));
  base().push_back(Atom(Vector3R(R(0.5), R(0.0), R(1.0) - z4), 8));
}

Gold::Gold()
{
  Real a = R(4.065);
  setLatticeVec1(Vector3R(a         , R(0.0)    , R(0.0)));
  setLatticeVec2(Vector3R(R(0.5) * a, R(0.5) * a, R(0.0)));
  setLatticeVec3(Vector3R(R(0.0)    , R(0.5) * a, R(0.5) * a));

  base().push_back(Atom(Vector3R(0, 0, 0), 79));
}

Crystal::Creators& Crystal::creators()
{
  static Creators creators;
  return creators;
}

class RegisterCrystals
{
public:
  RegisterCrystals()
  {
    Crystal::registerClass<Graphene>("Graphene");
    Crystal::registerClass<Silicon>("Silicon");
    Crystal::registerClass<SmBaCuO>("SmBaCuO");
    Crystal::registerClass<Gold>("Gold");
  }
} registerCrystals;

} // namespace Aurora
