//--- Aurora/AuroraLib/Atom.hpp ------------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_ATOM_HPP
#define AURORA_AURORA_LIB_ATOM_HPP

#include "AuroraLib/Math.hpp"

namespace Aurora
{

/// This class describes one atom of a @ref Sample.
class AURORA_API Atom
{
public:
  /// Creates a new Atom.
  /// @param position
  ///  The position of the new Atom in the sample described with three
  ///  cartesian coordinates.
  /// @param atomicNumber
  ///  Atomic number of the new atom.
  /// @param charge
  ///  The charge of the atom in units of the elementary charge.
  Atom(const Vector3R& position, size_t atomicNumber, int charge = 0) :
    mPosition(position), mAtomicNumber(atomicNumber), mCharge(charge),
    mBFactor(0)
  { }

  size_t atomicNumber() const
  {
    return mAtomicNumber;
  }

  Real bFactor() const
  {
    return mBFactor;
  }

  int charge() const
  {
    return mCharge;
  }

  /// Returns the position of the atom.
  Vector3R position() const
  {
    return mPosition;
  }

  void setAtomicNumber(size_t atomicNumber)
  {
    mAtomicNumber = atomicNumber;
  }

  void setBFactor(Real bFactor)
  {
    mBFactor = bFactor;
  }

  void setCharge(int charge)
  {
    mCharge = charge;
  }

  void setPosition(const Vector3R& position)
  {
    mPosition = position;
  }

private:
  Vector3R mPosition;
  size_t   mAtomicNumber;
  int      mCharge;
  Real     mBFactor;
};

} // namespace Aurora

#endif
