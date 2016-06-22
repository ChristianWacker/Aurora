//--- Aurora/AuroraLib/Prerequisites.hpp ---------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------
//
/// @file
/// @brief This header should be the first header to be included in every other
/// header. It contains common definitions and configuration options.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_PREREQUISITES_HPP
#define AURORA_AURORA_LIB_PREREQUISITES_HPP

#include "AuroraLib/Config.hpp"

// compiler-dependent stuff
#include "AuroraLib/CompilerAbstraction.hpp"
// the assert macro
#include <cassert>

#ifndef NDEBUG
# define AURORA_ASSUME(expression) assert((expression) && "Assumption failed.")
#else
# define AURORA_ASSUME(expression) AURORA_ASSUME_INTERN(expression)
#endif

/// Macro to mark code paths that are logically unreachable
#define AURORA_UNREACHABLE AURORA_ASSUME(false)

namespace Aurora
{
  typedef unsigned char uchar;
}

#if (1 == AURORA_USE_LAPACK)
// Use the 64-bit interface for the Intel MKL
# define MKL_ILP64
// Eigen should use Intel MKL as BLAS backend
# define EIGEN_USE_BLAS
#endif

#include <cstdint>
#include <memory>

// Declaration of types and classes to reduce the header dependencies.
namespace Aurora
{

#if (1 == AURORA_SINGLE_PRECISION)
  typedef float Real;
# define R(x) x ## f
#else
  typedef double Real;
# define R(x) x
#endif

class Atom;
enum class Bandlimit;
class Buffer2DInfo;
class Buffer2DInfoBase;
class CBandlimit;
class CBuffer2DCodec;

template<int64_t margin>
class CMarginBuffer2D;
typedef CMarginBuffer2D<0> CBuffer2D;

class Crystal;
class CtfCoherent;
class ElasticParameters;
class Element;
class Energy;
class Fft;
class FourierRingCorrelation;
class Histogram;
class Image;
class Memory;
class MicroscopeSimulation;
class FormFactorParametrization;
class PhaseShift;
class Potential;
class RBandlimit;
class RBuffer2D;
class StemDetector;
class StemSimulation;
class TemSimulation;
class Propagator;
class MultisliceParameters;

class Sample;
typedef std::shared_ptr<Sample> SamplePtr;
typedef std::shared_ptr<const Sample> ConstSamplePtr;

class Stem;

} // namespace Aurora

#endif
