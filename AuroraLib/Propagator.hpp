//--- Aurora/AuroraLib/Propagator.hpp ------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_PROPAGATOR_HPP
#define AURORA_AURORA_LIB_PROPAGATOR_HPP

#include "AuroraLib/Factory.hpp"
#include "AuroraLib/Laplace.hpp"
#include "AuroraLib/MultisliceParameters.hpp"
#include "AuroraLib/PhaseShift.hpp"

#include <Eigen/Sparse>
#include <vector>

namespace Aurora
{

/// Creates a Fresnel propagator for a slice of thickness deltaZ
AURORA_API CBuffer2D fresnelPropagator(const Buffer2DInfo& info, Real deltaZ,
                                       Bandlimit bandlimit);

/// Base class for the description of the wave propagation through a sample.
/// This class uses the factory pattern. Hence, concrete implementations can be
/// created with a call to @ref create().
class AURORA_API Propagator :
  public Factory<Propagator, CBuffer2D&, const MultisliceParameters&,
                 const ConstSamplePtr&>
{
public:
  virtual ~Propagator() {}
  /// Starts the propagation
  virtual void propagate() = 0;

  typedef std::function<void (size_t, size_t, const CBuffer2D&)> OnWave;
  OnWave onWave;

  static Creators& creators();

protected:
  /// Constructor
  /// @param wave
  ///  The wave that should be propagated.
  /// @param params
  ///  The simulation parameters
  Propagator(CBuffer2D& wave, const MultisliceParameters& params);

  void reportDivergence(const std::string& additionalText = "");

  const MultisliceParameters& mParams;
  // The wave we are propagating.
  CBuffer2D& mWave;

  Real mZ;

  int64_t mNumDivergences = 0;
};

class PropagatorPhaseShift : public Propagator
{
public:
  PropagatorPhaseShift(CBuffer2D& wave, const MultisliceParameters& params,
                       const ConstSamplePtr& sample, PhaseShiftMode mode);

protected:
  std::unique_ptr<PhaseShift> mPhaseShift;
};

class PropagatorPotential : public Propagator
{
public:
  PropagatorPotential(CBuffer2D& wave, const MultisliceParameters& params,
                      const ConstSamplePtr& sample);

protected:
  std::unique_ptr<Potential> mPotential;
};

/// This class implements the propagation in Fourier space.
class PropagatorClassicalFT : public PropagatorPhaseShift
{
public:
  /// Constructor
  /// @param wave
  ///  The wavefunction that should be propagated.
  /// @param params
  ///  The parameters describing the TEM simulation.
  PropagatorClassicalFT(CBuffer2D& wave, const MultisliceParameters& params,
                        const ConstSamplePtr& sample);
  void propagate() final;

private:
  // we will use two buffers for the wave: one in direct space and one in
  // Fourier space
  CBuffer2D mWaveFft;
  // The necessary FFTs to convert between the representations in real and
  // Fourier space
  const std::unique_ptr<Fft> mFft;
  const std::unique_ptr<Fft> mFftInv;

  // Fresnel propagator for a half slice (deltaZ / 2). This used for the first
  // and the last slice.
  const CBuffer2D mPropagatorHalf;
  // Fresnel propagator for a full slice (deltaZ).
  const CBuffer2D mPropagatorFull;
};

/// Propagator in real space using a series expansion of the exponential
/// function of the Laplacian and a separate calculation of the transfer
/// function.
template<class Method>
class PropagatorClassical : public PropagatorPhaseShift
{
public:
  PropagatorClassical(CBuffer2D& wave, const MultisliceParameters& params,
                      const ConstSamplePtr& sample);
  void propagate() final;

private:
  /// Propagates the wave through one slice.
  void propagateOneSlice(CBuffer2D& wave, Real deltaZ);

  CBuffer2D mTemp1;
  // we need a temporary buffer with a margin to calculate the discrete
  // Laplacian
  CMarginBuffer2D<Method::margin> mTemp2;

  std::unique_ptr<CBandlimit> mBandlimit;
};

/// Propagator in real space using an series expansion of the exponential
/// function and the discrete Laplacian
template<class Method>
class PropagatorSingle : public PropagatorPhaseShift
{
public:
  PropagatorSingle(CBuffer2D& wave, const MultisliceParameters& params,
                   const ConstSamplePtr& sample);
  void propagate() final;

private:
  /// Propagates the wave through one slice. This steps combines the
  /// propagation and the application of the phase shift.
  void propagateOneSlice(CBuffer2D& waveIn, CBuffer2D& waveOut,
                         const RBuffer2D& phaseShift);

  CBuffer2D mWaveAfter;
  // we need a temporary buffer with on pixel margin to calculate the discrete
  // Laplacian
  CMarginBuffer2D<Method::margin> mTemp;

  std::unique_ptr<CBandlimit> mBandlimit;
};

class PropagatorSingleFT : public PropagatorPhaseShift
{
public:
  PropagatorSingleFT(CBuffer2D& wave, const MultisliceParameters& params,
                     const ConstSamplePtr& sample);
  void propagate() final;

private:
  /// Propagate the wave through one slice. This steps combines the
  /// propagation and the application of the phase shift.
  void propagateOneSlice(CBuffer2D& waveIn, CBuffer2D& waveOut,
                         const RBuffer2D& phaseShift);

  CBuffer2D mWaveAfter;

  CBuffer2D mTemp;
  LaplaceFT mLaplaceFT;
  std::unique_ptr<CBandlimit> mBandlimit;
};

/// Propagator partially corrected (called PCMS in [Ming2013])
/// Implemented in real space.
template<class Method>
class PropagatorPartial : public PropagatorPhaseShift
{
public:
  PropagatorPartial(CBuffer2D& wave, const MultisliceParameters& params,
                    const ConstSamplePtr& sample);
  void propagate() final;

private:
  /// Propagates the wave through one slice.
  void propagateOneSlice(CBuffer2D& waveIn, CBuffer2D& waveOut,
                         const RBuffer2D& phaseShift);

  /// Calculates the argument of the exponential function.
  void helper(CBuffer2D& wave, const RBuffer2D& phaseShift);

  CBuffer2D mWaveAfter;
  // we need a temporary buffer with a margin to calculate the discrete
  // Laplacian
  CMarginBuffer2D<Method::margin> mHelp;
  CBuffer2D mOperatorPower;

  std::unique_ptr<CBandlimit> mBandlimit;
};

/// Propagator fully corrected [Chen1997] (called FCMS in [Ming2013])
/// Implementeded in real space.
template<class Method>
class PropagatorFull : public PropagatorPhaseShift
{
public:
  PropagatorFull(CBuffer2D& wave, const MultisliceParameters& params,
                 const ConstSamplePtr& sample);
  void propagate() final;

private:
  /// Propagates the wave through one slice.
  void propagateOneSlice(CBuffer2D& waveIn, CBuffer2D& waveOut,
                         const RBuffer2D& phaseShift);

  /// Calculates the argument of the exponential function.
  void helper(CBuffer2D& wave, const RBuffer2D& phaseShift);

  CBuffer2D mWaveAfter;
  // We need a temporary buffer with on pixel margin to calculate the discrete
  // Laplacian
  CMarginBuffer2D<Method::margin> mHelp;
  // Temporary buffer to store the power of the inner operator
  CBuffer2D mOperatorPower;

  std::unique_ptr<CBandlimit> mBandlimit;
};

/// Propagator fully corrected [Chen1997] (called FCMS in [Ming2013])
/// Implemented in Fourier space.
class PropagatorFullFT : public PropagatorPhaseShift
{
public:
  PropagatorFullFT(CBuffer2D& wave, const MultisliceParameters& params,
                   const ConstSamplePtr& sample);
  void propagate() final;

private:
  /// Propagates the wave through one slice.
  void propagateOneSlice(CBuffer2D& waveIn, CBuffer2D& waveOut,
                         const RBuffer2D& phaseShift);

  /// Calculates the argument of the exponential function.
  void helper(CBuffer2D& wave, const RBuffer2D& phaseShift);

  CBuffer2D mWaveAfter;

  CBuffer2D mHelp;
  CBuffer2D mTemp;
  LaplaceFT mLaplaceFT;
  // Temporary buffer to store the power of the inner operator
  CBuffer2D mOperatorPower;

  const std::unique_ptr<Fft> mFft;
  const std::unique_ptr<Fft> mFftInv;

  std::unique_ptr<CBandlimit> mBandlimit;
};

/// Propagator fully corrected [Chen1997] (called FCMS in [Ming2013]).
/// Simpler series expansion according to [Spiegelberg2015]
/// Implemented in Fourier space.
class PropagatorFullFT2 : public PropagatorPhaseShift
{
public:
  PropagatorFullFT2(CBuffer2D& wave, const MultisliceParameters& params,
                    const ConstSamplePtr& sample);
  void propagate() final;

private:
  /// Propagates the wave through one slice.
  void propagateOneSlice(const RBuffer2D& phaseShift);

  /// Calculates the argument of the exponential function.
  void helper(CBuffer2D& wave, const RBuffer2D& phaseShift);

  CBuffer2D mWaveAfter;

  // Temporary buffer to store the power of the operator
  CBuffer2D mOperatorPower;
  LaplaceFT mLaplaceFT;

  const std::unique_ptr<Fft> mFft;
  const std::unique_ptr<Fft> mFftInv;

  std::unique_ptr<CBandlimit> mBandlimit;

  Complex coefficient(unsigned n);
};

/// Propagator using the Runge-Kutta method and the paraxial approximation
/// of the Schrödinger equation, leaving out the second derivative in the
/// Schrödinger equation.
template<class Pde>
class PropagatorRungeKutta1 : public PropagatorPotential
{
public:
  PropagatorRungeKutta1(CBuffer2D& wave, const MultisliceParameters& params,
                        const ConstSamplePtr& sample);
  void propagate() final;

private:
  // we need a temporary buffers with a margin to safely calculate the
  // discrete Laplace operator
  typename Pde::MarginBuffer mA;
  typename Pde::MarginBuffer mATemp;
  CBuffer2D mAP0;
  CBuffer2D mAPA;
  CBuffer2D mAPB;
  CBuffer2D mAPC;

  std::unique_ptr<CBandlimit> mBandlimitA;

  Pde mPde;
};

/// Propagator using the Runge-Kutta method and the paraxial approximation
/// of the Schrödinger equation, leaving out the second derivative in the
/// Schrödinger equation.
class PropagatorRungeKutta1FT : public PropagatorPotential
{
public:
  PropagatorRungeKutta1FT(CBuffer2D& wave, const MultisliceParameters& params,
                          const ConstSamplePtr& sample);
  void propagate() final;

private:
  void f(CBuffer2D& inA, const RBuffer2D& potential, CBuffer2D& out);
  void diagnostic(size_t slice, Real z, const CBuffer2D& a) const;

  // we need a temporary buffers with a margin to safely calculate the
  // discrete Laplace operator
  CBuffer2D mA;
  CBuffer2D mATemp;
  CBuffer2D mAP0;
  CBuffer2D mAPA;
  CBuffer2D mAPB;
  CBuffer2D mAPC;
  CBuffer2D mTemp;
  LaplaceFT mLaplaceFT;

  const Complex mPotentialFactor;

  std::unique_ptr<CBandlimit> mBandlimitA;
};

/// Propagator using the Runge-Kutta method without approximating the
/// Schrödinger equation. Hence, we solve a second order differential
/// equation.
template<class Pde>
class PropagatorRungeKutta2 : public PropagatorPotential
{
public:
  PropagatorRungeKutta2(CBuffer2D& wave, const MultisliceParameters& params,
                        const ConstSamplePtr& sample);

  void propagate() final;

private:
  // we need temporary buffers with one pixel margin to safely calculate the
  // discrete Laplacian
  typename Pde::MarginBuffer mA;
  typename Pde::MarginBuffer mATemp;
  CBuffer2D mAP0;
  CBuffer2D mAPA;
  CBuffer2D mAPB;
  CBuffer2D mAPC;
  CBuffer2D mB;
  CBuffer2D mBTemp;
  CBuffer2D mBP0;
  CBuffer2D mBPA;
  CBuffer2D mBPB;
  CBuffer2D mBPC;

  std::unique_ptr<CBandlimit> mBandlimitA;
  std::unique_ptr<CBandlimit> mBandlimitB;

  Pde mPde;
};

class PropagatorRungeKutta2FT : public PropagatorPotential
{
public:
  PropagatorRungeKutta2FT(CBuffer2D& wave, const MultisliceParameters& params,
                          const ConstSamplePtr& sample);

  void propagate() final;
  void f(const CBuffer2D& inA, const CBuffer2D& inB,
         const RBuffer2D& potential, CBuffer2D& outBP);

  void diagnostic(size_t slice, Real z,
                  const CBuffer2D& a, const CBuffer2D& b) const;

private:
  // we need temporary buffers with one pixel margin to safely calculate the
  // discrete Laplacian
  CBuffer2D& mA;
  CBuffer2D mATemp;
  CBuffer2D mAP0;
  CBuffer2D mAPA;
  CBuffer2D mAPB;
  CBuffer2D mAPC;
  CBuffer2D mB;
  CBuffer2D mBP0;
  CBuffer2D mBPA;
  CBuffer2D mBPB;
  CBuffer2D mBPC;
  CBuffer2D mTemp;
  LaplaceFT mLaplaceFT;

  Real mWaveNum;
  Real mPotentialFactor;

  std::unique_ptr<CBandlimit> mBandlimitA;
  std::unique_ptr<CBandlimit> mBandlimitB;
};

/// Base class for all partial differential equation (PDE) that can be solved
/// with the methods declared in this header.
class PdeBase
{
public:
  /// Triplets are needed by Eigen to construct a sparse matrix. A triplet
  /// consists of the row, the column and the value of a single matrix element.
  typedef Eigen::Triplet<Complex> Triplet;
  typedef std::vector<Triplet> Triplets;

  typedef Eigen::Matrix<Complex, Eigen::Dynamic, 1> Vector;
  typedef Eigen::SparseMatrix<Complex, 0, int64_t> SparseMatrix;

  /// Constructor.
  PdeBase(const MultisliceParameters& params);
  /// This function is needed for the linearization of our two-dimensional
  /// buffers.
  /// @par When we use linear algebra to calculate the propagation, we have to
  /// convert out two-dimensional buffers into one-dimensional vectors. This
  /// functions return the index of an element in the array given its @p column
  /// and @p row.
  /// @par Furthermore, it supports the notion of block matrices/vectors.
  /// These are needed for the PDEs of the second order (@ref Pde2 and @ref
  /// Pde3) and the higher order algorithms (@ref PropagatorLobatto1,
  /// @ref PropagatorLobatto2).
  int64_t index(int64_t column, int64_t row, int64_t block = 0) const;
  /// This function does the same as @ref index, but it wraps around @p column
  /// and @p row. This needed for the calculation of the discrete laplace o
  /// perator.
  int64_t indexS(int64_t column, int64_t row, int64_t block = 0) const;

  void bufferToEigenVector(const CBuffer2D& buffer, Vector& vector,
                           int64_t part = 0);
  void eigenVectorToBuffer(const Vector& vector, CBuffer2D& buffer,
                           int64_t part = 0);

protected:
  // Number of columns
  int64_t mCols;
  // Number of rows
  int64_t mRows;
  // Number of pixels
  int64_t mNumPixels;
};

/// Class to represent the PDE
/// @f{eqnarray*}{
///   \partial_z\,A &=& \frac{\text{i}}{2\,k}\Delta_\perp\,A -
///                     \frac{\text{i}\,m\,\gamma}{k\,\hbar^2}\,V\,A
/// @f}
template<class Method>
class Pde1 : public PdeBase
{
public:
  typedef CMarginBuffer2D<Method::margin> MarginBuffer;

  Pde1(const MultisliceParameters& params);

  /// Calculates the derivative @f$\partial_z\,A@f$ and stores it in @p out.
  void f(MarginBuffer& inA, const RBuffer2D& potential,
         CBuffer2D& out) const;
  /// Calculates
  /// @f{eqnarray*}{
  ///   A + \text{weight}\times f(A,\,V) \to \text{out}
  /// @f}
  void g(MarginBuffer& inA, const RBuffer2D& potential, Real weight,
         CBuffer2D& out) const;
  void diagnostic(size_t slice, Real z, const MarginBuffer& a) const;
  /// Creates the triplets that are needed by Eigen to construct the matrix
  /// describing the PDE. This matrix is multiplied by @p weight and an identity
  /// matrix is added, if @p leftBlock equals @p rightBlock. The block variables
  /// are used to construct the matrix for Lobatto algorithms.
  void fillTriplets(Triplets& triplets, int64_t rowBlock, int64_t columnBlock,
                    Real weight, const RBuffer2D& potential) const;

private:
  Complex mLaplaceFactor;
  Complex mPotentialFactor;
};

/// Class to represent the PDE
/// @f{eqnarray*}{
///   \partial_z\,A &=& B \\
///   \partial_z\,B &=& -2\,\text{i}\,k\,B-\Delta_\perp\,A + \frac{2\,m\,\gamma}
///                     {\hbar^2}\,V\,A = f(A,\,B)
/// @f}
template<class Method>
class Pde2 : public PdeBase
{
public:
  typedef CMarginBuffer2D<Method::margin> MarginBuffer;

  Pde2(const MultisliceParameters& params);

  /// Initializes the buffer @p b to an appropriate initial value. For this
  /// PDE @p b is initialized to zero.
  void initB(const MarginBuffer& a, CBuffer2D& b,
             const std::unique_ptr<Potential>& potential);
  /// Calculates the derivative @f$\partial_z\,B=f(A,\,B,\,V)@f$ and stores it
  /// in @p outBP.
  void f(MarginBuffer& inA, const CBuffer2D& inB,
         const RBuffer2D& potential, CBuffer2D& outBP) const;
  /// Calculates
  /// @f{eqnarray*}{
  ///   B + \text{weight}\times f(A,\,B,\,V) \to \text{out}
  /// @f}
  void g(MarginBuffer& inA, const CBuffer2D& inB,
         const RBuffer2D& potential, Real weight, CBuffer2D& out) const;
  void diagnostic(size_t slice, Real z,
                  const MarginBuffer& a, const CBuffer2D& b) const;
  /// @copydoc Pde1::fillTriplets()
  void fillTriplets(Triplets& triplets, int64_t rowBlock, int64_t columnBlock,
                    Real weight, const RBuffer2D& potential) const;

private:
//void fStart(MarginBuffer& in, const RBuffer2D& potential, MarginBuffer& out);

  Real mWaveNum;
  Real mLaplaceFactor;
  Real mPotentialFactor;
};

/// Class to represent the PDE
/// @f{eqnarray*}{
///   \partial_z\,A &=& B \\
///   \partial_z\,B &=& -\Delta_\perp\,A + \frac{2\,m\,\gamma}{\hbar^2}\,V\,A
///                     -k^2\,A
/// @f}
template<class Method>
class Pde3 : public PdeBase
{
public:
  typedef CMarginBuffer2D<Method::margin> MarginBuffer;

  Pde3(const MultisliceParameters& params);

  /// Initializes the buffer @p b to an appropriate initial value. For this
  /// PDE @p b is initialized to @f$\text{i}\,k@f$.
  void initB(const MarginBuffer& a, CBuffer2D& b,
             const std::unique_ptr<Potential>& potential);
  /// Calculates the derivative @f$\partial_z\,B=f(A,\,B,\,V)@f$ and stores it
  /// in @p outBP.
  void f(MarginBuffer& inA, const CBuffer2D& inB,
         const RBuffer2D& potential, CBuffer2D& outBP) const;
  /// @copydoc Pde2::g()
  void g(MarginBuffer& inA, const CBuffer2D& inB,
         const RBuffer2D& potential, Real weight, CBuffer2D& out) const;
  void diagnostic(size_t slice, Real z,
                  const MarginBuffer& a, const CBuffer2D& b) const;
  /// @copydoc Pde1::fillTriplets()
  void fillTriplets(Triplets& triplets, int64_t rowBlock, int64_t columnBlock,
                    Real weight, const RBuffer2D& potential) const;

private:
  Real mWaveNum;
  Real mEnergyFactor;
  Real mLaplaceFactor;
  Real mPotentialFactor;
  Real mInitialBAbsSqr;
};

/// Base class for solving partial differential equation of the first order
/// using implicit schemes.
template<class Pde>
class PropagatorImplicit1 : public PropagatorPotential
{
public:
  PropagatorImplicit1(CBuffer2D& wave, const MultisliceParameters& params,
                      const ConstSamplePtr& sample);

protected:
  typename Pde::MarginBuffer mA;
  CBuffer2D mTemp;

  std::unique_ptr<CBandlimit> mBandlimitWave;

  Pde mPde;
};

template<class Pde>
class PropagatorTrapezoidal1 : public PropagatorImplicit1<Pde>
{
public:
  using PropagatorImplicit1<Pde>::PropagatorImplicit1;
  void propagate() final;
};

template<class Pde>
class PropagatorLobatto1 : public PropagatorImplicit1<Pde>
{
public:
  using PropagatorImplicit1<Pde>::PropagatorImplicit1;
  void propagate() final;
};

template<class Pde>
class PropagatorImplicit2 : public PropagatorPotential
{
public:
  PropagatorImplicit2(CBuffer2D& wave, const MultisliceParameters& params,
                      const ConstSamplePtr& sample);

protected:
  typename Pde::MarginBuffer mA;
  CBuffer2D mB;
  CBuffer2D mTemp;

  std::unique_ptr<CBandlimit> mBandlimitWave;
  std::unique_ptr<CBandlimit> mBandlimitB;

  Pde mPde;
};

template<class Pde>
class PropagatorTrapezoidal2 : public PropagatorImplicit2<Pde>
{
public:
  using PropagatorImplicit2<Pde>::PropagatorImplicit2;
  void propagate() final;
};

template<class Pde>
class PropagatorLobatto2 : public PropagatorImplicit2<Pde>
{
public:
  using PropagatorImplicit2<Pde>::PropagatorImplicit2;
  void propagate() final;
};

} // namespace Aurora

#endif
