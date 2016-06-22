//--- Aurora/AuroraLib/Fft.hpp -------------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_FFT_HPP
#define AURORA_AURORA_LIB_FFT_HPP

#include "AuroraLib/Factory.hpp"

namespace Aurora
{

enum class Direction
{
  forward, backward
};

/// Abstract base class representing a Fast Fourier Transformation (FFT) in
/// arbitrary dimensions
/// @remarks
///  The FFT is just an efficient way to implement a discrete Fourier
///  transformation (DFT). There are many different definitions of the DFT,
///  that only differ in sign and/or normalization. Therefore, one has to
///  specify the used convention. The discrete Fourier transformation (forward
///  transformation) shall be calculated as
///    @f[
///           y_k = \sum_{j=0}^{n-1} x_j \exp(+2\,\pi\,\mathrm{i}\,j\,k / n)
///    @f]
///  with @f$k\in \{0,\,1,\,\ldots,\,n-1\}@f$ and @f$n@f$ specifying the
///  number of sampling points. @f$n@f$ is sometimes called the dimension of
///  the DFT. But we will refrain from this to avoid confusion with the other
///  meanings of the word "dimension". The inverse discrete Fourier
///  transformation (backward transformation) only differs in the sign of the
///  exponent
///    @f[
///           y_k = \sum_{j=0}^{n-1} x_j \exp(-2\,\pi\,\mathrm{i}\,j\,k / n)
///    @f]
///  Because of the chosen normalization a FFT followed by an inverse FFT
///  results in a scaling of the original values with the number of elements,
///  i.e. @f$n@f$.
/// @remarks
///  The multi-dimensional FFT is the repeated application of the one
///  dimensional DFT to each dimension.
/// @remarks
///  The FFT is always out-of place, i.e. an input and output buffer are
///  needed. Currently only complex-to-complex transformations are supported.
class AURORA_API Fft
{
public:
  /// Release the data structures that describe the FFT.
  virtual ~Fft() {}

  /// Executes the actual discrete Fourier transformation
  virtual void transform() = 0;

  /// Returns the input buffer.
  const CBuffer2D& inBuffer() const
  {
    return mInBuffer;
  }

  /// Returns the output buffer.
  CBuffer2D& outBuffer() const
  {
    return mOutBuffer;
  }

  Direction direction() const
  {
    return mDirection;
  }

  static std::unique_ptr<Fft> create(const CBuffer2D& bufferIn,
                                     CBuffer2D& bufferOut, Direction direction,
                                     bool once);

protected:
  /// Creates a Fft.
  /// @param inBuffer
  ///  Pointer to the input buffer for the FFT.
  /// @param outBuffer
  ///  Pointer to the output buffer for the FFT.
  /// @param direction
  ///  The direction of the Fft see @ref Direction.
  Fft(const CBuffer2D& inBuffer, CBuffer2D& outBuffer, Direction direction) :
    mInBuffer(inBuffer), mOutBuffer(outBuffer), mDirection(direction) { }

  const CBuffer2D& mInBuffer;
  CBuffer2D& mOutBuffer;
  const Direction mDirection;
};

} // namespace Aurora

#endif
