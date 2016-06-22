//--- Aurora/AuroraLib/FftFftw.hpp ---------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_FFT_FFTW_HPP
#define AURORA_FFT_FFTW_HPP

#include "AuroraLib/Fft.hpp"

#include "AuroraLib/CBuffer2D.hpp"

#include <fftw3.h>
#include <mutex>

#if (1 == AURORA_SINGLE_PRECISION)
# define fftw_complex fftwf_complex
# define fftw_plan fftwf_plan
# define fftw_destroy_plan fftwf_destroy_plan
# define fftw_execute fftwf_execute
# define fftw_plan_dft_1d fftwf_plan_dft_1d
# define fftw_plan_dft_2d fftwf_plan_dft_2d
# define fftw_init_threads fftwf_init_threads
# define fftw_plan_with_nthreads fftwf_plan_with_nthreads
# define fftw_cleanup_threads fftwf_cleanup_threads
# define fftw_plan_guru64_dft fftwf_plan_guru64_dft
#endif

namespace Aurora
{

namespace Fftw
{

AURORA_API void init();
/// Mutex to serialize access to the FFTW planner routines
AURORA_API std::mutex& getMutex();

} // namespace fftw

/// Implementation of @ref Fft using the FFTW library.
/// @remarks
///  The FFTW library (http://www.fftw.org/) is among the fastest libraries
///  performing a FFT on the CPU, though it has not been optimized by a
///  hardware vendor. This class should not directly be instantiated. Use
///  instead @ref Fft::create().
/// @tparam fftwFlags
///  These flags will be used to create a plan. Their main purpose is to
///  choose the planning mode.
template<unsigned fftwFlags>
class FftFftw : public Fft
{
public:
  /// Creates a FFT for the transformation of iBuffer to outBuffer. It calls the
  /// FFTW library to generate a plan that is optimized for the given input and
  /// output buffers.
  /// @param inBuffer
  ///  Pointer to the @ref CBuffer2D that will contain the input values.
  /// @param outBuffer
  ///  Pointer to the @ref CBuffer2D that will receive the output values.
  /// @param direction
  ///  Direction of the Fourier transformation.
  FftFftw(const CBuffer2D& inBuffer, CBuffer2D& outBuffer,
          Direction direction) :
    Fft(inBuffer, outBuffer, direction)
  {
    Fftw::init();

    // We need a mutex, as all accesses to the FFTW planning routines must be
    // serialized
    std::unique_lock<std::mutex> lock(Fftw::getMutex());

    // This part is a little tricky. The input data are constant, but FFTW needs
    // to have write access to both buffers for the planning. Therefore, we will
    // cast away the const qualifier, copy the data of both buffers in temporary
    // buffers, do the planning and restore the buffer. Afterwards there will be
    // no further write accesses to the input buffer.
    auto& inBufferNC = const_cast<CBuffer2D&>(inBuffer);

    const int cols = static_cast<int>(inBufferNC.cols());
    const int rows = static_cast<int>(inBufferNC.rows());

    // Secure the input and output
    auto inBackup  = inBufferNC.clone();
    auto outBackup = outBuffer.clone();

    // convert to FFTW format (the sign convention of FFTW is opposite)
    int sign = (Direction::forward == direction) ? FFTW_BACKWARD : FFTW_FORWARD;
    fftw_complex* inData = reinterpret_cast<fftw_complex*>(inBufferNC.data());
    fftw_complex* outData = reinterpret_cast<fftw_complex*>(outBuffer.data());

    // To calculate a FFT with FFTW one has to create a plan that adapts to the
    // current hardware. The one-dimensional and two-dimensional FFT have
    // different planning routines.
    if (1 == rows)
    {
      mPlan = fftw_plan_dft_1d(cols, inData, outData, sign, fftwFlags);
    }
    else
    {
      // We need to swap the rows and columns as we use column-major order while
      // FFTW uses the row-major order
      mPlan = fftw_plan_dft_2d(rows, cols, inData, outData, sign, fftwFlags);
    }

    // Restore the buffers.
    inBufferNC.assign(inBackup);
    outBuffer.assign(outBackup);
  }

  /// @copydoc Fft::~Fft
  ~FftFftw() override
  {
    std::unique_lock<std::mutex> lock(Fftw::getMutex());
    fftw_destroy_plan(mPlan);
  }

  /// @copydoc Fft::transform
  void transform() override
  {
    fftw_execute(mPlan);
  }

private:
  void initFftw();

  /// Plan describing FFT. It contains information about the buffers, the
  /// direction of the FFT, data layout and so on.
  fftw_plan mPlan;
};

}

#endif
