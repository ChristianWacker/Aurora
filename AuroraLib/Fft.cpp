//--- Aurora/AuroraLib/Fft.cpp -------------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Fft.hpp"

#if (1 == AURORA_USE_FFTW)
#include "AuroraLib/FftFftw.hpp"
#endif

#include "AuroraLib/FftKiss.hpp"

namespace Aurora
{

/*static*/ std::unique_ptr<Fft> Fft::create(const CBuffer2D& bufferIn,
                                            CBuffer2D& bufferOut,
                                            Direction direction, bool once)
{
#if (1 == AURORA_USE_FFTW)
  if (once)
  {
    return std::make_unique<FftFftw<FFTW_ESTIMATE>>(bufferIn, bufferOut,
                                                    direction);
  }
  else
  {
    return std::make_unique<FftFftw<FFTW_MEASURE>>(bufferIn, bufferOut,
                                                   direction);
  }
#else
  (void)once;
  return std::make_unique<FftKiss>(bufferIn, bufferOut, direction);
#endif
}

} // namespace Aurora
