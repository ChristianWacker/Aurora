
//--- Aurora/AuroraLib/FftFftw.cpp ---------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/FftFftw.hpp"

#include "AuroraLib/CommandLine.hpp"

#include <iostream>

#ifdef _OPENMP
# include <omp.h>
#endif

#include <thread>

namespace Aurora
{

namespace Fftw
{

std::mutex mutex;

class AURORA_API Init
{
public:
  Init();
  ~Init();

private:
  size_t mNumThreads;
};

Init::Init() :
  #ifdef _OPENMP
  mNumThreads(omp_get_max_threads())
  #else
  mNumThreads(std::thread::hardware_concurrency())
  #endif
{
  assert(1 <= mNumThreads);

  if (1 == mNumThreads)
  {
    if (verbosity >= 1)
      std::cout << "Using FFTW without threads." << std::endl;
    return;
  }

  if (0 == fftw_init_threads())
    AURORA_THROW(EUnknown, "fftw_init_threads");

  // Use all the cores for the FFT
  fftw_plan_with_nthreads(static_cast<int>(mNumThreads));

  if (verbosity >= 1)
  {
    std::cout << "Initialized FFTW with " << mNumThreads << " threads."
              << std::endl;
  }
}

Init::~Init()
{
  if (mNumThreads > 1)
    fftw_cleanup_threads();

  if (verbosity >= 1)
    std::cout << "FFTW has been cleaned up.\n";
}

void init()
{
  static Init i;
}

std::mutex& getMutex()
{
  return mutex;
}

} // namespace Fftw

} // namespace Aurora
